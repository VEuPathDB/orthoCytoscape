#!/usr/bin/perl
use warnings;
use strict;

use JSON;
use Data::Dumper;


# OG6_112934 is good, not too big, but has EC and domains

my $originalJsonFile = $ARGV[0];
my $outJsonFile = $ARGV[1];

my $maxSlices = 16; # max number of node pie slices for ec numbers and pfam domains

my $originalJsonText = getJsonTextFromFile($originalJsonFile);
my $originalHash = decode_json($originalJsonText);

my $cytoscapeHash = makeCytoscapeFormat($originalHash,$maxSlices);

open(OUT,">",$outJsonFile) || die "Cannot open $outJsonFile for reading\n";
print OUT encode_json($cytoscapeHash);
close OUT;


#  foreach my $gene (keys %{$originalHash->{group}->{genes}}) {
#      print Dumper $originalHash->{group}->{genes}->{$gene};
#      exit;
#}
#   print Dumper $originalHash->{taxonCounts};
#     exit;
# }



exit;












sub getJsonTextFromFile {
    my ($file)=@_;
    local $/;
    open(IN,$file) || die "Cannot open $file\n";
    my $jsonText = <IN>;
    close IN;
    return $jsonText;
}

sub makeCytoscapeFormat {
    my ($oldFormat,$maxSlices) = @_;
    
    my $cytoscapeFormat = [];
    addNodes($cytoscapeFormat,$oldFormat,$maxSlices);
    addEdges($cytoscapeFormat,$oldFormat);

    return $cytoscapeFormat;
}

sub addNodes {
    my ($cytoscapeFormat,$oldFormat,$maxSlices) = @_;

    my $groupEcNums = getGroupEcNums($oldFormat);
    my $groupPfams = getGroupPfams($oldFormat);
    my $abbrevToGroupColor = getGroupColor($oldFormat);

    foreach my $node (keys %{$oldFormat->{nodes}}) {

	my $nodeObject;
	$nodeObject->{classes} = ["taxa"];
	
	my ($x,$y) = getXandY($oldFormat->{nodes}->{$node});
	$nodeObject->{position}->{x} = $x;
	$nodeObject->{position}->{y} = $y;

	my $id = $oldFormat->{nodes}->{$node}->{id};
	$nodeObject->{data}->{id} = $id;

	my $speciesColor =  $oldFormat->{group}->{genes}->{$id}->{taxon}->{color};
	$nodeObject->{data}->{speciesColor} = $speciesColor;

	my %nodeEcHash = map { $_ => 1 } @{$oldFormat->{group}->{genes}->{$id}->{ecNumbers}};
	createPie($nodeObject,$groupEcNums,$maxSlices,"ec");
	finalizePie($nodeObject,\%nodeEcHash,$groupEcNums,"ec");

	my %nodePfamHash = map { $_ => 1 } keys %{$oldFormat->{group}->{genes}->{$id}->{pfamDomains}};
        createPie($nodeObject,$groupPfams,$maxSlices,"pfam");
        finalizePie($nodeObject,\%nodePfamHash,$groupPfams,"pfam");

	my $abbrev = $oldFormat->{group}->{genes}->{$id}->{taxon}->{abbrev};
	my $taxonGroupColor = $abbrevToGroupColor->{$abbrev};
	$nodeObject->{data}->{groupColor} = $taxonGroupColor;

	push @{$cytoscapeFormat}, $nodeObject;

    }
}


sub addEdges {
    my ($cytoscapeFormat,$oldFormat) = @_;

    my $edgeNum=1;
    foreach my $edge (keys %{$oldFormat->{edges}}) {

	my $edgeObject;

	$edgeObject->{data}->{id} = "edge".$edgeNum;
	$edgeObject->{data}->{source} = $oldFormat->{edges}->{$edge}->{queryId};
	$edgeObject->{data}->{target} = $oldFormat->{edges}->{$edge}->{subjectId};
	$edgeObject->{data}->{eValue} = $oldFormat->{edges}->{$edge}->{E};

	my $edgeTypeSingleLetter = $oldFormat->{edges}->{$edge}->{T};
	$edgeObject->{data}->{type} = getEdgeType($edgeTypeSingleLetter);
	
	$edgeNum++;

	push @{$cytoscapeFormat}, $edgeObject;

    }
}


sub getGroupEcNums {
    my ($oldFormat) = @_;

    #create hash ec->count, so that can sort by decreasing occurence of EC number
    my $ecToCount;
    foreach my $ecNum (keys %{$oldFormat->{group}->{ecNumbers}}) {
	$ecToCount->{$ecNum} = $oldFormat->{group}->{ecNumbers}->{$ecNum}->{count};
    }

    my $ecNums = [];
    foreach my $ecNum (sort { $ecToCount->{$b} <=> $ecToCount->{$a} } keys %{$ecToCount}) {
	push @{$ecNums}, {
			   accession =>  $ecNum,
			   color     =>  $oldFormat->{group}->{ecNumbers}->{$ecNum}->{color}
	                 };
    }
    return $ecNums;
}


sub getGroupPfams {
    my ($oldFormat) = @_;

    #create hash pfam->count, so that can sort by decreasing occurence of pfam domain
    my $pfamToCount;
    foreach my $pfam (keys %{$oldFormat->{group}->{pfamDomains}}) {
	$pfamToCount->{$pfam} = $oldFormat->{group}->{pfamDomains}->{$pfam}->{count};
    }
    
    my $pfamDomains = [];
    foreach my $pfam (sort { $pfamToCount->{$b} <=> $pfamToCount->{$a} } keys %{$pfamToCount}) {
	push @{$pfamDomains}, { 
				accession =>  $pfam,
				color     =>  $oldFormat->{group}->{pfamDomains}->{$pfam}->{color}
	                      };
    }
    return $pfamDomains;
}


sub getXandY {
    my ($node) = @_;
	
    my $x = $node->{x};
    $x =~ s/"//g;
    $x = int($x);
    my $y =  $node->{y};
    $y =~ s/"//g;
    $y = int($y);

    return ($x,$y);
}


	

sub createPie {
    my ($node,$groupSliceInfo,$maxSlices,$type) = @_;

    my $numInputSlices = scalar @{$groupSliceInfo};
    my $numSlices = ($numInputSlices,$maxSlices)[$numInputSlices>$maxSlices];
    my $pct = 0;
    $pct = sprintf("%.0f",100/$numSlices) if ($numSlices > 0);

    # 16 slices is the cytoscape limit; some of these 16 will be empty depending on number of slices needed or maxSlices
    for (my $i=0; $i<16; $i++) {
	my $sliceNum=$i+1;
	my $pctField = $type.$sliceNum."pct";
	my $colorField = $type.$sliceNum."color";
	$pct = 0 if ($i >= $numSlices);
	my $color = "white";
	$color = $groupSliceInfo->[$i]->{color} if (exists $groupSliceInfo->[$i]);
	$node->{data}->{$colorField} = $color;
	$node->{data}->{$pctField} = $pct;
    }
}


sub finalizePie {
    my ($node,$nodeInfo,$groupSliceInfo,$type) = @_;

    return if (scalar @{$groupSliceInfo} == 0);

    for (my $i=0; $i<16; $i++) {
	if (exists $groupSliceInfo->[$i] && ! exists $nodeInfo->{$groupSliceInfo->[$i]->{accession}} ) {
	    my $sliceNum=$i+1;
	    my $colorField = $type.$sliceNum."color";
	    $node->{data}->{$colorField} = "white";
	}
    }
}


sub getGroupColor {
    my ($oldFormat) = @_;

    my $abbrevToGroupColor={};

    foreach my $taxon (keys %{$oldFormat->{taxons}}) {
	if (exists $oldFormat->{taxons}->{$taxon}->{groupColor}) {
	    my $groupColor = $oldFormat->{taxons}->{$taxon}->{groupColor};
	    getSpecies($oldFormat->{taxons}->{$taxon}->{children},$abbrevToGroupColor,$groupColor);	       
	}
    }
    return $abbrevToGroupColor;
}


sub getSpecies {
    my ($children,$abbrevToGroupColor,$groupColor) = @_;

    foreach my $child (keys %{$children}) {
	my $numChildren = keys %{$children->{$child}->{children}};
	if ($numChildren == 0) {
	    $abbrevToGroupColor->{$children->{$child}->{abbrev}} = $groupColor;
	} else {
	    getSpecies($children->{$child}->{children},$abbrevToGroupColor,$groupColor);
	}
    }
}


sub getEdgeType {
    my ($edgeTypeSingleLetter) = @_;

    my %edgeType = (
	    "O" => "Ortholog",
	    "C" => "Coortholog",
	    "P" => "Inparalog",
	    "L" => "PeripheralCore",
	    "M" => "PeripheralPeripheral",
	    "N" => "Other"
                   );

    die "edge type not known" if (! exists  $edgeType{$edgeTypeSingleLetter});
    
    return $edgeType{$edgeTypeSingleLetter};

}

