#!/usr/bin/perl -w
use strict;
use FileHandle;
use Getopt::Long;

# Morag Lewis
# morag.lewis@kcl.ac.uk

#Prediction of Potential Causative Regulatory Networks
#Perl script to find links between a master regulator and misregulated genes (eg from a microarray)

##############################################################################
## Copyright 2017 Morag Lewis                                               ##
##                                                                          ##
## Licensed under the Apache License, Version 2.0 (the "License");          ##
## you may not use this file except in compliance with the License.         ##
## You may obtain a copy of the License at                                  ##
##                                                                          ##
##   http://www.apache.org/licenses/LICENSE-2.0                             ##
##                                                                          ##
## Unless required by applicable law or agreed to in writing, software      ##
## distributed under the License is distributed on an "AS IS" BASIS,        ##
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. ##
## See the License for the specific language governing permissions and      ##
## limitations under the License.                                           ##
##                                                                          ##
##############################################################################

# This script takes as input:
# 1. A network of known interactions. This can be generated however appropriate, but the input format should be a sif file, tab-separated, one interaction per line:
#   eg. ENSMUSG00000000078      upregulates     ENSMUSG00000020377	source
#   The only interactions the script recognises are "upregulates" and "downregulates", but it can handle multiple gene combinations, eg microarray data from a double
#   gene knockout. In that case, the IDs should be separated by ;
# 2. A list of misregulated gene IDs and the misregulation value eg a fold change, one per line, tab-separated.
# 3. A list of genes at which the algorithm stops, with an indication of whether those genes are up- or downregulated (ie for a knockout, the value would be -1). 
#    These genes will be not be checked for upstream regulators. Use one gene per line, with the value separated by a tab (same format as the misregulated genes).
#    It's fine if the stop genes are also in the misregulation list, but their expression values should be in concordance - ie if you have a knockout, but the misregulation
#    is positive (an upregulation), something about your data needs checking.
# Source is simply where the data came from (eg ArrayExpress, or a PubMed ID)
#
# Optional inputs:
# --prefix=TMP (default) - this prefixes all the output files with the specified prefix
# --short - if this is included, the output network will consist of only the shortest paths between each misregulated gene and the stoplist. Default is off.
# --debug - outputs what the script is doing. Off by default.
# --targets - hardwires the direction of regulation of immediate downstream targets of stoplist in appropriate direction, rather than allowing their misregulation to be determined by the algorithm taking other genes into account. Default is off. Does not override known misregulation.

# The output is an internally consistent network linking one or more regulators with their known misregulated genes, as experimentally observed, in sif format.
# A list of weightings is also output, containing both the known misregulation values and predictions for genes whose misregulation is not known. For example:
# X upregulates Y.
# X is known to be downregulated, with a fold change of -4.
# Y is therefore also predicted to be downregulated, and will be given a predicted weighting of -1.
# The file will contain the following lines:
# X	-4	known
# Y	-1	predicted
# This can be used with Cytoscape to colour known and predicted nodes as required.
# A list of interactions with sources supporting each - this can be used with Cytoscape to add information to the edges
# The final output is a list of unused IDs, which are the genes from the microarray which couldn't be included in the network.


($#ARGV > -1) or die("Usage: popcorn.pl [ options ] known interactions sif file, microarray file, stop list\n");

my $debug=''; #optional variable for outputting details of what script is doing
my $short=''; #optional variable for shortest path, default value false
my $prefix='TMP'; #optional prefix for output files, default value TMP
my $targets=''; #optional variable for hardwiring weights of immediate targets of genes on stoplist
GetOptions ('short' => \$short, 'debug' => \$debug, 'prefix=s' => \$prefix, 'targets' => \$targets );

my ($sif, $microarray, $stopin);

($sif, $microarray, $stopin) = @ARGV;

print STDERR "Compiling network from input files $sif and $microarray, with stoplist $stopin and prefix $prefix, with shortest path set to $short, and targets set to $targets\n";

###############
## Variables ##
###############

my %weights = ( ); #this is for holding the predicted misregulation values
my %foldchange = ( ); #this is for storing the fold changes so they can be kept for final output
my %times = ( ); #this is for monitoring how many times something has been weighed
my %relevant = (); #this is for tracking which genes are relevant for connecting misregulated genes to stop list
my %source = (); #this holds the source of the data for each regulation
my $stoplist = ""; #list of master regulators (from stopin)
my $unused = "Unused IDS:\n"; #string to hold microarray genes which aren't in the network to start with
my %direct; #This is for holding direct targets of stoplist genes, only used when targets is set

#########################################################
## Reading in the gene(s) at which the algorithm stops ##
#########################################################

print STDERR "Reading in starting gene(s) from $stopin\n";
if ($debug) { print "Reading in starting gene(s) from $stopin\n"; }

my $fhsl = new FileHandle('< ' . $stopin);
while(<$fhsl>) {
	my ( $newstop, $stopweight) = split(/\t/);
	chomp($stopweight);
	$newstop =~ s/\s+//g;
	$stopweight =~ s/\s+//g;
	$stopweight =~ s/\r//g;

        $stoplist = $stoplist.$newstop."\t"; #add each entry to the stoplist, which is just one text string for easy searching, separated by tabs
	$foldchange{ $newstop } = $stopweight; #set foldchange for this gene
	$relevant{ $newstop } = 1; #set the relevant hash
	if ($stopweight > 0) {     #and set the weights hash value accordingly
		$weights{ $newstop } = 1;
	} elsif ( $stopweight < 0 ) {
		$weights{ $newstop } = -1;
	}
}
$fhsl->close();

if ($debug) { print "Stoplist is $stoplist\n"; }

########################
## Reading in network ##
########################

print STDERR "Reading in network from $sif\n";
if ($debug) { print "Reading in network from $sif\n"; }

#hashes for holding regulation
my %upby = ( );
my %downby = ( );

my $locked = ""; #string holding genes whose regulation is not to be changed (ie it's already known)

#Read in the sif file with all known links.
my $fh1 = new FileHandle('< ' . $sif);
while(<$fh1>) {
	my($molecule1, $activity, $molecule2, $ref) = split(/\t/);
	chomp($ref);
	$molecule1 =~ s/\s+//g;
	$activity =~ s/\s+//g;
	$molecule2 =~ s/\s+//g;
	$ref =~ s/\r//g;

	if ($debug) { print "Molecule 1 is $molecule1, molecule 2 is $molecule2, activity is $activity and reference is $ref\n"; }

        if ( $activity =~ m/upregulates/ ) {
                if ( exists( $upby{ $molecule2 } ) ) { #if there's already an entry for molecule 2
                        if ( $upby{ $molecule2 } !~/$molecule1/ ) { #and if molecule 1 isn't there already
                                $upby{ $molecule2 } = $upby{ $molecule2 }."\t".$molecule1; #add molecule 1, the regulator
                        }
                }
                else {
                        $upby{ $molecule2 } = $molecule1; #else, make the hash entry
                }
        }
        elsif ( $activity =~ m/downregulates/ ) {
                if ( exists( $downby{ $molecule2 } ) ) { #if there's already an entry for molecule 2
                        if ( $downby{ $molecule2 } !~/$molecule1/ ) { #and if molecule 1 isn't there already
                                $downby{ $molecule2 } = $downby{ $molecule2 }."\t".$molecule1; #add molecule 1, the regulator
                        }
                }
                else {
                        $downby{ $molecule2 } = $molecule1; #else, make the hash entry
                }
        }

	#if the targets option is set, and molecule 1 is on the stoplist, add molecule1-molecule 2 to the square hash of direct targets and make the entry the activity
	if ( ( $targets ) && ( $stoplist =~ /$molecule1/ ) ) {
		if ( $debug ) { print "Direct target of $molecule1 found: $molecule2\n"; }

		#if $molecule2 is already in the direct target hash
		if ( exists ($direct{ $molecule1 }{ $molecule2} ) ) {
			#if the recorded regulation is different to the current regulation, and is not "both", set it t"both"
			#it will no longer be considered as a direct target, but will be determined by the context of the network like all the other genes
			#in this case, the only way to force the target is to make it part of the stoplist
			#if the regulation is the same, or it is already "both", don't do anything
			if ( ( $direct{ $molecule1 }{ $molecule2 } !~ /$activity/ ) && ( $direct{ $molecule1 }{ $molecule2 } !~ m/both/ ) ) {
				$direct{ $molecule1 }{ $molecule2 } = "both";
			}
		} else {
			$direct{ $molecule1 }{ $molecule2 } = $activity;
		}
	}

	#deal with double knockouts (and they should only occur in $molecule1) - add their constituent parts to weights and times hashes
	if ( $molecule1 =~ /;/ ) {
		my ( @mols ) = split(/;/, $molecule1);
		my $jim;
		for $jim ( 0 .. $#mols ) {
			if (!(exists($weights{ $mols[ $jim ] }) ) ) { $weights{ $mols[ $jim ] } = 0; }
			$times{ $mols[ $jim ] } = 0;
		}
	}

	#set the weight hash, if weight does not already exist
	if ( !(exists($weights{ $molecule1 })) ) { $weights{ $molecule1 } = 0; }
	if ( !(exists($weights{ $molecule2 })) ) { $weights{ $molecule2 } = 0; }
	$times{ $molecule1 } = 0; #set the number of times tested to 0
	$times{ $molecule2 } = 0;


	#Finally, add the source of the data to the source hash
	my $sourcekey = $molecule1." (".$activity.") ".$molecule2;
	if( exists( $source{ $sourcekey } ) ) { #if there's already an entry for this regulation
		if ($source{ $sourcekey } =~ /$ref/) { #check to see if this particular source has already been noted (there will be duplicate entries)
			$source{ $sourcekey } = $source{ $sourcekey }.";".$ref; #if not, add it
		}
	} else { #if there isn't an entry, make one
		$source{ $sourcekey } = $ref;
	}
}
$fh1->close();

##########################################################
## Reading in the misregulated genes eg microarray data ##
##########################################################

print STDERR "Reading in microarray data from $microarray\n";
if ($debug) { print "Reading in microarray data from $microarray\n"; }

my @marray = []; #array of strings of genes to be put into network

my $fh3 = new FileHandle('< ' . $microarray);
while(<$fh3>) {
	my ( $gene, $exp ) = split/\s+/; #get the gene and the expression value
	$gene =~ s/\s+//g;
	$exp =~ s/\s+//g;
	$exp =~ s/\r//g;
	if ( exists( $weights{ $gene } )) { #only consider the genes which are present in the network
		#set the weight hash
		if ($exp > 0) {     
			$weights{ $gene } = 1;
		} elsif ( $exp < 0 ) {
			$weights{ $gene } = -1;
		} else { print STDERR "Warning! Value is neither <0 nor >0: $exp!\n"; }

		#put the gene into the locked list, if it's not already there
		if ( $locked !~ /$gene/ ) {
			$locked = $locked.$gene."\t"; 
		}

		#set the fold change hash
		$foldchange{ $gene } = $exp;

		#add the gene to the relevant hash
		$relevant{ $gene } = 1;	

	} else { #record which genes weren't used at all
		$unused = $unused.$gene."\n";
	}
}
$fh3->close();

#put the locked gene list into the first position of the microarray array:
$marray[ 0 ] = $locked;

if ($debug) { print "List of locked genes is $locked\n"; }

###########################################################
## If "targets" is set, go through the hash and add each ##
## direct target with only one direction of regulation   ##
## to the stoplist, and its weight to the weights hash   ##                                     
###########################################################

if ( $targets ) {
	foreach my $stopgene ( keys %direct  ) {
		foreach my $gene (keys %{ $direct{ $stopgene } } ) {
			if ($direct{ $stopgene }{ $gene } =~ m/upregulates/ ) {
				$weights{ $gene } = $weights{ $stopgene };
				if ( $stoplist !~ /$gene/ ) {
					$stoplist = $stoplist.$gene."\t";
				}
			} elsif ( $direct{ $stopgene }{ $gene } =~ m/downregulates/ ) {
				$weights{ $gene } = $weights{ $stopgene } * -1;
				if ( $stoplist !~ /$gene/ ) {
					$stoplist = $stoplist.$gene."\t";
				}
			}
		}
	}
	if ( $debug ) { print "Direct targets added; stoplist is now $stoplist\n"; }
}


###########################################################
## This is where the expression values are predicted.    ##
## Genes directly connected to a known misregulated gene ##
## get a heavier weighting in that direction.            ##
###########################################################

print STDERR "Assigning expression values\n";
if ($debug) { print "Assigning expression values\n"; }

my $indexcount = $#marray; #set the indexcount to the length of the microarray array
my $topper = 0; #iteration counter

while ( $topper <= $indexcount ) { #iterate through the array of gene lists
	my ( @trgts ) = split( /\t/, $marray[ $topper ]); #get the gene list from each array entry

	my $tobelocked = ""; #this is the list of genes to be locked after this round
	my $whizzer;

	#go through the gene list and find each gene's partners, assigning their values in the weighting hash
	for $whizzer ( 0 .. $#trgts ) { 
		if ( $stoplist =~/$trgts[ $whizzer ]/ ) { next; } #skip if the gene is in the stoplist		

		if ( exists( $upby{ $trgts[ $whizzer ] } ) ) { #if the gene has any regulators upregulating it
			if ($debug) { print "Upregulators of $trgts[ $whizzer ] are $upby{ $trgts[ $whizzer] }\n"; }

			my ( @TFs ) = split( /\t/, $upby{ $trgts[ $whizzer ] } ); #get the list
			my $bunty;
			for $bunty ( 0 .. $#TFs ) {
				if ( $locked !~ /$TFs[ $bunty ]/ ) { #if the regulator is not yet locked, increment its weight by the weight of its target
					$weights{ $TFs[ $bunty ] } = $weights{ $TFs[ $bunty ] } + $weights{ $trgts[ $whizzer ] };

					#if the regulator regulates a gene with a known expression difference - ie, a gene in the fold-change hash - increment its weight again
					if ( defined( $foldchange{ $trgts[ $whizzer ] } ) ) {
						$weights{ $TFs[ $bunty ] } = $weights{ $TFs[ $bunty ] } + $weights{ $trgts[ $whizzer ] };
					}

					if ( $tobelocked !~/$TFs[ $bunty ]/ ) { $tobelocked = $tobelocked.$TFs[ $bunty ]."\t"; } #if regulator not already in $tobelocked, add it
				}
			}
		} 

		if ( exists( $downby{ $trgts[ $whizzer ] } ) ) { #if the gene has any regulators downregulating it
			if ($debug) { print "Downregulators of $trgts[ $whizzer ] are $downby{ $trgts[ $whizzer] }\n"; }

			my ( @TFs ) = split( /\t/, $downby{ $trgts[ $whizzer ] } ); #get the list
			my $bunty;
			for $bunty ( 0 .. $#TFs ) {
				if ( $locked !~/$TFs[ $bunty ]/ ) { #if the regulator isn't yet locked, increment its weight by the weight of its target

					$weights{ $TFs[ $bunty ] } = $weights{ $TFs[ $bunty ] } - $weights{ $trgts[ $whizzer ] };

					#if the regulator regulates a gene with a known expression difference - ie, a gene in the fold-change hash - increment its weight again
					if ( defined( $foldchange{ $trgts[ $whizzer ] } ) ) {
						$weights{ $TFs[ $bunty ] } = $weights{ $TFs[ $bunty ] } - $weights{ $trgts[ $whizzer ] };
					}

					if ( $tobelocked !~/$TFs[ $bunty ]/ ) { $tobelocked = $tobelocked.$TFs[ $bunty ]."\t"; } #if regulator not already in $tobelocked, add it
				}
			}
		} 
	}

	#make array from tobelocked
	my ( @abouttobelocked ) = split( /\t/, $tobelocked );

	#Go through abouttobelocked array and check for multiple knockouts
	my $ranma;
	for $ranma ( 0 .. $#abouttobelocked) {
		if ( $abouttobelocked[ $ranma ] =~ /;/ ) { #if there is a multiple KO...
			my ( @dkolist ) = split(/;/, $abouttobelocked[ $ranma ]); #get its component genes
			my $akane;
			for $akane ( 0 .. $#dkolist) {

				if ( ( $locked !~ /$dkolist[ $akane ]/ ) && ( $stoplist !~ /$dkolist[ $akane ]/ ) ) { #if the component genes are not already locked or in the stoplist
					if ($debug) { print "Checking component genes of $abouttobelocked[ $ranma ]; $dkolist[ $akane ] is not already locked\n"; }
					#... then add the weight of the multiple knockout to each component gene's weight
					if ( exists( $weights{ $dkolist[ $akane ] } ) && exists( $weights{ $abouttobelocked[ $ranma ] } ) ) {
						$weights{ $dkolist[ $akane ] } = $weights{ $dkolist[ $akane ] } + $weights{ $abouttobelocked[ $ranma ] };
					} else { #Component genes should be in the weights array - catch and warn if they are not
						print STDERR "Warning! $abouttobelocked[ $ranma ] has component gene $dkolist[ $akane ] and either $abouttobelocked[ $ranma ] or $dkolist[ $akane ] does not have an entry in the weights hash\n";
					}
					# and add the components to the abouttobelocked array
					push( @abouttobelocked, $dkolist[ $akane]);
					# and add the ID to the tobelocked string
					$tobelocked = $tobelocked.$dkolist[ $akane]."\t";
				}
			}
		}
	}

	#Go through tobelocked array again, and reset the weight hash values; remove anything that's at 0
	my $misty;
	for $misty ( 0 .. $#abouttobelocked ) {
		if ( $weights{ $abouttobelocked[ $misty ] } < 0 ) { #if the overall consensus is that the gene is downregulated, set it to -1
			$weights{ $abouttobelocked[ $misty ] } = -1;
			$relevant{ $abouttobelocked[ $misty ] } = 1; #add the gene to the relevant hash
			if ($debug) { print "Weight of $abouttobelocked[ $misty ] is -1\n"; }
		}
		elsif ( $weights{ $abouttobelocked[ $misty ] } > 0 ) { #if the overall consensus is that the gene is upregulated, set it to 1
			$weights{ $abouttobelocked[ $misty ] } = 1;
			$relevant{ $abouttobelocked[ $misty ] } = 1; #add the gene to the relevant hash
			if ($debug) { print "Weight of $abouttobelocked[ $misty ] is +1\n"; }
		}
		elsif ( $weights{ $abouttobelocked[ $misty ] } == 0 ) { #if there is no overall consensus, remove it from the abouttobelocked array and the tobelocked string
			if ( $times{ $abouttobelocked[ $misty ] } <= 2 ) { #... and if it's been tested less than three times, leave it at 0 and don't lock it
				$tobelocked =~ s/$abouttobelocked[ $misty ]//; #delete the gene ID from $tobelocked
				$tobelocked =~ s/\t\t/\t/g; #replace any double tabs left over from deleting the gene
				$times{ $abouttobelocked[ $misty ] } = $times{ $abouttobelocked[ $misty ] } + 1; #increment the times hash for that gene
				$relevant{ $abouttobelocked[ $misty ] } = 1; #add it to the relevant hash
			} #if it has been tested 3 times, and there is still no consensus, allow it to be locked at 0.
		}
	}


	#If there's anything in $tobelocked...
	if ( $tobelocked ne "" ) {
		push (@marray, $tobelocked);	#Add tobelocked list to array as a new entry
		$locked = $locked.$tobelocked;  #Add IDs to $locked string
	}

        $indexcount = $#marray; #reset the indexcounter in case a new set of genes was added to the array
        $topper++; #increment iteration counter
}

######################################################################################
## This is where each link is checked for consistency with respect to the weighting ##
######################################################################################

print STDERR "Calculating weights and removing inconsistent links\n";
if ($debug) { print "Calculating weights and removing inconsistent links\n"; }

#square hashes for confirmed regulations - first key acts on second key
#ie if $squpreg{A}{B}=1, A upregulates B
my %squpreg;
my %sqdownreg;

#square hashes for confirmed regulations - first key is acted on by second key
#ie if $squpby{A}{B}=1, A is upregulated by B
my %squpby;
my %sqdownby;

#square hash for multiple genes
my %sqparts;

#For upregulation, the weights of regulator and target should be the same; ie, if x upregulates y and x is 1, y should be 1, while if x is -1, y should be -1. x*y = 1 (or 0, if either has a weight of 0).
#For downregulation, the weights of regulator and target should be opposite; ie, if x downregulates y and x is 1, y should be -1, while if x is -1, y should be 1. x*y = -1 (or 0, if either has a weight of 0).
#0 weights are not accepted. 

#Run through the relevant hash and find the regulators from the upby and downby hashes
foreach my $gene ( keys %relevant ) { 

	if (exists( $upby{ $gene} ) ) {
		my $TF2 = $upby{ $gene }; #rename for ease of use
		if ($debug) { print "Gene is $gene, found upregulators $upby{ $gene}\n"; }

		my ( @tflist ) = split( /\t/, $TF2 ); #get the regulators
		my $minnie;
		for $minnie ( 0 .. $#tflist ) {
			#if the weight of the regulator multiplied by the weight of the regulated gene is >= 0, then they are regulated in the same direction; add to square hash 
			if ( ( $weights{ $tflist[ $minnie ] } * $weights{ $gene } ) >= 0 ) {
				$squpreg{ $tflist[ $minnie ] }{$gene} = 1;
				$squpby{ $gene }{ $tflist[ $minnie ] } = 1;
				if ($debug) { print "$gene is upregulated by $TF2; Squpby is $squpby{ $gene }{ $tflist[ $minnie ] } for $gene and { $tflist[ $minnie ] }\n"; }
			}
		}
	}
	if (exists( $downby{ $gene } ) ) {
		my $TF4 = $downby{ $gene }; #rename for ease of use
		if ($debug) { print "Gene is $gene, found downregulators $downby{ $gene}\n"; }

		my ( @tflist ) = split( /\t/, $TF4 ); #get the tfs
		my $minnie;
		for $minnie ( 0 .. $#tflist ) {
			#if the weight of the regulator multiplied by the weight of the regulated gene is <= 0, then they are regulated in opposite directions; add to square hash 
			if ( ( $weights{ $tflist[ $minnie ] } * $weights{ $gene } ) <= 0 ) {
				$sqdownreg{ $tflist[ $minnie ] }{ $gene } = 1;
				$sqdownby{ $gene }{ $tflist[ $minnie ] } = 1;
				if ($debug) { print "$gene is downregulated by $TF4; Sqdownby is $sqdownby{ $gene }{ $tflist[ $minnie ] } for $gene and { $tflist[ $minnie ] }\n"; }
			}
		}
	}
}

#Finally, run through the weights hash and look for multiple IDs (the ones separated by ";") which have a predicted weight
foreach my $entry ( keys %weights ) {
	if ( $entry =~ /;/ ) {
		my ( @bits ) = split(/;/, $entry );
		my $eva;
		#For each component of a multiple ID, output only those links where the component's weight matches the multiple's weight; ie if x * x;y is positive, or if it is 0
		for $eva ( 0 .. $#bits ) {
			if ( ( $weights{ $entry } * $weights{ $bits[ $eva] } ) >= 0 )  {
				#Add the multiple ID and the matching component to the sqparts hash
				$sqparts{ $bits[ $eva ] }{ $entry } = 1;
			}
		}
	}
}

#print the number of genes found to be linked between stoplist and input microarray
my $relsize = keys %relevant;
print STDERR "Found $relsize linked genes\n";
if ($debug) { 
		print "Found $relsize linked genes\n"; 
		foreach my $link (keys( %weights )) {
			print "Predicted weight of $link is $weights{ $link }\n";
		}
}


########################################################################################################################
## If shortest paths only are desired (default is not), now go through the microarray and for each gene, follow its   ##
## its regulators "up" until the stoplist is reached. Only those genes on a shortest path to the stoplist will be     ##
## output, along with paths connecting misregulated genes to each other if detected during the process.               ##
########################################################################################################################

my %toprint; #hash of genes to print at the end

if ( $short ) { #if shortest paths desired
	print STDERR "Finding shortest paths\n";
	if ($debug) { print "Finding shortest paths\n"; }
	my @findpath;
	foreach my $magenes ( keys %foldchange ) { #go through the foldchange hash and put the genes into the findpath array. This does include the stoplist genes but they are not processed further
		push( @findpath, $magenes);
	}		

	my $confirmed = $stoplist; #paths already connected to stoplist (start with just the stoplist genes)

	my $moose;
	for $moose ( 0 .. $#findpath ) { #for each gene of known regulation
		my @paths; #array to hold paths to look through
		$paths [ 0 ] = $findpath[ $moose ]; #the first gene of the path is the microarray regulated gene
		my @pathstoprint; #array to hold paths for inclusion
		my $stopnow = 0; #counter to signal stop	

		while ( $stopnow < 1 ) { #until a gene on the stoplist or confirmed list is reached
			my @pathstoadd; #array to hold new paths

			my $buzz;
			for $buzz ( 0 .. $#paths ) { #go through the paths array
				my $path = $paths[ $buzz ]; #rename path string for ease of reference
				my ( $gene, $rest ) = split(/\t/, $path ); 	#get the highest-level gene on the path
				if ($debug) { print "Considering $gene for $path\n"; }

				if (( $#paths == 0 ) && ($confirmed =~ /$gene/ )) { #if there is only one path being investigated, and the gene has already been linked, set stopnow and add path to pathstoprint
					push(@pathstoprint, $path);
					$stopnow = 1;
					if ($debug) { print "One path ($path) being investigated, $gene already linked\n"; }
					last; #break the for loop - no need to go any higher
				} 

				if ( exists( $sqdownby{ $gene } ) ) { #if there are downregulators of the gene under consideration
					foreach my $tf ( keys %{ $sqdownby{ $gene } } ) {
						my $newpath = $tf."\t".$path; #make the new path from each regulator of the gene under consideration
						if ( $stoplist =~ /$tf/ ) { #if the regulator is on the stoplist
							push(@pathstoprint, $newpath); #put the new path into the paths to print array
							$stopnow = 1; 
							if ($debug) { print "$tf on stoplist\n"; }
							$confirmed = $confirmed."\t".$newpath; #add the pathway to the confirmed list
						} elsif ( exists($foldchange{$tf} ) ) { #if the TF is itself a misregulated gene
							if ($debug) { print "Adding path to pathstoprint: $newpath\n"; }
							push(@pathstoprint, $newpath); #put the new path into the paths to print array but don't stop the loop
						} else {
							if ($debug) { print "Adding path to pathstoadd: $newpath\n"; }
							push(@pathstoadd, $newpath); #otherwise put the new path into the paths to add array
						} 
					}
				}
				if ( exists( $squpby{ $gene } ) ) { #if there are upregulators of the gene under consideration
					foreach my $tf ( keys %{ $squpby{ $gene } } ) {
						my $newpath = $tf."\t".$path; #make the new path from each regulator of the gene under consideration
						if ( $stoplist =~ /$tf/ ) { #if the regulator is on the stoplist
							push(@pathstoprint, $newpath); #put the new path into the paths to print array
							$stopnow = 1;
							$confirmed = $confirmed."\t".$newpath; #add the pathway to the confirmed list
							if ($debug) { print "$tf on stoplist\n"; }
						} elsif ( exists($foldchange{$tf} ) ) { #if the TF is itself a misregulated gene
							if ($debug) { print "Adding path to pathstoprint: $newpath\n"; }
							push(@pathstoprint, $newpath); #put the new path into the paths to print array but don't stop the loop
						} else {
							if ($debug) { print "Adding path to pathstoadd: $newpath\n"; }
							push(@pathstoadd, $newpath); #otherwise put the new path into the paths to add array
						} 
					}
				}
				if ( exists( $sqparts{ $gene } ) ) { #if the "gene" is a multiple set
					foreach my $tf ( keys %{ $sqparts{ $gene } } ) {
						my $newpath = $tf."\t".$path; #make the new path from each part gene under consideration
						if ( $stoplist =~ /$tf/ ) { #if the part is on the stoplist
							push(@pathstoprint, $newpath); #put the new path into the paths to print array
							$confirmed = $confirmed."\t".$newpath; #add the pathway to the confirmed list
							$stopnow = 1; 
						} elsif ( exists($foldchange{$tf} ) ) { #if the TF is itself a misregulated gene
							push(@pathstoprint, $newpath); #put the new path into the paths to print array but don't stop the loop
						} else {
							push(@pathstoadd, $newpath); #otherwise put the new path into the paths to add array
						} 
					}
				}
			}
				
			if ( @pathstoadd ) { #if pathstoadd array is not empty
				@paths = @pathstoadd; #set the paths array to be pathstoadd, to carry on the iteration
				if ($debug) { print "Pathstoadd not empty; looping\n"; }
			} else { #else set stopnow to 1 to end the loop
				$stopnow = 1;
				if ($debug) { print "Pathstoadd empty; stopping\n"; }
			}
		}

		#now iterate through the paths to print array and put the genes into the toprint hash and add them to the confirmed string
		my $jump;
		for $jump ( 0 .. $#pathstoprint ) {
			my ( @genes ) = split( /\t/, $pathstoprint[ $jump ]); #get the genes in an array
			#go through the array and add each one to toprint hash
			my $mym;
			for $mym ( 0 .. $#genes ) {
				$toprint{ $genes[ $mym ] } = 1;
				if ($debug) { print "Adding $genes[ $mym ] to toprint\n"; }
			}
		}	
	}
	#print the number of genes found to be linked in the shortest paths between stoplist and input microarray
	my $printsize = keys %toprint;
	print STDERR "Found $printsize linked genes in shortest paths\n";
} else {
	foreach my $magenes ( keys %weights ) { #if "short" not set, mark every gene to be printed
		$toprint{ $magenes } = 1;
	}
}

#################################################################################################
## Finally, output the files. Start with the stoplist gene(s) and follow the links downstream. ##
## Only print regulations where both genes are in the toprint hash                             ##
## Print the references from the sourcetoprint hash                                            ##
## Then print the weights, known and predicted, and finally the unused genes.                  ##
#################################################################################################

print STDERR "Checking network and outputting sif file\n";
if ($debug) { print "Checking network and outputting sif file\n"; }

my $sif_output = $prefix.".network_checked.sif";
my $weights_output = $prefix.".final_weights.tab";
my $reference_output = $prefix.".references.attrs"; #tab-separated file, with specific suffix for a Cytoscape edge table
my $unused_output = $prefix.".unused_genes.tab";

#Array to hold the genes to be checked
my @output;
#Start with the stoplist
$output[ 0 ] = $stoplist;

#String to hold the genes already checked
my $considered = "";

#Hash to hold regulations to get sources for
my %sourcetoprint;

#Variables for iterating through the output array
my $iterator = $#output; #set the iterator to the length of the output array
my $mboy = 0; #iteration counter


open (SIFOUT, ">>$sif_output");

while ( $mboy <= $iterator ) { #go through the array
	if ($debug) { print "Gene list is $output[ $mboy ]\n"; }
	my ( @glist ) = split( /\t/, $output[ $mboy ]); #get the gene list from the array entry

	#string for holding genes about to be added to "considered" list
	my $tobeconsidered = ""; 

	#run through the genes in the gene list
	my $miaka; 
	for $miaka ( 0 .. $#glist ) {

		#variable for the gene ID for ease of use
		my $gene = $glist[ $miaka];
		if ($debug) { print "Considering $gene\n"; }

		#if the gene has not already been considered, and if it is in toprint
		if ( ( $considered !~ /$gene/ ) && (exists( $toprint{ $gene } ) ) ) {
			
			#if the gene upregulates downstream genes, print them if they are in the toprint hash
			if ( defined( $squpreg{ $gene }) ) {
				foreach my $sami ( keys( %{ $squpreg{ $gene } } ) ) {
					if ( exists( $toprint{ $sami } ) ) {
						if ($debug) { print "Found targets! $sami\n"; }
						print SIFOUT $gene."\tupregulates\t".$sami."\n";
						$sourcetoprint{ $gene." (upregulates) ".$sami } = 1;

						#if they are not already in the considered or tobeconsidered list, add them to "tobeconsidered"
						if ( ( $tobeconsidered !~ /$sami/ ) && ( $considered !~ /$sami/ ) ) {
							$tobeconsidered = $tobeconsidered.$sami."\t";
						}
					}
				}
			}

			#if the gene downregulates downstream genes, print them if they are in the toprint hash
			if ( defined( $sqdownreg{ $gene })) {
				foreach my $sami ( keys( %{ $sqdownreg{ $gene } } ) ) {
					if ( exists( $toprint{ $sami } ) ) {
						if ($debug) { print "Found targets! $sami\n"; }
						print SIFOUT $gene."\tdownregulates\t".$sami."\n";
						$sourcetoprint{ $gene." (downregulates) ".$sami } = 1;

						#if they are not already in the considered or tobeconsidered list, add them to "tobeconsidered"
						if ( ( $tobeconsidered !~ /$sami/ ) && ( $considered !~ /$sami/ ) ) {
							$tobeconsidered = $tobeconsidered.$sami."\t";
						}
					}
				}
			}

			#if the gene is part of a multiple ID list, print them if they're in the toprint hash
			if ( defined( $sqparts{ $gene } ) ) {
				foreach my $sami ( keys( %{ $sqparts{ $gene } } ) ) {
					if ( exists( $toprint{ $sami } ) ) {
						print SIFOUT $gene."\tis part of\t".$sami."\n";

						#if they are not already in the considered or tobeconsidered list, add them to "tobeconsidered"
						if ( ( $tobeconsidered !~ /$sami/ ) && ( $considered !~ /$sami/ ) ) {
							$tobeconsidered = $tobeconsidered.$sami."\t";
						}
					}
				}
			}

			#Add gene to $considered list
			$considered = $considered.$gene."\t";
			if ($debug) { print "Considered list is now $considered\n"; }
		}
	}

	#if there's anything in the $tobeconsidered string, add it to the output array as a new entry
	if ( $tobeconsidered ne "" ) {
		push (@output, $tobeconsidered);	
		if ($debug) { print "New string added to array: $tobeconsidered\n"; }
	}

        $iterator = $#output; #reset the indexcounter in case a new set of genes was added to the array
        $mboy++; #increment iteration counter
}

#if targets set, there will be direct targets in the stoplist which need to be considered. These will be in the
#direct hash
if ( $targets ) {
        foreach my $stopgene ( keys %direct  ) {
                foreach my $gene (keys %{ $direct{ $stopgene } } ) {
                        if ($direct{ $stopgene }{ $gene } =~ m/upregulates/ ) {
                                print SIFOUT $stopgene."\tupregulates\t".$gene."\n";
                        } elsif ( $direct{ $stopgene }{ $gene } =~ m/downregulates/ ) {
				print SIFOUT $stopgene."\tdownregulates\t".$gene."\n";
                        }
		}
	}
}

close (SIFOUT);

#Print the weights
open (WGTOUT, ">>$weights_output");
print WGTOUT "Gene ID\tfold change\tstatus\n";

#Go through the weights hash and print everything which is in the toprint hash
foreach my $geneentry ( keys %toprint ) {
	my $status = "predicted"; #status of gene - known or predicted expression change

	#reset predicted weights so they are all either 1 or -1
	if ( $weights{ $geneentry } > 1 ) { $weights{ $geneentry } = 1; } #if the weight of the gene is greater than 1, reset to 1
	if ( $weights{ $geneentry } < -1 ) { $weights{ $geneentry } = -1; } #if the weight of the gene is less than -1, reset to -1

	#if the gene has an entry in the fold change hash - ie if it's one of the original misregulated genes with known regulation...
	if ( exists( $foldchange { $geneentry } ) ) { 
		$weights{ $geneentry } = $foldchange{ $geneentry }; #reset the weights hash entry to the known regulation
		$status = "known"; 				    #set its status to "known"
	}

	#print gene ID, weight and status, tab-separated
	print WGTOUT $geneentry."\t".$weights{ $geneentry }."\t".$status."\n";
}
close (WGTOUT);

#Print the references
open(REFOUT, ">>$reference_output");
print REFOUT "References\n";

#Run through the sourcetoprint hash, get the reference from the source hash, and print it
foreach my $refentry ( keys %sourcetoprint ) {
	if ( exists( $source{ $refentry } ) ) {
		print REFOUT $refentry." = ".$source{ $refentry }."\n";
	}
}
close (REFOUT);

#Finally, print the unused genes
open (UNOUT, ">>$unused_output");
print UNOUT $unused;
close (UNOUT);
