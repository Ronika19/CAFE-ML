#!/usr/bin/env perl 

# Run: ./mlcafe.pl --phylogenetic --genus Escherichia --fasta example.fna example.ptt --verbose

use File::Copy;
#use Bio::SeqIO;
use List::Util 'max';
use List::Util 'min';
use List::MoreUtils qw(uniq);
use Getopt::Long qw(GetOptions);

OPTIONS();

# If output filename is provided by user
if (defined $outfile) {
	$output=shift(@ARGV); #print("$output\n");
	$name = $output; 
	#$name =~ s/NC_004431\///;
}

#If sequence file is not provided
if (!defined $fasta) {
	if (scalar(@ARGV)<1){
	print("Sequence file not provided\n");
	exit;
	}
}

#If sequence file is provided
else {
	$f = shift @ARGV;
	$in_filename = $f;
	$name = $f if(!defined $outfile);
	$name =~ s/\..*//; 
	#$name =~ s/NC_004431\///; #print("$name\n");
	#Check if sequence file exists
	unless (-e $f) {
		print("Sequence file does not exist\n");
		exit;
	} 
	$nce = 1;
	#Check if sequence file is in fasta format	
	open(F,$f);
	$line = <F>;
	if ($line !~ /^>/){
		print("Incorrect sequence file. Please provide sequence file in fasta format\n");
		exit;
	}
	#Get sequence from file
	while($line = <F>){
		$line =~ s/\n//;
		if($line !~ /^>/){
			$dna .= uc($line);
		}
	}
}

#Check only ATGC bases are present in the genome
@bases = ('A','T','G','C');
@genome = split("",$dna);
foreach(@genome){
	$s = $bases[rand@bases];
	$_ =~ s/[^ATGC]/$s/g;
}

open (OUT, ">Output_Files\/$name\_CAFE_1"); #print("$name\n");
$m = join("",@genome);
print(OUT "$m");
$lenm = length($m);
print("Bacterial genome loaded\nLength of the genome is $lenm bp\n") if ($verb == 1); 

#If annotation file is provided
if (!defined $annotation){
	if (scalar(@ARGV) < 1){
		print("Annotation file does not exist\n");
		exit; 
	}

	$f1 = shift @ARGV;

	#Check if annotation file exists
	unless (-e $f1) {
		print("Annotation file does not exist\n");
		exit;
	} 
	
	open(F1,$f1);
	@annot_ptt = <F1>;
	
	#Check if annotation file is in ptt format
	$a = 0; $b = 0; $c = 0;
	foreach $annot_ptt(@annot_ptt) {
		@C = split(/\t/,$annot_ptt);
		$coords = $C[0];
		if (($coords =~ m/\.\./) && ($#C == 8)) {
			($start,$end) = split('\.\.',$coords); #print("$start\t$end\n");
			$arr1[$a] = $start;
			$arr2[$b] = $end;
			$arr3[$c] = $C[8]; # Gene Product
			$arr4[$d] = $C[3]; # PID
			$a++; $b++; $c++; $d++;
		}
	}
	
	@identifiers = split('\t',$annot_ptt[2]);

	if (defined $identifiers[0] && defined $identifiers[1] && defined $identifiers[8]) {
		if (($identifiers[0] !~ /Location/) | ($identifiers[1] !~ /Strand/) | ($identifiers[8] !~ /Product/)) {
			print("Incorrect annotation file. Please provide annotation file in ptt format\n");
			exit;
		}
	}

	else {
		print("Incorrect annotation file. Please provide annotation file in ptt format\n");
		exit;
	}

	#Get annotations from ptt file
	open (PTTOUT, "Output_Files\/>$name\_gene_coord");
	$gene_counter = 0;
	for($i = 0; $i < $a; $i++){
		chomp $arr3[$i];
		print(PTTOUT "$arr3[$i]\t$arr1[$i]\t$arr2[$i]\t$arr4[$i]\n");
		$st = $arr1[$i]; $en = $arr2[$i]; $fn = $arr3[$i]; $prot_id = $arr4[$i];
		$gene_start[$gene_counter] = $st; $gene_end[$gene_counter] = $en; $gene_func[$gene_counter] = $fn;
		$pid[$gene_counter] = $prot_id;
		$gene_counter += 1;
	}
}

#If phylo module is used
if ($phylo == 1){
	print("Checking phylogenetic distribution of genes\n") if ($verb == 1);
	print("$input_genus\n"); print("\n\n");
	#shift(@ARGV);
	system("cat ./faa/* > Output_Files/faa_database");
	@faa_files = <./faa/*>;
	$expect = scalar(@faa_files); #print("$expect\n");
	if ($expect < 4) {
		print("Phylo module requires atleast 5 genomes for comparison\n");
		exit;
	}

	if(defined $annotation) {$infile = $name;}
	##make database
	system("makeblastdb -dbtype prot -in Output_Files/faa_database");#faa_database
	system("blastp -db Output_Files/faa_database -query $infile\.faa -outfmt \"6 qseqid sseqid stittle salltitles qcovs pident\" -out Output_Files/$name\_blast_output -num_threads 4") if(defined $annotation); 
	system("blastp -db Output_Files/faa_database -query $name\.faa -outfmt \"6 qseqid sseqid stittle salltitles qcovs pident\" -out Output_Files/$name\_blast_output -num_threads 4") if(!defined $annotation);  
	$blastf = "Output_Files\/$name\_blast_output";
	open(BL, $blastf);
	open (BLOUT, ">Output_Files\/$name\_phy"); open(PHYOUT, ">Output_Files\/$name\_phyout");
	
	$bf = 0; open(QACC, ">Output_Files\/$name\_query_accession");
	while (<BL>){
		chomp;
		@BFC = split('\t', $_);
		$queryacc[$bf] = $BFC[0]; print(QACC "$BFC[0]\n");
		#$subacc[$bf] = $BFC[1];
		$blsub[$bf] = $BFC[2];
		$coverage[$bf] = $BFC[3];
		$identities[$bf] = $BFC[4];
		$bf++;
	}
	$phprev="";
	$phcount=0;
	$phcum_count=0;
	$pi=0;
	foreach(@queryacc) {
		if($phprev ne $_) {
			if($phcount) {
				#printf("%s:%d \n",$_,$phcount);
				$phash{$phprev} = $phcount;
				$phash_cum{$phprev} = $phcum_count;
				#print("$_ $phcount $cum_count\n");
			}
			$ptemp_count = $phcount; #save values for last query
			$ptemp_cum_count = $phcum_count;
			$phprev = $_;
			$phcount = 0;
		}
		$phcount++;
		$phcum_count++;
		$pi++;
	}

	for ($i=0; $i<$bf; $i++) {
		if ($blsub[$i] =~ m/\[(.*)\]/){ #get blast subject names
			$blastname[$i]=$1; #print("$blastname[$i]\n");		
		}
	}

	for ($i=0; $i<$bf; $i++){
		$cur = $blastname[$i];
		if (defined $cur){
			@blast = split('\s+',$cur);
		}
		else { $blastname[$i] = ''; }
		$blast_genus[$i] = $blast[0];
		#$blast_sp[$i] = $blast[1];
		$blast_strain[$i] = join('', @blast[2..$#blast]);
		#print("$i $blast_strain[$i]\n");
		#print("$blast_genus[$i]\t$blast_sp[$i]\t$blast_strain[$i]\n");
	}
	
	@uniq_queryacc = uniq(@queryacc);

	# Extrating Protein Ids from faa file
	for ($i=0; $i<scalar(@uniq_queryacc); $i++) {
		$protein_identifier[$i] = (split('\|', $uniq_queryacc[$i]))[1]; #print("$uniq_queryacc[$i]\t$protein_identifier[$i]\n");
	}

	$last_element = $uniq_queryacc[-1];
	$phash{$last_element} = $ptemp_count;
	$phash_cum{$last_element} = $ptemp_cum_count;

	$spcount = 0;
	@all_strains = '';
	$phv = 0;
	
	for ($k=0; $k<scalar(@uniq_queryacc); $k++){
		$cur_query = $uniq_queryacc[$k];
		$next_query = $uniq_queryacc[$k+1]; #print("$cur_query\t$next_query\n");
		pop(@all_strains);
		$start = $phash_cum{$cur_query}-$phash{$cur_query};
		$end = $phash_cum{$cur_query};
		#print("$start\t$end\t$input_genus\n");
		for ($i=$start; $i<$end; $i++){
			
			if ($blast_genus[$i] =~ m/$input_genus/i){
				print(PHYOUT "$i\t$cur_query\t$phash{$cur_query}\t$coverage[$i]\t$identities[$i]\t$blast_strain[$i]\n");
				if (($coverage[$i] >= 70) && ($identities[$i] >= 70)) { #if identities and coverage are greater than 70
					$spcount += 1;
					push(@all_strains,$blast_strain[$i]);
					#print("$i\t$cur_query\t$spcount\t$phash{$cur_query}\t$blast_strain[$i]\n");
				}
			}
		}
		#print("@all_strains\n");
		$actual_count = scalar(uniq(@all_strains)); #Do not count multiple times if query matches multiple genes in same genome/strain
		
		if (!defined $expect ) {$pc=10000;}
		else { $pc = $actual_count/$expect; }
		if ($pc >= 0.5){
			$value = 0; #typical phyletic pattern
			$phyvalue[$phv] = 0;
			$phv++;
		}
		else {
			$value = 1; # atypical phyletic pattern
			$phyvalue[$phv] = 1;
			$phv++;
		} 
		
		if (!defined $expect) {
			print(BLOUT "$cur_query\t0\t$actual_count\t$phash{$cur_query}\t$value\n");
		}
		
		else {
			print(BLOUT "$cur_query\t$expect\t$actual_count\t$phash{$cur_query}\t$value\n");
		}
		
		$spcount=0;
		@all_strains='';
	}
}


open (PYO1, ">Output_Files\/$name\_PhyGenes"); #print "gene counter $gene_counter\n";
for ($l=0;$l<$gene_counter;$l++) {
	if (grep(/^$pid[$l]$/, @protein_identifier)) {
		for ($x=0; $x<scalar(@protein_identifier); $x++) {
			if ($pid[$l] == $protein_identifier[$x]) {
				print(PYO1 "$gene_func[$l]\t$gene_start[$l]\t$gene_end[$l]\t$phyvalue[$x]\n");
			}
		}
	}
	else {
		print(PYO1 "$gene_func[$l]\t$gene_start[$l]\t$gene_end[$l]\t0\n");
	}
}
close PYO1;

@Dlibrary=qw(transposase transposon integrase integration phage prophage bacteriophage mobile mobility insertion recombinase plasmid resolvase);
open (PYOUT1, "Output_Files\/$name\_PhyGenes");
@pyout1 = <PYOUT1>; $marker_value = 0;
for ($i=0; $i<scalar(@pyout1); $i++) {
	chomp $pyout1[$i];
	@split_pyout1 = split("\t", $pyout1[$i]);
	$gfunc[$i] = $split_pyout1[0]; $k = $i+1;
	$phyval[$i] = $split_pyout1[3];
	$gstart[$i] = $split_pyout1[1]; $gend[$i] = $split_pyout1[2];

	$gfunc[$i] =~ s/[^a-zA-Z0-9]/ /g; #@uniq_genes = ();
	for $Dlibrary (@Dlibrary) {
	
		if ($gfunc[$i] =~ /$Dlibrary/ig) {
			$marker_value = 1;
			$string = "$k\t$gstart[$i]\t$gend[$i]\t$gfunc[$i]\t$phyval[$i]\t$marker_value\n";
			push (@genes, $string); #print("$string\n"); 
			@uniq_genes = uniq(@genes); 
			}
		else {
			$marker_value = 0;
			$string = "$k\t$gstart[$i]\t$gend[$i]\t$gfunc[$i]\t$phyval[$i]\t$marker_value\n"; 
			push (@genes, $string); #print("$string");
			@uniq_genes = uniq(@genes); 
			}
		}
	}
close PYOUT1; close PYOUT2;
#print(scalar(@uniq_genes), "\n"); 
open (PYOUT3, ">Output_Files\/$name\_AllGenes.txt");
for ($i=0; $i<scalar(@uniq_genes); $i++) {
	chomp $uniq_genes[$i];
	@split_pyout3 = split("\t", $uniq_genes[$i]);
	$gfuncn[$i] = $split_pyout3[6]; $ks[$i] = $split_pyout3[0];
	if ($ks[$i] == $ks[$i-1]) { 
		$uniq_genes[$i-1] = "";
		}
	}
#print(scalar(@uniq_genes), "\n"); #print("@uniq_genes\n"); 
for ($i=0; $i<scalar(@uniq_genes); $i++) {
	chomp $uniq_genes[$i];
	print(PYOUT3 "$uniq_genes[$i]\n") unless ($uniq_genes[$i] == "");
	}
close PYOUT3;
=c
if ($expert==0){
	system ("rm Output_Files\/$name\_CAFE_1");	
	#system ("rm Output_Files\/$name\_CAFE.ptt"); 
	#system ("rm Output_Files\/$name\_CAFE.faa");
	system ("rm Output_Files\/$name\_phy") if($phylo==1); 
	system ("rm Output_Files\/$name\_blast_output") if($phylo==1);	
	system ("rm Output_Files\/$name\_phyout") if($phylo==1);	
	system ("rm Output_Files\/>$name\_gene_coord");
	system ("rm Output_Files\/faa_database") if (defined $phylo);	
	system ("rm Output_Files\/$name\_PhyGenes") if (defined $annotation);
	system ("rm Output_Files\/$name\_query_accession");
}
=cut

if (defined $annotation) {
	print("Completed scanning for genomic island specific marker genes!\n"); 
}

sub information {
	print("Aberrant Phyletic Pattern & Marker Enrichment Analysis\n");
	exit;
}

sub usage {
	print(STDERR "\nUsage:\n  $0 [options] --phylogenetic --genus [genus name] --fasta example.fna example.ptt --verbose\n\n");
  	exit;
}

sub verbose {
	$verb=1; print("Verbose Activated!\n");
}

sub expert {
	$expert=1; print("Keeping temporary files!\n");
}

sub phylogenetic {
	$phylo=1; print("Phylogenetic Module Activated!\n");
}

sub genus {
	$input_genus=shift(@ARGV);
	#print "Input genus $input_genus\n";
	return $input_genus;
}

sub OPTIONS {
	GetOptions(
	"help" =>\&usage,
	"information" =>\&information,
	"annotation" =>\$annotation,
	"fasta" =>\$fasta,
	"phylogenetic" =>\&phylogenetic,
	"genus" =>\&genus,
	"out" =>\$outfile,
	"verbose" =>\&verbose,
	"expert" =>\&expert,
	);
}

