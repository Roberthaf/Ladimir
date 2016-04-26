##### Program description
#
# Authors:
#	 Robert Anton Hafthorsson
#	
# Description:
#	This program is intended to compaire LADs and miR's and see if they(miR's)
# 	are always outside the LADs or with in.
#
# Modules:
#	none
#		
# Procedure:
#	1. Read in the LADs reagion and miRNA data
#	2. Process the inputfiles
#		1) Procces cLAD files. Find chromosome and positions of cLADs.
#		2) Procces miRNA files. Identify miR´s
#	3. Search for miRNA that lie within a cLAD region.
#	4. Analizes and create statistics and output files 
#
#
#		ATH: Hvernig output? Hvernig á að skoða þeta.
#
#
# usage: Ladimir_v3.pl [LAD file] [miRNA file] [Gene file] "name of output"

use strict;
use Data::Dumper;
use Array::Utils qw(:all);
use Getopt::Std;
use DateTime;
use experimental 'smartmatch';

## Read in the data
my %options=();
getopts("c:l:m:g:o:hhelp", \%options);
print "-l $options{l}\n" if defined  $options{l};
print "-m $options{m}\n" if defined  $options{m};
print "-g $options{g}\n" if defined  $options{g};
print "-o $options{o}\n" if defined  $options{o};
print "-c $options{c}\n" if defined  $options{c};

print "\n"; 

print "Other command options found on the command like :\n" if $ARGV[0];
my $cutoff = $options{c};
if (!defined $options{c}){
	$cutoff= 10000000;
}

my $lad_file = $options{l};
my $mir_file = $options{m};
my $gene_file = $options{g};
my $animal = $options{o};

my $outfile;
my $outfile2;
my $outfile3 = "testfile.txt";

if (defined $options{o}){
	$outfile = "$animal.txt";
	$outfile2 = "$animal.html";
}else{
	print "outfile is not specified\n";
	$outfile= "outfile.txt";
}
if ($options{h} ) {
print <<EOF;



				LADIMIR - A search tool for Lamina acosiated domain and miRNA research.
				
SYNOPSIS
		ladimir cLADfile miRNAfile GeneFile outputfile
	
	ladimir [-l cLAD file] [-m miRNA file] [-g Gene file] [-o Outputfile]
	
DESCRIPTION
	
	Ladimir reads in cLAD,miRNA and gene files and from them analizes if miRNA's 
	lie within the cLAD regions and what lie up and down stream of each miRNA. 	
	
	
	
	For the [name of outputfiles] no file extentions are needed.
	Ladimir will create a file containing miRNA's that are within cLADS
	And a .gff3 file with a list of miRNA that where used to do the comparison. 

	The options for the program are as follows:
	
	-h -help	Print this help file and exit
	
	-o 			Create all output files in a specified directory.
				The directory must exist as the program will not
				create a new one. If not output name is give files
				will be created wtih no name.
	
	-l 			Specifies a non-default file which must containing
				the cordinates for lamina associated domain (LAD's)
				Valid format is bed.
				
	-m			Specifies a non-default file which must contain the
				miRNA information. Valid formats are gff3, gff2 and gtf.

	-g 			Specifies a non-default file which must contain the
				gene information. Valid formats are gff3,gff2 and gtf. 
	
	-c			Select a cut off distance for upstream and downstream genes.
				If no number is specified then range is 
EOF
exit;
## Exit after help information have been displayed
}else{

## Arrys for input
my @input_lad;
my @input_mir;
my @input_genes;

## Open the file's
open ( my $LADFILE, '<', $lad_file )
	or die "No file LAD file found: $lad_file";
	
open ( my $MIRFILE, '<', $mir_file )
	or die "No miRNA file found: $mir_file";

open ( my $GENEFILE, '<', $gene_file )
	or die "No gene file found: $gene_file";	
	
## Files to write out
open my $OUT, '>', $outfile or die "Cannot open: $outfile";
open my $OUT2, '>', $outfile2 or die "Cannot open: $outfile2";

##########PRINT OUT FILE USED FOR TROUBLE SHOOTING. REMOVE LATER
open my $OUT3, '>', $outfile3 or die "Cannot open: $outfile3";
##########


## Print Out information to command line for user.
print "processing LAD information\n"; 


################################################## 
#### Start parsing input files. 
## LAD file.
while ( my $line = <$LADFILE> ){
	chomp $line;
	my ($chr, $start, $end) = split (/\s/, $line);
	push @input_lad, [$chr, $start, $end];
}

## Parse Gene file
print "processing Gene File\n";
while ( my $line = <$GENEFILE> ){
	chomp $line;
	if ( $line =~ /biotype=protein_coding/ ){
		if ( $line =~ /(.{1,3}?)\t.+?\tgene\t(\d+)\t(\d+).+Name=(.+?);\b/gm ) {
				#push into the array the gene name and start and end positions
				push @input_genes, ["chr$1", $2, $3, $4];
		}
	}else{
	#Skip over lines which have no 
	next;
	}
}
print "Found\t", scalar @input_genes ," Genes in file\n\n";
close $GENEFILE;

my @indata = <$MIRFILE>;
my $species;
#Start processing miRNA information.

print "processing miRNA information\n";
if( (@indata[3]) =~/.+Mus musculus.+/ ) {
   $species = 'Mus musculus';
   print "File is from Mus musculus\n";
    foreach  my $line (@indata ){
		chomp $line;
		if ($line =~ /(\d{1,2}|\w)\t.\tmiRNA_primary_transcript/){
			# This is for printing out the names of the miRNA that are inside LADS
			$line =~ /(\d{1,2}|\w)\t.\tmiRNA_primary_transcript\t(\d+)\s(\d+).+Name=(.+\b)/;
			my $mir_name = substr $4,4,(length $4);
			push @input_mir, ["chr$1", $2, $3, $mir_name];
		}
	}
}
elsif( (@indata[3]) =~/.+Drosophila melanogaster.+/ ){
      $species = 'Drosophila melanogaster'; 
   print "File is from Drosophila M.\n";
    foreach  my $line (@indata ){
		chomp $line;
		if ($line =~ /(.{1,5})\t.\tmiRNA_primary_transcript/){
		# This is for printing out the names of the miRNA that are inside LADS
			$line =~ /(.{1,5})\t.\tmiRNA_primary_transcript\t(\d+)\s(\d+).+Name=(.+\b)/;
			my $mir_name = substr $4,4,(length $4);
			my $chr_name = $1;
			$chr_name =~ s/2RHet/2R/;	
			$chr_name =~ s/2LHet/2L/;	
			push @input_mir, ["chr$chr_name", $2, $3, $mir_name];
		}
	}
}
elsif( (@indata[3]) =~/.+Homo sapiens.+/ ){
   $species = 'Homo sapiens'; 
   print "File is from Human\n";
   foreach  my $line (@indata ){
		chomp $line;
		if ($line =~ /(.{1,6})\t.\tmiRNA_primary_transcript/){
			# This is for printing out the names of the miRNA that are inside LADS
			$line =~ /(.{1,6})\t.\tmiRNA_primary_transcript\t(\d+)\s(\d+).+Name=(.+\b)/;
			# push @input_mir, [$1, $2, $3, $4];
			my $mir_name = substr $4,4,(length $4);
			push @input_mir, [$1, $2, $3, $mir_name];
		}
	}
}
elsif( (@indata[3]) =~/.+Caenorhabditis elegans.+/ ){
   $species = 'Caenorhabditis elegans'; 
   print "File is from C.Elegans\n";
    foreach  my $line (@indata ){
		chomp $line;
		if ($line =~ /(.{1,6})\t.\tmiRNA_primary_transcript/){
			# This is for printing out the names of the miRNA that are inside LADS
			$line =~ /(.{1,6})\t.\tmiRNA_primary_transcript\t(\d+)\s(\d+).+Name=(.+\b)/;
			# push @input_mir, [$1, $2, $3, $4];
			my $mir_name = substr $4,4,(length $4);
			push @input_mir, [$1, $2, $3, $mir_name];
		}
	}
}else{
print "Species is unkown. Please make sure that species name is in latin and one of four accepted
 species\nHomo sapiens, Mus musculus, Caenorhabditis elegans,Drosophila melanogaster\n
 In a correctly formated file species information is on the 4th line in the miRNA file.";
}
close $MIRFILE;
#Start the comparison algorithm
my $lad_scalar = scalar @input_lad;
my $mir_scalar = scalar @input_mir;
my @mir_os_lad;
my @mir_wi_lad;
my @final_array;

# Itterate through all the miR's
for ( my $i = 0;$i < $mir_scalar;$i++){
	# Itterate through all the LADs
	for (my $j = 0;$j < $lad_scalar;$j++){
		my $temp_length = $input_lad[$j][2] - $input_lad[$j][1]; #Find current LADs Length
	#Check if chromosome is the same nr.
		if ($input_lad[$j][0] eq $input_mir[$i][0] ){
		#Build Multi dimentional array that has the information on miR's inside of LADS.
			if ( ( ( $input_mir[$i][1] - $input_lad[$j][1]) >= 0 ) and ( ($input_mir[$i][1] - $input_lad[$j][1]) <= $temp_length ) ){
				push @mir_wi_lad,[ $input_mir[$i][0],$input_mir[$i][1],$input_mir[$i][2],$input_mir[$i][3] ];			
			}
			##Check if miRNA lies on partially inside of a cLAD region. Start of miRNA is outside and end is inside. If yes
			##its counted as inside LAD region.
			if ( ( ( $input_mir[$i][1] - $input_lad[$j][1]) <= 0 ) and  ( ( $input_mir[$i][2] - $input_lad[$j][1] ) >= 0 ) ){
				push @mir_wi_lad,[ $input_mir[$i][0],$input_mir[$i][1],$input_mir[$i][2],$input_mir[$i][3] ];
			}
			## New iterations check if miRNA is partially within in cLAD region near the end.
			if ( ( ( $input_mir[$i][1] - $input_lad[$j][2]) <= 0 ) and  ( ( $input_mir[$i][2] - $input_lad[$j][2] ) >= 0 ) ){
				push @mir_wi_lad,[ $input_mir[$i][0],$input_mir[$i][1],$input_mir[$i][2],$input_mir[$i][3] ];
			}
		}
	}
	#Start processing the Gene file.
	my @temp_arrayH;
	my @temp_arrayL;
	for ( my $j = 0;$j < scalar @input_genes;$j++ ){
		 if ( $input_mir[$i][0] eq $input_genes[$j][0] ){			
			my $temp_number = 0;
			$temp_number = ( $input_mir[$i][1]-$input_genes[$j][1] );		
			
			if ( $temp_number <= 0 ){
				push @temp_arrayL, [$temp_number, $input_mir[$i][3],$input_genes[$j][3],$input_mir[$i][0] ];
			}
			if ( $temp_number >=0 ){
				push @temp_arrayH, [ $temp_number, $input_mir[$i][3],$input_genes[$j][3],$input_mir[$i][0] ];
				
			}
		}
	}
	push @final_array , [ $temp_arrayH[-1],$temp_arrayL[0] ];
	#printing the percentage of the progress
	my $per = int(($i/$mir_scalar)*100);
	print "\b\b\b\b$per%";
}
	print $OUT3 @input_mir;
	print "\b\b\b\bComputation is complete\n";
	#Create a array's from known mir inside lad containing only names
	my @mir_wi_lad_name;
	my @mir_wi_lad_chr;
	#Take out the miRNA names.
	for my $i (@mir_wi_lad){
		push @mir_wi_lad_name, @$i[3];
	}
	my @mir_all_lad_name;
	for my $i (@input_mir){
		push @mir_all_lad_name, @$i[3];
	}
	#Check to see which miRNA are on the outside or between catagory
	my @mir_os_lad_new = array_minus (@mir_all_lad_name,@mir_wi_lad_name); 
	#Create an array which contains the miR's that are outside of LAD
	my @mir_os_lad_temp;
	foreach my $line (@input_mir){
		for (my $i=0; $i < scalar @mir_os_lad_new; $i++){
			if ( @$line[3] eq $mir_os_lad_new[$i] ){	
				push @mir_os_lad_temp, ("@$line[0]:@$line[1]:@$line[2]:@$line[3]");
			}
		}
	}
	#Remove Dublicate miR's.
	my %dublicate;
	my @outside_mir_unique = grep { !$dublicate{$_}++ } @mir_os_lad_temp;
	foreach my $line (@outside_mir_unique){
		my ($chr,$spst,$epst,$name) = split /:/, $line;
		push @mir_os_lad, [ $chr, $spst, $epst, $name];
	}
	#Start printing out
	print "Number of miR's outside of LADS regions: ", scalar @mir_os_lad,"\n";
	print "Number of miR's within LADS: ", scalar @mir_wi_lad ,"\n";
	print "miRNA total numbers: ", scalar @mir_all_lad_name, "\n";
	print $OUT "Number of miR's outside of LADS regions: ", scalar @mir_os_lad,"\n";
	print $OUT "Number of miR's within LADS: ", scalar @mir_wi_lad ,"\n";
	print $OUT "miRNA total numbers: ", scalar @mir_all_lad_name, "\n";
	##Here I will start the stastistical analysis
	my %mir_count; 
	my %mir_count_out;
	##Here I check to see the number of miR outside
	foreach my $line (@mir_os_lad){
		$mir_count_out{$$line[0]}++;
	}
	my @chrom_names = keys %mir_count_out;
	##Here I check to see the number of miR within
	foreach my $line (@mir_wi_lad){	
		$mir_count{$$line[0]}++;
	}
	#Here i save every chromosome name to an array 
	foreach my $k (@chrom_names){
		if ( defined $mir_count{$k} ){
		next;
		}else{
		$mir_count{$k} = 0; 
		}
	}

	#Do calculation for miRNA within cLADS
	my @chr_number_wi;
	my @chr_number_ot;
	print "\nThis tables show how many miRNA are inside cLADs each chr\n";
	print $OUT "This tables show how many miRNA are inside cLADs each chr\n";
	#Subroutine to fix a comma at the end of an element when writing out a javascript variable.

	sub endfix{
		my $fix = $_[0];
		$fix = substr $fix, 0,(length($fix)-1); 
		return $fix;
	}
#Variables for printing out various information
my $print_out_wi_name;my $print_out_wi_numb;
my $print_out_ot_name;my $print_out_ot_numb;
#Start sorting miRNA information
foreach my $k (sort {substr($a, 3) <=> substr($b, 3) || substr($a, 3) cmp substr($b, 3)} keys %mir_count){
	push @chr_number_wi, $mir_count{$k};
	print $OUT "$k $mir_count{$k}\n";
	$print_out_wi_name .= "'$k',";
	$print_out_wi_numb .= "$mir_count{$k},";
	print "$k = $mir_count{$k}\n";
}
# my $test= endfix($print_out_wi_numb);

$print_out_wi_name = endfix($print_out_wi_name);
$print_out_wi_numb = endfix($print_out_wi_numb);

# Do calculation for miRNA outside cLADS
print "\nThis tables show how many miRNA are outside cLADs each chr\n";
print $OUT "\nThis tables show how many miRNA are outside cLADs each chr\n";
foreach my $k (sort {substr($a, 3) <=> substr($b, 3) || substr($a, 3) cmp substr($b, 3)} keys %mir_count_out){
	push @chr_number_ot, $mir_count_out{$k};
	print $OUT "$k $mir_count_out{$k}\n";
	$print_out_ot_name .= "'$k',";
	$print_out_ot_numb .= "$mir_count_out{$k},";
	print "$k = $mir_count_out{$k}\n";	
}
# foreach my $k (sort {substr($a, 3) <=> substr($b, 3) || substr($a, 3) cmp substr($b, 3)} keys %gene_count){	
# }
$print_out_ot_name = endfix ($print_out_ot_name);
$print_out_ot_numb = endfix ($print_out_ot_numb);

#Start printing out for text file
print $OUT "\nA list of miR's that Ladimir identidies as within cLADS"; 
print $OUT "\nChr name\tStart pst\tEnd pst\tmiR name\n";

#Starting creating javascript variables.
my $Datatable_info;	
my $Datatable_info_ot;
foreach my $i(@mir_wi_lad){
	print $OUT "@$i\n";
	$Datatable_info .= "\[\"$$i[3]\",\"$$i[1]\",\"$$i[2]\",\"$$i[0]\"\],";	
}
$Datatable_info = endfix ($Datatable_info);
$Datatable_info_ot = endfix ($Datatable_info_ot);

print $OUT "\nmiRNA outside of lads\n";
foreach my $i(@mir_os_lad){
	print $OUT "@$i\n";
	$Datatable_info_ot .= "\[\"$$i[3]\",\"$$i[1]\",\"$$i[2]\",\"$$i[0]\"\],";
}
my $print_genes;
print $OUT "\nExamin which genes are next to miRNA's";
my $genereport=0;

foreach my $g(@final_array){
	my $mir_name_temp;
	my $temp_chromo;
	#make sure there is always a name for the mir or the chromosome
	if(defined $$g[0][1]){
	$mir_name_temp = $$g[0][1];
	}else{
	$mir_name_temp = $$g[1][1];
	}
	if(defined $$g[0][3]){
	$temp_chromo = $$g[0][3];
	}else{
	$temp_chromo = $$g[1][3];
	}
	
	my $within;
	if ( $mir_name_temp ~~ @mir_wi_lad){
		$within = "Within LAD";
	}else{
		$within = "Outside LAD";
	}
	my $cutoff_down = -1*$cutoff;
	
	if ( $$g[0][0] <= $cutoff and $$g[1][0]>= $cutoff_down) {
		if ($mir_name_temp eq ''){
		# next;
		}else{
		my $temp_Number = -1*$$g[1][0]; #Change the number to positive.
		$print_genes .= "\[\"$mir_name_temp\",\"$$g[0][2]\",\"$$g[0][0]\",\"$$g[1][2]\",\"$temp_Number\",\"$temp_chromo\"\,\"$within\"],";
		print $OUT "miRNA name ",$mir_name_temp,"\t Upstream Gene ", $$g[0][2],"\t", $$g[0][0],"\t Downstream Gene ", $$g[1][2],"\t", $temp_Number,"\tChromosome ", $temp_chromo ,"\tLAD localazation\s", $within,"\n";
		$genereport++;
		}
	}
}
$print_genes = endfix ($print_genes);

###### START GENE DENSITY
	my %gene_count;
	my $print_genedensity_names;
	my $print_genedensity_numbers;
	
	foreach my $k ( @input_genes){
		if ( $$k[0] ~~ @chrom_names){
			$gene_count{$$k[0]}++;
		}
	}
	# Sorting and printing out geneomic information
	my @gene_quantity;
	foreach my $k (sort {substr($a, 3) <=> substr($b, 3) || substr($a, 3) cmp substr($b, 3)} keys %gene_count){
		$print_genedensity_names .= "\"$k\","; 
		push @gene_quantity, $gene_count{$k};
	}
	$print_genedensity_names = endfix ($print_genedensity_names);

	##### Gene Density Calculations
	my @chr_size;
	##known chromosome sizes in Mega base pairs.
	if ($species eq 'Homo sapiens'){
		@chr_size= ( 155,59,249,243,198,191,180,171,159,146,141,135,135,133,114,107,102,90,81,78,59,63,48,51 );
	};
	if($species eq 'Mus musculus'){
		@chr_size = ( 162,195,182,160,157,152,150,145,129,125,131,122,120,120,125,104,91,95,91,61 );
	};
	if($species eq 'Drosophila melanogaster'){
		@chr_size = ( 22.4,23.0,21.2,24.4,27.9,1.4 );
	};
	if($species eq 'Caenorhabditis elegans'){
		@chr_size = ( 14.9,15.3,13.8,18.5,20.9,17.1 );
	};
	my @gene_density;
	for ( my $i=0; $i < scalar @gene_quantity; $i++){
		my $gene_density_temp;
		$gene_density_temp =  sprintf ("%.1f", $gene_quantity[$i]/$chr_size[$i] );
		$print_genedensity_numbers .= "$gene_density_temp,";
	}
	$print_genedensity_numbers = endfix ($print_genedensity_numbers);

##### START HTML document printing.
my $dt = DateTime->now;
print $dt,"\n";

my $mir_ot_number = scalar @mir_os_lad;my $mir_wi_number = scalar @mir_wi_lad;
my $mir_all_number = scalar @mir_all_lad_name;
my $genes_number = scalar @input_genes;
my $precent = int(($mir_ot_number/$mir_all_number)*100);

print $OUT2 "<!DOCTYPE html>
<html>
  <head>
  <meta content=\"text/html;charset=utf-8\" http-equiv=\"Content-Type\">
  <meta content=\"utf-8\" http-equiv=\"encoding\">
  <script src=\"https://ajax.googleapis.com/ajax/libs/jquery/1.12.0/jquery.min.js\"></script>
  <script src=\"https://code.highcharts.com/highcharts.js\"></script>
  <script src=\"https://cdn.datatables.net/1.10.11/js/jquery.dataTables.min.js\"></script>
  <link rel=\"stylesheet\" href=\"https://cdn.datatables.net/1.10.11/css/jquery.dataTables.min.css\">
  <script src=\"http://code.highcharts.com/modules/exporting.js\"></script>

  </head>
	<style type=\"text/css\">
	body {font-family: \"Times New Roman\";font-size: 20px;}
	header{position:fixed;top:0px;left:-5px;padding-right:5px;background-color:lightgrey;
	width:100%;height:60px;z-index:10;}
	header h1 {	position:absolute;top:-20px;left:355px;	}
	header p {position:absolute;top:0px;left:1200px}
	#left_side{position:fixed;float:left;left:0px;top:60px;width:350px;height:100%;
	background-color:lightgrey;}
	#left_side h2 {margin-top:8px;}
	#left_side h2,h3 {margin-left:15px;}
	#left_side p{margin-left:25px;}	
	.disp_windows{position:relative;left:360px;top:115px;height:450px;width:1000px;}
	.hidden_tab {display:none;} 
	nav{position:fixed;left:350px;top:60px;width:100%;z-index:9;}
	ul {list-style-type: none;margin: 0;padding: 0;	overflow: hidden;background-color: #333;}
	li {float: left;}
	li a {display: block;color: white;text-align: center;padding: 14px 16px;text-decoration: none;}
	li a:hover:not(.active) {background-color: #4CAF50;}	
	.active {background-color: #4CAF50;}
	</style>

	<body>
	<header>
	<h1>LADIMIR analysis report</h1>
	<p id=\"date\">$dt</p>
	</header>
	<nav id=\"nav_bar\" role=\"navigation\">
		<ul>
		<li id=\"one\" class=\"active\"><a href=\"#graph\">Colocalization</a></li>
		<li id=\"three\" ><a href=\"#mir_in\">Table 1: miRNA Within cLAD</a></li>
		<li id=\"four\"><a href=\"#mir_out\">Table 2: miRNA outside cLAD</a></li>
		<li id=\"five\"><a href=\"#genes\">miRNA Geneomic Context</a></li>
		</ul>
	</nav>
	
	<div id=\"graph\" class=\"disp_windows\" >
		<div id=\"container\" style=\"min-width: 600px; max-width: auto; height: 600px; width: 1200px; margin: 0 0\"></div>
	</div>
		
	<div id=\"mir_in\" class=\"disp_windows hidden_tab\" >
		<table id=\"Table\" class=\"display\" width=\"80%\"></table>
	</div>
	
	<div id=\"mir_out\" class=\"disp_windows hidden_tab\">
		<table id=\"Table2\" class=\"display\" width=\"80%\"></table>
	</div>
	
	<div id=\"genes\" class=\"disp_windows hidden_tab\">
		<table id=\"Table3\" class=\"display\" width=\"80%\"></table>
	</div>

	
	<div id=\"left_side\" class=\"\">
		<div id=\"results\"> 
			<h2>Results from analysis</h2>
			<h3>Species</h3>
			<p>$species</p>
			<h3>miRNA analysis</h3>
			<p>Nr of miR's within LADS: $mir_wi_number </br>
			Nr of miR's outside of LADS : $mir_ot_number</br>
			Total numbers of pre-miRNA: $mir_all_number</p>	
			<p>% of miRNA outside LADs: $precent%</p>
			<h3>Gene file analysis</h3>
			<p>Nr of genes found in file : $genes_number</p>
			<p>Maximum distance betwen miRNA and Genes: $cutoff bp</p>
			<p>Nr of Genes found with cutoff :$genereport </p>
			<h3></h3>
			<p></p>
		</div>
	</div>
  </body>
  <script>

 \$(function () {
    \$('#container').highcharts({

		    exporting: {
            sourceWidth: 1200,
            sourceHeight: 600,
            // scale: 2 (default)
            chartOptions: {
                subtitle: null
            }
        },
        title: {
            text: 'Number of miRNAs within , outside of cLADS and gene density.'
        },
        xAxis: {
            categories: [$print_out_ot_name],
            title: {
                text: null
            }
        },
        yAxis: [{ // Primary yAxis
			
			labels: {
                format: '{value} miRNAs',
                style: {
                    color: Highcharts.getOptions().colors[1]
                }
            },

        }, { // Secondary yAxis
			
            labels: {
                format: '{value} Genes/Mbp',
                style: {
                    color: Highcharts.getOptions().colors[0]
                }
            },
            opposite: true
        }],
       
        tooltip: {
            shared : false
        },
        plotOptions: {
            bar: {
                dataLabels: {
                    enabled: true
                }
            }
        },
        legend: {
            layout: 'vertical',
            align: 'right',
            verticalAlign: 'top',
            x: -120,
            y: 0,
            floating: true,
            borderWidth: 1,
            backgroundColor: ((Highcharts.theme && Highcharts.theme.legendBackgroundColor) || '#FFFFFF'),
            shadow: true
        },
        credits: {
            enabled: false
        },
        series: [{
            name: 'Within cLAD',
			type: 'column',
			data: [$print_out_wi_numb]            
        },{
            name: 'Outside cLAD',
			type: 'column',
			data: [$print_out_ot_numb]
           
        },{
            name: 'Gene Density',
            type: 'column',
			yAxis: 1,
            data: [$print_genedensity_numbers]        
		}]
    });
});
var dataSet = [
	$Datatable_info
];
\$(document).ready(function() {
    \$('#Table').DataTable( {
        data: dataSet,
        columns: [
            { title: \"Name\" },
            { title: \"Start position\" },
            { title: \"End position\" },
            { title: \"chromosome\" }
        ]
    } );
} );
var dataSet2 = [
	$Datatable_info_ot
];
\$(document).ready(function() {
    \$('#Table2').DataTable( {
        data: dataSet2,
        columns: [
            { title: \"Name\" },
            { title: \"Start position\" },
            { title: \"End position\" },
            { title: \"chromosome\" }
        ]
    } );
} );
var dataSet3 =[
	$print_genes
];
\$(document).ready(function() {
    \$('#Table3').DataTable( {
        data: dataSet3,
        columns: [
			{ title: \"miRNA Name\" },
            { title: \"Downstream Gene\" },
            { title: \"Distance from miRNA\" },
            { title: \"Upstream Gene\" },
			{ title: \"Distance from miRNA\" },
			{ title: \"chromosome\" },
			{ title: \"LAD localazation\" }
            
        ]
    } );
} );
\$( document ).ready(function() {
     \$(window).on(\"hashchange\", function(){
         var hash = window.location.hash.substring(1); // hash part of url without the first letter (#)
		switch(hash){
			case \"graph\":
				\$(\"#\"+hash).removeClass('hidden_tab'); 
				\$(\"#mir_in\").addClass('hidden_tab');
				\$(\"#mir_out\").addClass('hidden_tab');
				\$(\"#genes\").addClass('hidden_tab');
				
				\$(\"#one\").addClass('active');
				\$(\"#three\").removeClass('active');
				\$(\"#four\").removeClass('active');
				\$(\"#five\").removeClass('active');
			break;
			
			case \"mir_in\":
				\$(\"#\"+hash).removeClass('hidden_tab'); 
				\$(\"#graph\").addClass('hidden_tab');
				\$(\"#mir_out\").addClass('hidden_tab');
				\$(\"#genes\").addClass('hidden_tab');

				\$(\"#three\").addClass('active');
				\$(\"#one\").removeClass('active');
				\$(\"#four\").removeClass('active');
				\$(\"#five\").removeClass('active');
			break;
			
			case \"mir_out\":
				\$(\"#\"+hash).removeClass('hidden_tab'); 
				\$(\"#graph\").addClass('hidden_tab');
				\$(\"#mir_in\").addClass('hidden_tab');
				\$(\"#genes\").addClass('hidden_tab');

				\$(\"#four\").addClass('active');
				\$(\"#one\").removeClass('active');
				\$(\"#three\").removeClass('active');
				\$(\"#five\").removeClass('active');
			break;
			case \"genes\":
				\$(\"#\"+hash).removeClass('hidden_tab'); 
				\$(\"#graph\").addClass('hidden_tab');
				\$(\"#mir_out\").addClass('hidden_tab');
				\$(\"#mir_in\").addClass('hidden_tab');

				
				\$(\"#five\").addClass('active');
				\$(\"#one\").removeClass('active');
				\$(\"#three\").removeClass('active');
				\$(\"#four\").removeClass('active');
			break;
			
		}
	});
}); 
</script>
  
</html>  
";

close $LADFILE;
close $OUT;
close $OUT2;
close $OUT3; ######REMOVE LATER
}
