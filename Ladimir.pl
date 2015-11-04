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
# usage: Ladimir_v2.pl [LAD file] [miRNA file] "name of animal"

use strict;
use Data::Dumper;
use Array::Utils qw(:all);

## Read in the data
my $lad_file = $ARGV[0];
my $mir_file = $ARGV[1];
my $animal = $ARGV[2];
my $outfile = "$animal.txt";
my $outfile2 = "$animal.new.gff3";

if ($ARGV[0] eq "-h"){

print <<EOF;
Ladimir is used to search for micro RNA's within cLAD regions.
USAGE:
Ladimir_v2.pl [cLAD.bed] [miR.gff3] [name of outputfiles]

For the [name of outputfiles] no file extentions are needed.
Ladimir will create a file containing miRNA's that are within cLADS
And a .gff3 file with a list of miRNA that where used to do the comparison. 

EXAMPLE USAGE:
~perl Ladimir_v2.pl cLAD_mm10.bed mmu10.gff3 mouse_new

EOF
exit;

}else{
## Variables
my @input_lad;
my @input_mir;

## Open the file's
open ( my $LADFILE, '<', $lad_file )
	or die "No such file: $lad_file";
	
open ( my $MIRFILE, '<', $mir_file )
	or die "No such file: $mir_file";
	
	
## File to write out
open my $OUT, '>', $outfile or die "Cannot open: $outfile";
open my $OUT2, '>', $outfile2 or die "Cannot open: $outfile2";

print "processing LAD information\n";
## Start parsing the files. #Parsin LAD files.
while ( my $line = <$LADFILE> ){
	chomp $line;
	my ($chr, $start, $end) = split (/\s/, $line);
	push @input_lad, [$chr, $start, $end];
}

my @indata = <$MIRFILE>;
close $MIRFILE;

my $selection;
# print $OUT @indata;

print "processing miR information\n";
if( (@indata[3]) =~/.+Mus musculus.+/ ) {
   print "File is from Mus musculus\n";
   $selection=1;
    foreach  my $line (@indata ){
		chomp $line;
		if ($line =~ /(\d{1,2}|\w)\t.\tmiRNA_primary_transcript/){
			print $OUT2 $line,"\n";
			$line =~ /(\d{1,2}|\w)\t.\tmiRNA_primary_transcript\t(\d+)\s(\d+).+Name=(.+\b)/;
			push @input_mir, ["chr$1", $2, $3, $4];
		}
	}
}
if( (@indata[3]) =~/.+Drosophila melanogaster.+/ ){
   print "File is from Drosophila M.\n";
   $selection=1;
    foreach  my $line (@indata ){
		chomp $line;
		if ($line =~ /(.{1,5})\t.\tmiRNA_primary_transcript/){
			print $OUT2 $line,"\n";
			$line =~ /(.{1,5})\t.\tmiRNA_primary_transcript\t(\d+)\s(\d+).+Name=(.+\b)/;
			push @input_mir, ["chr$1", $2, $3, $4];
		}
	}
}
if( (@indata[3]) =~/.+Homo sapiens.+/ ){
   print "File is from Human\n";
   $selection=1;
   foreach  my $line (@indata ){
		chomp $line;
		if ($line =~ /(.{1,6})\t.\tmiRNA_primary_transcript/){
			print $OUT2 $line,"\n";
			$line =~ /(.{1,6})\t.\tmiRNA_primary_transcript\t(\d+)\s(\d+).+Name=(.+\b)/;
			push @input_mir, [$1, $2, $3, $4];
		}
	}
}
if( (@indata[3]) =~/.+Caenorhabditis elegans.+/ ){
   print "File is from C.Elegans\n";
   $selection = 2;
    foreach  my $line (@indata ){
		chomp $line;
		if ($line =~ /(.{1,6})\t.\tmiRNA_primary_transcript/){
			print $OUT2 $line,"\n";
			$line =~ /(.{1,6})\t.\tmiRNA_primary_transcript\t(\d+)\s(\d+).+Name=(.+\b)/;
			push @input_mir, [$1, $2, $3, $4];
		}
	}
}

#Start the comparison algorithm
my $lad_scalar = scalar @input_lad;
my $mir_scalar = scalar @input_mir;
my @mir_os_lad;
my @mir_wi_lad;

# Itterate through all the miR's
for ( my $i = 0;$i < $mir_scalar;$i++){
	# Itterate through all the LADs
	for (my $j = 0;$j < $lad_scalar;$j++){
		my $temp_length = $input_lad[$j][2] - $input_lad[$j][1]; #Find current LADs Length
	#Check if chromosome is the same nr.
	if ($input_lad[$j][0] eq $input_mir[$i][0] ){
		#Build Multi dimentional array that has the information on miR's inside of LADS.
		if( ( ($input_mir[$i][1] - $input_lad[$j][1]) >= 0 ) and ( ($input_mir[$i][1] - $input_lad[$j][1]) <= $temp_length ) ){
			push @mir_wi_lad,[ $input_mir[$i][0],$input_mir[$i][1],$input_mir[$i][2],$input_mir[$i][3] ];
			$i++;
			}
		}
	}
	
	#printing the percentage of the progress
	my $per = int(($i/$mir_scalar)*100);
	print "\b\b\b\b$per%";
}
print "\b\b\b\bComputation is complete\n";
# Create a array's from known mir inside lad containing only names
my @mir_wi_lad_name;
my @mir_wi_lad_chr;
for my $i (@mir_wi_lad){
	push @mir_wi_lad_name, @$i[3];
}
my @mir_all_lad_name;
for my $i (@input_mir){
	push @mir_all_lad_name, @$i[3];
}

#Check to see which miRNA are on the outside or between catagory
my @mir_os_lad_new = array_minus (@mir_all_lad_name,@mir_wi_lad_name);

#Create an array which contains the miR's that are out side of LAD
my @mir_os_lad_temp;
foreach my $line (@input_mir){
	for (my $i=0; $i < scalar @mir_os_lad_new; $i++){
		if ( @$line[3] eq $mir_os_lad_new[$i] ){	
			my $name = @$line[3];
			my $pst = @$line[1];
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

print "This are miR's outside of LADS ", scalar @mir_os_lad,"\n";
print $OUT "This are miR's outside of LADS ", scalar @mir_os_lad,"\n";

print "This is the nunmber of miR's within LADS ", scalar @mir_wi_lad ,"\n";
print $OUT "This is the nunmber of miR's within LADS ", scalar @mir_wi_lad ,"\n";

print "This is the number of ALL miR's ", scalar @mir_all_lad_name, "\n";
print $OUT "This is the number of ALL miR's ", scalar @mir_all_lad_name, "\n";

##Here I will start the stastistical analysis
my %mir_count;
my %mir_count_out;
my @mir_inside_stat;
my @mir_outside_stat;
my %mir_name_stat;

foreach my $line (@mir_wi_lad){
	##Here I check to see the number of miR within
	$mir_count{$$line[0]}++;
	# push @mir_inside_stat, $$line[0];
}
foreach my $line (@mir_os_lad){
	$mir_count_out{$$line[0]}++;
	# push @mir_outside_stat, $line;
}

# # Do calculation for miRNA within cLADS
print "\nThis tables show how many miRNA are inside cLADs each chr\n";
print $OUT "This tables show how many miRNA are inside cLADs each chr\n";

if ($selection == 1){
	foreach my $k (sort {substr($a, 3) <=> substr($b, 3)} keys %mir_count){
	print $OUT "$k $mir_count{$k}\n";
	# push(@mir_wi_values,$mir_count{$k});
	# push(@mir_wi_keys,$k);
	print "$k = $mir_count{$k}\n";
	}
}elsif($selection == 2){
	foreach my $k (sort {$a cmp $b} keys %mir_count){
	print $OUT "$k $mir_count{$k}\n";
	# push(@mir_wi_values,$mir_count{$k});
	# push(@mir_wi_keys,$k);
	print "$k = $mir_count{$k}\n";
	}
}

# Do calculation for miRNA outside cLADS
print "\nThis tables show how many miRNA are outside cLADs each chr\n";
print $OUT "\nThis tables show how many miRNA are outside cLADs each chr\n";
if ($selection == 1){
	foreach my $k (sort {substr($a,3) <=> substr($b,3)} keys %mir_count_out){
		print $OUT "$k $mir_count_out{$k}\n";
		# push(@mir_wi_values,$mir_count_out{$k});
		# push(@mir_wi_keys,$k);
		print "$k = $mir_count_out{$k}\n";
	}
}
if ($selection == 2){
	foreach my $k (sort {$a cmp $b} keys %mir_count_out){
		print $OUT "$k $mir_count_out{$k}\n";
		# push(@mir_wi_values,$mir_count_out{$k});
		# push(@mir_wi_keys,$k);
		print "$k = $mir_count_out{$k}\n";
	}
}

print $OUT "\nA list of miR's that Ladimir identidies as within cLADS"; 
print $OUT "\nChr name\tStart pst\tEnd pst\tmiR name\n";
foreach my $i(@mir_wi_lad){
	print $OUT "@$i\n";
}
# print Dumper @indata;
close $LADFILE;
close $OUT;
close $OUT2;
}
