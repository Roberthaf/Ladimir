Ladimir - A search tool for Lamina acosiated domain and miRNA research
======================================================================

Installation
------------
The software is a stand alone Perl. That requirs Perl version 5.2 or above to run.

Usage
-----
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

