# Ladimir
Ladimir is a script used to search for miRNA that lie within a cLAD boundry.
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
