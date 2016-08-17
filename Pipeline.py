import subprocess as sp
import time, sys
import gzip

def cutAdapt(IdxIn,StartTime):
	
	# Define the command for Popen

	CutAdapt     = ["cutadapt","-a","CTGTAGGCACCATCAATAGATCGGAA","-o","2-CutAdapt/","1-Data/"]
	CutAdapt[4] += IdxIn + ".fastq"
	CutAdapt[5] += IdxIn + ".fastq.gz"

	# Perform cutadapt adapter removal

	Process = sp.Popen(CutAdapt, stdout=sp.PIPE, stderr=sp.PIPE)

	# Take the output of the CutAdapt call and print it to file!

	FileOut = open("2-CutAdapt/Reports/" + IdxIn + ".txt","w")
	FileOut.write(Process.communicate()[0])
	FileOut.close()

	print IdxIn, " - CutAdapt - ", time.time() - StartTime

def qualFilter(IdxIn,StartTime):

	PHREDDict = {
		"!":9.999999e-01, "\"":7.943282e-01, "#":6.309573e-01, "$":5.011872e-01, "%":3.981072e-01,
		"&":3.162278e-01, "\'":2.511886e-01, "(":1.995262e-01, ")":1.584893e-01, "*":1.258925e-01,
		"+":1.000000e-01, ",":7.943282e-02, "-":6.309573e-02, ".":5.011872e-02, "/":3.981072e-02,
		"0":3.162278e-02, "1":2.511886e-02, "2":1.995262e-02, "3":1.584893e-02, "4":1.258925e-02,
		"5":1.000000e-02, "6":7.943282e-03, "7":6.309573e-03, "8":5.011872e-03, "9":3.981072e-03,
		":":3.162278e-03, ";":2.511886e-03, "<":1.995262e-03, "=":1.584893e-03, ">":1.258925e-03,
		"?":1.000000e-03, "@":7.943282e-04, "A":6.309573e-04, "B":5.011872e-04, "C":3.981072e-04,
		"D":3.162278e-04, "E":2.511886e-04, "F":1.995262e-04, "G":1.584893e-04, "H":1.258925e-04,
		"I":1.000000e-04, "J":7.943282e-05
	}

	# Get number of lines and divide by four to get the number of individual entries. Then, calculate
	# the quality score for each line and trim everything with less than a 99.5% confidence in sequence.
	# The confidence is defined as the product of (1 - the quality scores of the individual nucleotides).
	
	File    = open("2-CutAdapt/Reports/" + IdxIn + ".txt")

	for Idx in xrange(0,5):
		Line = File.readline()

	File.close()

	Length  = int(Line.rstrip("\n").split(" ")[-1])

	File    = open("2-CutAdapt/" + IdxIn + ".fastq")
	FileOut = open("3-Filtered/" + IdxIn + ".fastq", "w")
	
	for IdxN in range(0,Length):
		Identifier  = File.readline().rstrip("\n")
		Sequence    = File.readline().rstrip("\n")
		QIdentifier = File.readline().rstrip("\n")
		PHRED       = File.readline().rstrip("\n")
		Score       = 1.0
		Len         = len(PHRED)
		
		if Len < 4:
			pass
		else:
			for IdxL in range(0,Len):
				Score   = Score*(1 - PHREDDict[PHRED[IdxL]])	
			
			if (Score > 0.99):
				FileOut.write(Identifier + "\n" + Sequence + "\n" + QIdentifier + "\n" + PHRED + "\n")
			
	File.close()
	FileOut.close()
	
	print IdxIn, " - QualFilt - ", time.time() - StartTime

def ncRNAAlign(IdxIn,StartTime):

	# Use Bowtie to align to the non-coding RNAs/tRNAs...
	# Modify line 86 to use different settings as documented in the Bowtie manual

	# The bowtie settings being used are...
	# -S to output as a SAM file
	# -v to output as 
	# -p 4 to use all four cores on my processor
	# --best to ensure that the reported alignments are the best in terms of mismatches
	# --un to generate the FASTQ file of all the unmapped reads (i. e. those that map to the genome)
	
	Input    = "3-Filtered/" + IdxIn + ".fastq"
	Output   = "4-NC-Subtracted/SAM/" + IdxIn + ".sam"
	Unmapped = "4-NC-Subtracted/" + IdxIn + ".fastq"
	Bowtie  = ["bowtie", "-Sv", "3", "-p", "4", "--best", "--un", Unmapped, "0-Reference/Indexes/ncRNA", Input, Output]
	Process = sp.Popen(Bowtie, stdout=sp.PIPE, stderr=sp.PIPE)

	FileOut = open("4-NC-Subtracted/Reports-" + IdxIn + ".txt","w")
	FileOut.write(Process.communicate()[1])
	FileOut.close()

	print IdxIn, " - BowtieNC - ", time.time() - StartTime

def genomeAlign(IdxIn,StartTime):

	# Take the the processed FASTQ with reads mapping to ncRNA removed and map to the genome!
	
	# The bowtie settings being used are...
	# -S to output as a SAM file
	# -m to output only reads that map unambiguously
	# -p 4 to use all four cores on my processor
	# --best to ensure that the reported alignments are the best in terms of mismatches
	# --strata to ensure that the reported alignments are only in the best strata

	Input   = "4-NC-Subtracted/" + IdxIn + ".fastq"
	Output  = "5-Aligned/" + IdxIn + ".sam"
	Bowtie  = ["bowtie","-Sm","1","-p","4","--best","--strata","0-Reference/Indexes/Hybrid", Input, Output]
	Process = sp.Popen(Bowtie, stdout=sp.PIPE, stderr=sp.PIPE)
	
	FileOut = open("5-Aligned/Reports/" + IdxIn + ".txt","w")
	FileOut.write(Process.communicate()[1])
	FileOut.close()
	
	# sp.call(["gzip","-f","2-CutAdapt/" + IdxIn + ".fastq"])
	# sp.call(["gzip","-f","3-Filtered/" + IdxIn + ".fastq"])
	# sp.call(["gzip","-f","4-NC-Subtracted/" + IdxIn + ".fastq"])
	# sp.call(["gzip","-f","4-NC-Subtracted/SAM/" + IdxIn + ".sam"])
	
	print Idx, " - BowtieGe - ", time.time() - StartTime

def makeWIGs(IdxIn,StartTime):

	# Take a given SAM file Bowtie output and use it to generate a WIG file.
	
	# If you only want to generate WIG files of a specific size or want to use 3' or 5' end mapping of the positions
	# this is where you have to modify the code accordingly. Remember that the output of Bowtie leads to mapping of
	# the read to the nucleotide left of where mapping begins. This is true for both strands!

	LengthDict = {
		  "chrI":230218,   "chrII":813184,  "chrIII":316620, "chrIV":1531933,   "chrV":576874,  
		 "chrVI":270161, "chrVII":1090940, "chrVIII":562643,  "chrIX":439888,   "chrX":745751, 
		 "chrXI":666816, "chrXII":1078177, "chrXIII":924431, "chrXIV":784333, "chrXV":1091291,
		"chrXVI":948066,     "chrM":85779, "pAG425-NL-HA-Dhh1": 9470, "pYES52GAL-mCherry":6911
	}


	# Chose your mapping for the WIG file and apply and length filters and set the counter of
	# mapped nucleotides to zero. Finally, read through the 19 lines of header in a SAM file 
	# that occur when mapping to the yeast genome.

	Mapping = "3"    # "3", "Center"
	Lengths = "All"  # 28, 21, "All"
	Count   = 0      # Number of mapped nucleotides

	File       = open("5-Aligned/" + IdxIn + ".sam")
	
	# Create dictionaries for the plus and minus strands, each with a different chromosome
	# as a key which links to a list of '0' values of the length of the chromosome.
	
	PlusDict   = {}
	MinusDict  = {}
	
	Burn    = [File.readline() for line in xrange(21)] 

	for iL in LengthDict:
		PlusDict[iL]  = [0]*LengthDict[iL]
		MinusDict[iL] = [0]*LengthDict[iL]

	# Read through the file and assign density based on the above selection criterion.

	for iL in File:
		Split = iL.split("\t")
		Len   = len(Split[9])

		if (Len == Lengths) or (Lengths == "All"):
			if Mapping == "5":
				if Split[1] == "0":
					PlusDict[Split[2]][int(Split[3]) - 1]                  += 1
					Count                                                  += 1
				if Split[1] == "16":
					MinusDict[Split[2]][int(Split[3]) + len(Split[9]) - 2] += 1
					Count                                                  += 1
			if Mapping == "3":
				if Split[1] == "0":
					PlusDict[Split[2]][int(Split[3]) + len(Split[9]) - 2]  += 1
					Count                                                  += 1
				if Split[1] == "16":
					MinusDict[Split[2]][int(Split[3]) - 1]                 += 1
					Count                                                  += 1
			if Mapping == "Center":
				if Split[1] == "0":
					for Idx in xrange(len(Split[9])):
						PlusDict[Split[2]][Idx + int(Split[3]) - 1]        += 1
						Count                                              += 1
				if Split[1] == "16":
					for Idx in xrange(len(Split[9])):
						MinusDict[Split[2]][Idx + int(Split[3]) - 1]       += 1
						Count                                              += 1
		else:
			pass

	print Count
	
	FileOutP   = open("6-WIGs/" + Mapping +  "-End/" + str(Lengths) + "/RPM/" + IdxIn + "-P.wig","w")
	FileOutM   = open("6-WIGs/" + Mapping +  "-End/" + str(Lengths) + "/RPM/" + IdxIn + "-M.wig","w")
	FileOutPC  = open("6-WIGs/" + Mapping +  "-End/" + str(Lengths) + "/Counts/" + IdxIn + "-P.wig","w")
	FileOutMC  = open("6-WIGs/" + Mapping +  "-End/" + str(Lengths) + "/Counts/" + IdxIn + "-M.wig","w")

	FileOutP.write("track name=tracklabel totalCounts=" + str(Count) + " viewLimits=-5:5 color=79,159,36\n")
	FileOutM.write("track name=tracklabel totalCounts=" + str(Count) + " viewLimits=-5:5 color=79,159,36\n")

	FileOutPC.write("track name=tracklabel totalCounts=" + str(Count) + " viewLimits=-5:5 color=79,159,36\n")
	FileOutMC.write("track name=tracklabel totalCounts=" + str(Count) + " viewLimits=-5:5 color=79,159,36\n")
	
	for IdxC in PlusDict:
		FileOutP.write("fixedStep  chrom=" + IdxC + "  start=1  step=1\n")
		FileOutM.write("fixedStep  chrom=" + IdxC + "  start=1  step=1\n")
		FileOutPC.write("fixedStep  chrom=" + IdxC + "  start=1  step=1\n")
		FileOutMC.write("fixedStep  chrom=" + IdxC + "  start=1  step=1\n")
		
		for IdxN in range(0, len(PlusDict[IdxC])):
			FileOutP.write(str(PlusDict[IdxC][IdxN]*1000000.0/Count) + "\n")
			FileOutM.write(str(MinusDict[IdxC][IdxN]*1000000.0/Count) + "\n")
			FileOutPC.write(str(PlusDict[IdxC][IdxN]) + "\n")
			FileOutMC.write(str(MinusDict[IdxC][IdxN]) + "\n")
	
	FileOutP.close()	
	FileOutM.close()
	FileOutPC.close()	
	FileOutMC.close()
	
	print IdxIn, " - WIGsMade - ", time.time() - StartTime	

def errorMessage(In):

	# Spits out an informative error message if the command line call is not of the right format.

	if In == 0:
		print "ERROR: You need to enter what step you want to start processing from!"
		print ""
		print "0 - From the original FASTQ files"
		print "1 - From the CutAdapt filtered FASTQ files"
		print "2 - From the Quality filtered FASTQ files"
		print "3 - From the ncRNA aligned SAM files"
		print "4 - From the ncRNA subtracted FASTQ files"
		print "5 - From the genome aligned SAM files"
		print ""
		print "In addition, provide a directory where the input files are found (no prefixes!)"

		sys.exit(0)	
	if In == 1:
		print "ERROR: No files found at that location!"

		sys.exit(0)

#################################################################################################################
#                                                                                                               #
#                    Here begins the non-function (but, still functional) part of the code                      #
#                                                                                                               #
#################################################################################################################

Start = time.time()

# Check to see that the input is of the proper format, if not, use errorMessage() to handle it.

if len(sys.argv) != 3:
	errorMessage(0)
else:
	try:
		Step = int(sys.argv[1])
	except ValueError:
		errorMessage(0)
		
	if (Step < 0) or (Step > 4):
		errorMessage(0)

	# Define the files over which to operate

	LS = ["ls", "-1", sys.argv[2]]
	Process  = sp.Popen(LS, stdout=sp.PIPE, stderr=sp.PIPE)
	List     = Process.communicate()[0].split("\n")
	List     = filter(None, List)

	if len(List) == 0:
		errorMessage(1)

	Files = [Idx.split(".")[0] for Idx in List]

# You can override the input from the command line by specifying the values below!
# This is great for debugging in case you want to make sure that you're not going to nuke
# you're previously calculated hard work!

# Files    = ["Ini-Prof-TetherM"] 
# Step     = 4

for Idx in Files:
	Options  = {0: cutAdapt, 1: qualFilter, 2: ncRNAAlign, 3: genomeAlign, 4: makeWIGs}

	for OptIdx in range(Step, 5):
		Options[OptIdx](Idx, Start)