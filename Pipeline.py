import subprocess as sp
import shutil
import gzip
import time

# This code requires the following programs to be installed on your computer
#
# 1) sra-tools (https://github.com/ncbi/sra-tools/wiki/Downloads)
# 2) cutadapt  (https://cutadapt.readthedocs.io/en/stable/)

def cleanFile(File, Condition):

	# gZip a file and delete the un-gZipped version!

	if Condition == "gzip":
		with open(File, 'rb') as FIn, gzip.open(File + ".gz", 'wb') as FOut:
			shutil.copyfileobj(FIn, FOut)

		sp.Popen(["rm", File])

	if Condition == "rm":
		sp.Popen(["rm", "4-Subtracted/" + iX + ".fastq"])

def makeDirectory(Path):
	
	# Check if a folder named exists at Path. If not, create it!

	Split = Path.split("/")

	if len(Split) == 1:
		List = sp.Popen(["ls"], stdout=sp.PIPE, stderr=sp.PIPE) 
	else:
		List = sp.Popen(["ls", "/".join(Split[:-1])], stdout=sp.PIPE, stderr=sp.PIPE)
	
	if Path in List.communicate()[0].split("\n"):
		pass
	else:
		Make = sp.Popen(["mkdir", Path])

def parseParams(Path):

	# Open the parameter file and read in parameters and files to operate on

	File      = open(Path)
	SRAList   = []
	NameList  = []
	ParamDict = {}

	for iL in File:
		if iL[0] == "#":
			Split = iL[2:-1].split(":")
			if len(Split) > 1:
				ParamDict[Split[0].strip()] = Split[1].strip()
		else:
			if len(iL) != 1:
				Split = iL[:-1].split("\t")
				SRAList.append(Split[0])
				NameList.append(Split[1])

	return ParamDict, SRAList, NameList

def downloadData(SRAList, NameList, Params):

	# Check to see if all the files 

	makeDirectory("1-Raw")

	for iX in xrange(len(SRAList)):
		Dump  = sp.Popen(["fastq-dump", SRAList[iX]], stdout=sp.PIPE, stderr=sp.PIPE)
		Dump.wait()

		Move  = sp.Popen(["mv", SRAList[iX] + ".fastq", "1-Raw/" + NameList[iX] + ".fastq"])
		Move.wait()

		cleanFile("1-Raw/" + NameList[iX] + ".fastq", Params["Clean"])

def cutAdapt(SRAList, NameList, Params):

	makeDirectory("2-Trimmed")
	makeDirectory("2-Trimmed/Reports")

	for iX in NameList:
		CutAdapt     = ["cutadapt","-a","CTGTAGGCACCATCAATAGATCGGAA","-o","2-Trimmed/","1-Raw/"]
		CutAdapt[4] += iX + ".fastq.gz"
		CutAdapt[5] += iX + ".fastq.gz"

		Trim         = sp.Popen(CutAdapt, stdout=sp.PIPE, stderr=sp.PIPE)
		Trim.wait()

		FileOut = open("2-Trimmed/Reports/" + iX + ".txt","w")
		FileOut.write(Trim.communicate()[0])
		FileOut.close()

def qualityFilter(SRAList, NameList, Params):

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

	makeDirectory("3-Filtered")

	for iX in NameList:
		File   = open("2-Trimmed/Reports/" + iX + ".txt")
		Burn   = [File.readline() for Idx in xrange(0,8)]
		Length = int(Burn[-1][:-1].split(" ")[-1].replace(",",""))
		File.close()

		File    = gzip.open("2-Trimmed/" + iX + ".fastq.gz", "rb")
		FileOut = open("3-Filtered/" + iX + ".fastq", "w")
	
		for iN in range(0,Length):
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

def ncRNASubtract(SRAList, NameList, Params):

	makeDirectory("4-Subtracted")
	makeDirectory("4-Subtracted/SAM")
	makeDirectory("4-Subtracted/Reports")

	for iX in NameList:
		Input    = "3-Filtered/" + iX + ".fastq"
		Output   = "4-Subtracted/SAM/" + iX + ".sam"
		Unmapped = "4-Subtracted/" + iX + ".fastq"
		Bowtie   = ["bowtie", "-Sv", "3", "-p", "4", "--best", "--un", Unmapped, "0-Reference/Indexes/ncRNA", Input, Output]
		Subtract = sp.Popen(Bowtie, stdout=sp.PIPE, stderr=sp.PIPE)
		Subtract.wait()

		FileOut  = open("4-Subtracted/Reports/" + iX+ ".txt","w")
		FileOut.write(Subtract.communicate()[1])
		FileOut.close()

		cleanFile("3-Filtered/" + iX + ".fastq", Params["Clean"])
		cleanFile("4-Subtracted/SAM/" + iX + ".sam", Params["Clean"])

def genomeAlign(SRAList, NameList, Params):

	makeDirectory("5-Aligned")
	makeDirectory("5-Aligned/Reports")

	for iX in NameList:
		Input   = "4-Subtracted/" + iX + ".fastq"
		Output  = "5-Aligned/" + iX + ".sam"
		Bowtie  = ["bowtie","-Sm","1","-p","4","--best","--strata","0-Reference/Indexes/Genome", Input, Output]
		Align   = sp.Popen(Bowtie, stdout=sp.PIPE, stderr=sp.PIPE)

		FileOut = open("5-Aligned/Reports/" + iX + ".txt", "w")
		FileOut.write(Align.communicate()[1])
		FileOut.close()

		cleanFile("4-Subtracted/" + iX + ".fastq", Params["Clean"])

def makeWIGs(SRAList, NameList, Params):

	Mapping = Params["Mapping"]    # "3", "Center"

	if Params["Lengths"] == "All":
		Lengths = "All"  # 28, 21, "All"
	else:
		Lengths = int(Params["Lengths"])

	makeDirectory("6-WIGs")
	makeDirectory("6-WIGs/" + Mapping + "-End")
	makeDirectory("6-WIGs/" + Mapping + "-End/" + str(Lengths))

	for iX in NameList:

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
		
		File    = open("5-Aligned/" + iX + ".sam")
		Count   = 0      # Number of mapped nucleotides

		
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
				if Mapping == "c":
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
		
		FileOutP   = open("6-WIGs/" + Mapping +  "-End/" + str(Lengths) + "/" + iX + "-P.wig","w")
		FileOutM   = open("6-WIGs/" + Mapping +  "-End/" + str(Lengths) + "/" + iX + "-M.wig","w")

		FileOutP.write("track name=tracklabel totalCounts=" + str(Count) + " viewLimits=-5:5 color=79,159,36\n")
		FileOutM.write("track name=tracklabel totalCounts=" + str(Count) + " viewLimits=-5:5 color=79,159,36\n")

		for IdxC in PlusDict:
			FileOutP.write("fixedStep  chrom=" + IdxC + "  start=1  step=1\n")
			FileOutM.write("fixedStep  chrom=" + IdxC + "  start=1  step=1\n")
			
			for IdxN in range(0, len(PlusDict[IdxC])):
				FileOutP.write(str(PlusDict[IdxC][IdxN]*1000000.0/Count) + "\n")
				FileOutM.write(str(MinusDict[IdxC][IdxN]*1000000.0/Count) + "\n")
		
		FileOutP.close()	
		FileOutM.close()

Params, SRAs, Names = parseParams("Param.in")
Options             = {1: downloadData, 2: cutAdapt, 3: qualityFilter, 4: ncRNASubtract, 5: genomeAlign, 6: makeWIGs}
Start               = time.time()

for iOpt in range(int(Params["Start"]), int(Params["Stop"]) + 1):
	Options[iOpt](SRAs, Names, Params)

	print "Step " + str(iOpt) + " completed! Time taken thus far: " + str(time.time() - Start)
