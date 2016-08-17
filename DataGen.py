from Gene import Gene

import numpy as np
import pandas as pd
import subprocess as sp
import cPickle as pickle
import time, sys

Start = time.time()

def LoadDict(inString, LengthDict):
	Density = {0:{},1:{}}
	
	String = ['-P.wig','-M.wig']
	
	for Idx in range(0,2):
		File = open(inString + String[Idx])
	
		for line in File:
			if line[0] == 'f':
				Chr               = line.strip('\n').split('  ')[1][6:]
				if Chr in LengthDict:
					Density[Idx][Chr] = np.array([0.0]*LengthDict[Chr])
					IdxD              = 0
			elif line[0] == 't':
				pass
			else:
				if Chr in LengthDict:
					Density[Idx][Chr][IdxD]  = float(line)
					IdxD                    += 1
				
	return Density

# Load the reference genome

PickleGenome  = open('../0-Reference/Chromosomes.pckl')
Genome        = pickle.load(PickleGenome)
PickleGenome.close()

# Load in the necessary WIG files:

LengthDict = {
      'chrI':230218,   'chrII':813184,  'chrIII':316620, 'chrIV':1531933,   'chrV':576874,  
	 'chrVI':270161, 'chrVII':1090940, 'chrVIII':562643,  'chrIX':439888,   'chrX':745751, 
     'chrXI':666816, 'chrXII':1078177, 'chrXIII':924431, 'chrXIV':784333, 'chrXV':1091291,
    'chrXVI':948066,     'chrM':85779
}

File  = open('../0-Reference/Coding.txt')
Genes = {}
    
for iF in File:
	Name        = iF.split('\t')[1]
	Genes[Name] = Gene(iF)
	Genes[Name].DefineSeq(Genome[Genes[Name].chr])

RNADensity   = {}
ProfDensity  = {}
RNACounts    = {}
ProfCounts   = {}
RNADSamples  = []
ProfDSamples = []
RNACSamples  = []
ProfCSamples = []

LS          = ["ls", "-1", "../5-Aligned/"]
Process     = sp.Popen(LS, stdout=sp.PIPE, stderr=sp.PIPE)
List        = Process.communicate()[0].split("\n")
List        = filter(None, List)

Files = [Idx.split(".")[0] for Idx in List]

Files = ["1-ACAGTG-2781-NALog"]

for iF in Files:
	if iF[4:8] == "Prof":
		ProfDSamples.append("../6-WIGs/3-End/All/RPM/" + iF)
		ProfCSamples.append("../6-WIGs/3-End/All/Counts/" + iF)
	elif iF[4:7] == "RNA":
		RNADSamples.append("../6-WIGs/3-End/All/RPM/" + iF)
		RNACSamples.append("../6-WIGs/3-End/All/Counts/" + iF)

for iR in RNADSamples:
	RNADensity[iR]  = LoadDict(iR, LengthDict)

for iR in RNACSamples:
	RNACounts[iR]   = LoadDict(iR, LengthDict)

for iP in ProfDSamples:
	ProfDensity[iP] = LoadDict(iP, LengthDict)

for iP in ProfCSamples:
	ProfCounts[iP] = LoadDict(iP, LengthDict)

OutD = open('../7-Analysis/Data/RPKM-All-Datasets.txt','w')
OutC = open('../7-Analysis/Data/Counts-All-Datasets.txt','w')

OutD.write("Gene")
OutC.write("Gene")

for iP in ProfDSamples:
	OutD.write("\t" + iP[14] + iP[18] + iP[23:])
for iR in RNADSamples:
	OutD.write("\t" + iR[14] + iR[18] + iR[22:])
for iP in ProfCSamples:
	OutC.write("\t" + iP[17] + iP[21] + iP[26:])
for iR in RNACSamples:
	OutC.write("\t" + iR[17] + iR[21] + iR[25:])

OutD.write("\n")
OutC.write("\n")

for IdxG in Genes:
	if IdxG[0] == 'Y':
		Genes[IdxG].ProfilingDensity(ProfDensity, ProfDSamples)
		Genes[IdxG].RNASeqDensity(RNADensity, RNADSamples)

		OutD.write(IdxG)

		for iP in ProfDSamples:
			OutD.write("\t" + str(Genes[IdxG].profRPKM[iP]))

		for iR in RNADSamples:
			OutD.write("\t" + str(Genes[IdxG].mRNARPKM[iR]))

		OutD.write("\n")
		OutC.write(IdxG)

		Genes[IdxG].ProfilingDensity(ProfCounts, ProfCSamples)
		Genes[IdxG].RNASeqDensity(RNACounts, RNACSamples)

		for iP in ProfCSamples:
			OutC.write("\t" + str(Genes[IdxG].profTotal[iP]))

		for iR in RNACSamples:
			OutC.write("\t" + str(Genes[IdxG].mRNATotal[iR]))

		OutC.write("\n")
