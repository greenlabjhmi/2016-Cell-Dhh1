# Given a UCSC Genome Browser Table Output and a genome
# generate the sequence of all listed transcripts.

from string import maketrans
import math
import numpy as np

class Gene(object):
	
	def __init__(self, InString):
		Split = InString.split('\t')
		
		self.name      = Split[1]
		self.chr       = Split[2]
		self.strand    = Split[3]
		self.numexon   = int(Split[8])
		self.exstart   = [int(x) for x in Split[9].split(',')[:-1]]
		self.exstop    = [int(x) for x in Split[10].split(',')[:-1]]
		self.exframe   = [int(x) for x in Split[15].split(',')[:-1]]
		
	def DefineSeq(self,InString):
		Sequence      = ''
		for Idx in xrange(0,len(self.exstart)):
	        	if self.strand == '+':
                		Sequence += InString[self.exstart[Idx]:self.exstop[Idx]]
            		else:
                		Sequence += InString[self.exstart[len(self.exstart) - Idx - 1]:self.exstop[len(self.exstart) - Idx - 1]][::-1].translate(maketrans("TAGC", "ATCG"))
		
		self.length    = len(Sequence)	
		self.sequence  = Sequence
		self.codonlist = map(''.join, zip(*[iter(Sequence)]*3))
		
	def ProfilingDensity(self, DensityDict, Samples):
		DensityWrapper = {}
		RPKMWrapper    = {}
		TotalWrapper   = {}
		
		for Idx in Samples:
			Density = np.array([])	
			for IdxN in xrange(0,len(self.exstart)):
				if self.strand == '+':
					Density = np.concatenate((Density, DensityDict[Idx][0][self.chr][self.exstart[IdxN]:self.exstop[IdxN]]))
				else:
					Density = np.concatenate((Density, DensityDict[Idx][1][self.chr][self.exstart[len(self.exstart) - IdxN - 1]:self.exstop[len(self.exstart) - IdxN - 1]][::-1]))
					
			DensityWrapper[Idx] = Density
			RPKMWrapper[Idx]    = np.sum(Density)/self.length*1000
			TotalWrapper[Idx]   = np.sum(Density)
			
		self.profdensity = DensityWrapper
		self.profRPKM    = RPKMWrapper	
		self.profTotal   = TotalWrapper
	
	def RNASeqDensity(self, DensityDict, Samples):
		DensityWrapper = {}
		RPKMWrapper    = {}
		TotalWrapper   = {}

		for Idx in Samples:
			Density = np.array([])	
			for IdxN in xrange(0,len(self.exstart)):
				if self.strand == '+':
					Density = np.concatenate((Density, DensityDict[Idx][0][self.chr][self.exstart[IdxN]:self.exstop[IdxN]]))
				else:
					Density = np.concatenate((Density, DensityDict[Idx][1][self.chr][self.exstart[len(self.exstart) - IdxN - 1]:self.exstop[len(self.exstart) - IdxN - 1]][::-1]))

			DensityWrapper[Idx] = Density
			RPKMWrapper[Idx]    = np.sum(Density)/self.length*1000
			TotalWrapper[Idx]   = np.sum(Density)

		self.mRNAdensity = DensityWrapper
		self.mRNARPKM    = RPKMWrapper	
		self.mRNATotal   = TotalWrapper	
				
	def CalculateTAI(self, TAIDict):
		TAI = 0.0

		for iC in self.codonlist[:-1]:
			if iC in TAIDict:
				TAI += math.log(TAIDict[iC],2)

		self.TAI = 2**(TAI/len(self.codonlist[:-1]))