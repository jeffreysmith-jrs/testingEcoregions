#Import libraries
import os,sys,random,csv,sys,time,shutil,fnmatch,scipy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sin, cos, sqrt, log, radians, fabs
from matplotlib import pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from pandas import DataFrame
import ecopy as ep
from collections import Counter
from scipy.optimize import curve_fit

#Start timer
start = time.time()

#Here's a function to remove occurences with too few pixels, I haven't used it in any analysis, but you're welcome to play around with it
def removeMinOccs(data, minOccs):
	data	=	data[np.logical_or.reduce([data[:,-1] >= minOccs])]
	return data

def reduceMaxOccs(data, maxOccs):	
	occupiedPixels	=	np.count_nonzero(data[:,-1])

	
	#Randomize overrepresented pixels
	npix	=	len(data)
	stats	=	np.zeros(shape=(npix, 8))

	#Loop through pixels
	x = 0
	for x in range(npix):

		#Identify pixels with too many occs
		if np.sum(data[x,5:-1]) > maxOccs:
			
			#Extract that line to play with
			subArray	=	np.reshape(data[x,5:-1], len(data[x,5:-1]))
			
			#Make list to hold the randomly sampled occs
			elements	=	[]
			
			#Loop through occs in that pixel
			for z in range(len(subArray)):
				
				#make a list with x entries for each species where x is the number of occs in that species in that pixel
				for s in range(int(subArray[z])):
					toAppend	=	['Species.' + str(z)]
					elements	=	elements + toAppend
			
			#Randomly subsample from the array of species (where the number of occs in that species = the number of occs in that pixel)
			chosen	=	[]
			
			#only choose up to the number of max occs
			for y in range(maxOccs):
				
				#ranomdly select an element
				select	=	random.randint(0,len(elements)-1)
				selected	=	elements[select]
				
				#Delete that element from the original array and add it to your array of selected occurences
				del elements[select]
				chosen	=	chosen + [selected]
				
			#Make an array equal in length  to the number of species in your whole dataset
			replacementRow	=	np.zeros(len(subArray))
			
			#Parse out your selected occs into the right format to put back into data array
			for p in range(len(chosen)):
				speciesNumber	=	int(chosen[p].split('.')[1])
				replacementRow[speciesNumber] += 1
			
			#Replace the original line with your now subsetted line of data
			data[x,5:-1]	=	replacementRow	

	return data, occupiedPixels

def setUpHoldingArray(data):
	#Set up for randomly shuffling while mantaining ecorgion bin width
	unique		=	list(np.unique(data[:,4]))
	numUniques	=	len(unique)
	holding		=	np.empty(len(unique))
	npix		=	len(data)
	maxSpp		=	np.count_nonzero(np.sum(data, axis = 0))

	#Set up yuor first holding array, x-y valude indicates how many pixels are in taht ecoregion (0) and the number of species after that ecoregion is sampled (1)
	x 	= 	0
	y	=	0
	for x in range(npix):
		if x > 0:
			if data[x,4] == oldBiome:
				oldBiome	=	data[x,4]
				x += 1
			else:
				holding[y]	=	x 
				oldBiome	=	data[x,4]
				x += 1
				y += 1
				
				
		else:
			oldBiome	=	data[x,4]
			x += 1

	
	holding[numUniques-1]	=	npix - 1
	
	
	oldSelect	=	0
	for y in range(len(holding)):
		holding[y]	=	holding[y] - oldSelect
		oldSelect 		+= holding[y]

	return holding, unique, numUniques, npix, maxSpp
	
	
def calculatePredictedGrid(data, holding):
	
	#Add function to calculate chao dissimilarity
	def chao(j, k):
		
		#Find species shared between the two sites
		shared	=	np.argwhere((j>=1) & (k>=1)).flatten()
		
		#If there are no shared species return dissimilarity of 1
		if len(shared > 0):
			
			
			#Calculate the total number of individuals in site J and K of all species
			Nj	=	np.sum(j)
			Nk	=	np.sum(k)
			
			#Calculate the number of individuals in site J and K of species in both sites
			Cj	=	np.sum(j[shared])
			Ck	=	np.sum(k[shared])
			

			#Find speices that only occur with 1 individual in site in either J or K that are also present in the other site		
			jones	=	np.argwhere((j>=1) & (k==1)).flatten()
			kones	=	np.argwhere((j==1) & (k>=1)).flatten()
			
			#Find speices that only occur with 2 individuals in site in either J or K that are also present in the other site		
			jtwos	=	np.argwhere((j>=1) & (k==2)).flatten()
			ktwos	=	np.argwhere((j==2) & (k>=1)).flatten()
			
			#Sum up the number of singletons
			a1j	=	len(jones)
			a1k	=	len(kones)
			
			#Sum up doubletons (set to one if there are none as specified in the original paper)
			a2j	=	max(len(jtwos), 1)
			a2k	=	max(len(ktwos), 1)
			
			#Calculate the total number of individuals in site J of species that occur with 1 individual in site K and vice versa
			Sj	=	np.sum(j[jones])
			Sk	=	np.sum(k[kones])
			
			#Calculate indices
			U	=	min(Cj/Nj + (Nk -1)/Nk * a1j/(2*a2j) * Sj/Nj, 1)
			V	=	min(Ck/Nk + (Nj -1)/Nj * a1k/(2*a2k) * Sk/Nk, 1)	

			z = U*V / (U + V - U*V)
				return z
		else:
			return 0

	#Calculate stattistical differences between pixels
	if metric == 'jaccard':
		#Subset out community matrix
		X_true	=	data[:,5:-1]
		
		#Turn into presence absence
		X_true[X_true> 0] = 1
		
		#Calculate jaccard
		statDistances	=	scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(X_true, metric='jaccard', p=2, w=None, V=None, VI=None)).astype('float')
		
		#Convert from dissimilarity to simmilarity
		statDistances	=	np.subtract(1, statDistances)
	
	#Calculate chao
	elif metric == 'chao':
		X_true	=	data[:,5:-1]
		statDistances	=	scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(X_true, chao)).astype('float')

	#Calculate geographic distances
	realDistances	=	scipy.spatial.distance.cdist(data[:,2:3].astype('float'), data[:,2:3].astype('float'), metric='euclidean', p=2, V=None, VI=None, w=None).astype('float')


				
	#Set pixels with no overlapping species or the pixels compared to themselves to nan
	statDistances[realDistances == 0] = 'nan'
	statDistances[statDistances == 0] = 'nan'
	realDistances[realDistances == 0] = 'nan'


	#log community similarity distances and normalize, removing no data values
	statDistances_log 	=	np.log(statDistances)
	statDistances_log 	= 	np.divide(np.subtract(statDistances_log, np.nanmean(statDistances_log)), np.nanstd(statDistances_log))	
	flatStatA 		= 	np.ndarray.flatten(statDistances_log)
	flatStat_log 	= 	flatStatA[~np.isnan(flatStatA)]	
		
	#square root community similarity distances and normalize, removing no data values
	statDistances_sqrt	=	np.sqrt(statDistances)
	statDistances_sqrt 	= 	np.divide(np.subtract(statDistances_sqrt, np.nanmean(statDistances_sqrt)), np.nanstd(statDistances_sqrt))	
	flatStatA 		= 	np.ndarray.flatten(statDistances_sqrt)
	flatStat_sqrt 	= 	flatStatA[~np.isnan(flatStatA)]	
		
	#Normalize linear community similarity distances, removing no data values
	statDistances 	= 	np.divide(np.subtract(statDistances, np.nanmean(statDistances)), np.nanstd(statDistances))	
	flatStatA 		= 	np.ndarray.flatten(statDistances)
	flatStat		= 	flatStatA[~np.isnan(flatStatA)]		
	
	#Noramlize geographic distance values, and remove no data values
	realDistances_linear	=	np.multiply(np.divide(np.subtract(realDistances_linear, np.nanmean(realDistances_linear)), np.nanstd(realDistances_linear)), -1)
	distFlat_linear		= 	np.ndarray.flatten(realDistances_linear)
	distFlat_linear		= 	distFlat_linear[~np.isnan(flatStatA)]	



	#Perform regression between log transformed community composition
	
	d1 = {'distance': distFlat_linear, 'jaccard': flatStat_log}
	df1 = pd.DataFrame(data=d1)
	

	nick	=	smf.ols(formula = 'jaccard~distance', data = df1)
	res		=	nick.fit()

	#Extract paramters
	intercept		=	res.params[0]
	distValue		=	res.params[1]


	#Make grid that predicts communtiy similarity on geographic distance alone
	predictedGrid	=	np.add(np.multiply(realDistances_log, distValue), intercept)
	
	#Subtract observed community similarity from predicted similarity
	observedMinusPredicted_log	=	np.subtract(statDistances,  predictedGrid)
	
	

	#Perform regression between square-root transformed community composition

	d2 = {'distance': distFlat_linear, 'jaccard': flatStat_sqrt}
	df2 = pd.DataFrame(data=d2)
	

	nick	=	smf.ols(formula = 'jaccard~distance', data = df2)
	res		=	nick.fit()
	
	#Extract paramters
	intercept		=	res.params[0]
	distValue		=	res.params[1]

	#Make grid that predicts communtiy similarity on geographic distance alone
	predictedGrid	=	np.add(np.multiply(realDistances_sqrt, distValue), intercept)
	
	#Subtract observed community similarity from predicted similarity
	observedMinusPredicted_sqrt	=	np.subtract(statDistances,  predictedGrid)
	


	#Perform regression between community composition

	d3 = {'distance': distFlat_linear, 'jaccard': flatStat}
	df3 = pd.DataFrame(data=d3)
	

	nick	=	smf.ols(formula = 'jaccard~distance', data = df3)
	res		=	nick.fit()

	#Extract paramters
	intercept		=	res.params[0]
	distValue		=	res.params[1]

	#Make grid that predicts communtiy similarity on geographic distance alone
	predictedGrid	=	np.add(np.multiply(realDistances_linear, distValue), intercept)
	
	#Subtract observed community similarity from predicted similarity
	observedMinusPredicted_linear	=	np.subtract(statDistances,  predictedGrid)
	

	
	return X_true, observedMinusPredicted_log, observedMinusPredicted_sqrt, observedMinusPredicted_linear, statDistances
"""Function to analyze results"""


def analyzeResults(z, output, toSkip):
	#Make final results array
	summaryStats	=	np.zeros(shape=(z, 4))		
	
	#Loop through your output files
	for l in range(z):
		#Bump this up if you wanna skip some files at the beginning
		tt = l + toSkip
			
		#Define your input file for the statistics portion
		pixelBySpeciesPivot	= str(tt) + 'a.csv'

		try:
			
			
			"""Start section on original data management"""
			#Read in species pivot table data
			dataO 	=	np.genfromtxt(pixelBySpeciesPivot, dtype= 'float' , delimiter=',', skip_header=1)

			#Change non-occurences to zeros (important for futrue analysis)
			data	=	np.nan_to_num(dataO)

			#Filter out extra occurances and pixels with too few occurences
			data	=	removeMinOccs(data, 0)
			data, occupiedPixels	=	reduceMaxOccs(data, 50)
			
			#Figure out where ecoregion boundaries fall
			holding, unique, numUniques, npix, maxSpp	=	setUpHoldingArray(data)
			
			"""Calculate your expected grid based on distance alone"""
			X_true, observedMinusPredicted_log, observedMinusPredicted_sqrt, observedMinusPredicted_linear, statDistances	=	calculatePredictedGrid(data, holding)

			
		
			"""Start section on premutations"""
			trials 	=	5000
			gridValues	=	np.zeros(shape=(trials, 4))
			
			#Loop through trials
			for i in range(trials):

				#Do modularity calculation
				a = 0
				b = 0
				j = 0
				
				
				#Loop through your ecoregions
				for j in range(len(unique)):
					b	+=	int(holding[j])
					
					#Add values in same ecoregion comparisons for exponential regression
					z 	=	0
					try:
						z	=	np.nansum(observedMinusPredicted_log[a:b,a:b])
					except:
						z 	= 	0
					
					if np.isfinite(z) == True:
						gridValues[i,0]	+= z
					
					#Add values in same ecoregion comparisons for quadratic regression
					z 	=	0
					try:
						z	=	np.nansum(observedMinusPredicted_sqrt[a:b,a:b])
					except:
						z 	= 	0
					
					
					if np.isfinite(z) == True:
						gridValues[i,1]	+= z
					
					#Add values in same ecoregion comparisons for linear regression
					z 	=	0
					try:
						z	=	np.nansum(observedMinusPredicted_linear[a:b,a:b])
					except:
						z 	= 	0
					
					
					if np.isfinite(z) == True:
						gridValues[i,2]	+= z
					
					
					#Add values in same ecoregion comparisons without correcting for distance
					z 	=	0
					try:
						z	=	np.nansum(statDistances[a:b,a:b])
					except:
						z 	= 	0
					
					if np.isfinite(z) == True:
						gridValues[i,3]	+= z


					 
					a	+=	int(holding[j])
					
				#Shuffle where ecoregion borders fall
				np.random.shuffle(holding)
				 
			#Calculate the percentile for which the actual data value (the first one you calculated) falls relative to the randomly permuted values)
			for k in range(trials):
				if gridValues[k,0] <= gridValues[0,0]:
					summaryStats[l,0] += 1 
				
				if gridValues[k,1] <= gridValues[0,1]:
					summaryStats[l,1] += 1 

				if gridValues[k,2] <= gridValues[0,2]:
					summaryStats[l,2] += 1 
				
				if gridValues[k,3] <= gridValues[0,3]:
					summaryStats[l,3] += 1 
			
			
			summaryStats[l,0]	=	float(summaryStats[l,0]) / float(trials)
			summaryStats[l,1]	=	float(summaryStats[l,1]) / float(trials)
			summaryStats[l,2]	=	float(summaryStats[l,2]) / float(trials)
			summaryStats[l,3]	=	float(summaryStats[l,3]) / float(trials)

			#Write to output files
			with open(output, 'ab') as csvfile:
				spamwriter = csv.writer(csvfile, delimiter=',')
				spamwriter.writerow([tt, maxSpp, np.sum(X_true), len(X_true), occupiedPixels, np.count_nonzero(~np.isnan(statDistances)), numUniques, unique, summaryStats[l,0], summaryStats[l,1],summaryStats[l,2],summaryStats[l,3]])
			print [tt, maxSpp, np.sum(X_true), len(X_true), occupiedPixels, np.count_nonzero(~np.isnan(statDistances)), numUniques, unique, summaryStats[l,0], summaryStats[l,1],summaryStats[l,2],summaryStats[l,3]]
			
		#If it doesn't work write a line of noData values
		except:
			print tt
			summaryStats[l,:] = np.nan
			with open(output, 'ab') as csvfile:
				spamwriter = csv.writer(csvfile, delimiter=',')
				spamwriter.writerow([tt, -9999,  -9999,  -9999,  -9999,  -9999, -9999, -9999,  -9999, -9999, -9999, -9999])
				
				
		
	
def runCode(dir, output):
	try:
		os.makedirs(dir + '/summary/')
	except:
		a = 1

	
	#Determine number of ouptut files
	os.chdir(dir)
	
	analyzeResults(5000, output, 0)		
	

metric = 'jaccard'

#Change directory and make temp & output 
dir		=	'/scratch/PI/gdaily/Jeff/GIS_Layers/GBIF/Amphibians/Outputs'
output 	=	'summary/ddm_summary_2001_jaccard_0.csv'
runCode(dir, output)




