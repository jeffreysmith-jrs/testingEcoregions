"""

This goes through a CSVs for a given taxa (or multiple taxa if looped through), and calculates statistics for each transect using our 2 main tests
The species-discoerty curves (fitting step-wise functions and calculating residuals against residuals from random permutations of step-wise function)


"""


#Written by Jeffrey R. Smith
#Last updated May 30, 2017



#Import libraries
import os,sys,random,csv,sys,time,shutil,fnmatch,scipy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sin, cos, sqrt, log, radians, fabs
from matplotlib import pyplot as plt
import seaborn as sns
import statsmodels.formula.api as smf

start = time.time()

"""Function to analyze results"""
def analyzeResults(w):
	
	#Get rid of existing output files
	try:
		os.remove('summary/summary2017.csv')
	except:
		a  = 1 
	
	#Loop through your output files
	for l in range(w):
		print l
		
		#Bump this up if you wanna skip some files at the beginning
		tt = l + 0
				
		#Define your input file for the statistics portion
		pixelBySpeciesPivot	= str(tt) + 'c.csv'
		
		if os.path.exists(pixelBySpeciesPivot):
		
			"""Start section on original data management"""
			#Read in species pivot table data
			dataO 	=	np.genfromtxt(pixelBySpeciesPivot, dtype= 'float' , delimiter=',', skip_header=1)

			#Change non-occurences to zeros (important for futrue analysis)
			data	=	np.nan_to_num(dataO)
		
						
			#Set max occs and count up occupied pixels
			maxOccs	=	50
			occupiedPixels	=	np.count_nonzero(data[:,-1])
			

			
			
			#Randomize overrepresented pixels
			npix	=	len(data)
			stats	=	np.zeros(shape=(npix, 8))
			
			#Loop through pixels
			x = 0
			for x in range(npix):
			
				#Identify pixels with too many occs
				if data[x,-1] > maxOccs:
					
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
			"""End section on original data management"""
			
			"""Start section on AUC"""

			#Holding for data frame
			dataR		=	data[:,5:-1].astype(int)
		
			#Set up for randomly shuffling while mantaining ecorgion bin width
			unique		=	list(np.unique(data[:,4]))
			numUniques	=	len(unique)
			holding		=	np.zeros(shape=(numUniques,2), dtype = np.int32)
			npix		=	len(dataR)
			maxSpp		=	np.count_nonzero(np.sum(dataR, axis = 0))


			#Calculate the number of total number of species seen as a function of the number of pixels you have crossed
			x = 0 
			numSpecies	=	np.zeros(npix)


			for x in range(npix-1):
				numSpecies[x]	=	np.count_nonzero(np.sum(dataR[0:x+1,:], axis = 0))
			numSpecies[npix-1] = maxSpp			
				
				

			#Set up your first holding array that indicates the pixel value for when you cross into a new ecoregion
			x 	= 	0
			y	=	0
			for x in range(npix):
				if x > 0:
					if data[x,4] == oldBiome:
						oldBiome	=	data[x,4]
						x += 1
					else:
						holding[y,0]	=	x 
						oldBiome	=	data[x,4]
						x += 1
						y += 1
						
						
				else:
					oldBiome	=	data[x,4]
					x += 1
		
			
			holding[-1,0]	=	npix - 1
			
			
						
			#Change it so the first column of the holding array conains the number of pixels in each ecoregion (from the cumulative number of pixels)
			oldSelect	=	0
			for y in range(len(holding)):
				holding[y,0]	=	holding[y,0] - oldSelect
				oldSelect 		+= holding[y,0]
		
			
			#Calculate the total number of species you have encountered once you finsih sampling an ecoregion 
			x	=	0
			y	=	0
				
				
			for y in range(numUniques):
					
				x += holding[y,0]
				holding[y,1]	=	numSpecies[x]
			holding[-1,1]	=	maxSpp	
			
			print holding
				

				
			
					
		
			#Make permutation arrays
			trials		=	5000
			aucValue		=	np.zeros(trials)

			#Loop through for each permutation
			#The first time we loop through we do it for the real ecoregion boundaries and then loop through 4999 random permutations
			for w in range(trials):
				
				#Calculate the stepwise function based on ecoregion width and how many species you've found cumulatively once you're done sampling each ecoregion 
				modelSpecies	=	np.zeros(npix)	
				u	=	0
				v	=	0
				x	=	0			
				for v in range(numUniques):
					for u in range(int(holding[v,0])):
						modelSpecies[x] = holding[v,1]
						x += 1
					u = 0
					v += 1
				modelSpecies[-1] = maxSpp

				#Calculate the residual between the actual number of species seen and the modelled stepwise function based on known or random ecoregion boundairs
				aucValue[w]	=	np.sum((numSpecies - modelSpecies)**2)				
		


								
				
				#Randomly shuffle by ecoregion								
				np.random.shuffle(holding)
				
				x	=	0
				y	=	0
				
				#Recalculate the number of species you would've found in each of the 
				for y in range(numUniques):
					
					x += holding[y,0]
					holding[y,1]	=	numSpecies[x]
				holding[-1,1]	=	maxSpp	
		
			"""End section on AUC"""
			
			
			"""Figure out p values"""
			count_SAC	=	0	
			
			#Check to see if the permutation residual is less than or greater than the residual from the known ecoregion bounds (first loop through) 
			x	=	0
			for x in range(len(aucValue)-1):
				x += 1
				if aucValue[x] < aucValue[0]:
					count_SAC	+= 1
				

			"""End p values section"""
			
			#Write to output file
			output	=	'summary/summary2017.csv'

					
			with open(output, 'ab') as csvfile:
				spamwriter = csv.writer(csvfile, delimiter=',')
				spamwriter.writerow([float(count_SAC)/float(len(aucValue)), '0', '0', '0', numUniques, np.nansum(dataR), npix, maxSpp, occupiedPixels, unique, tt])
			print float(count_SAC)/float(len(aucValue)), '0', '0', '0', numUniques, np.nansum(dataR), npix, maxSpp, occupiedPixels, unique, tt
			

		
		else:
			output	=	'summary/summary2017.csv'
					
			with open(output, 'ab') as csvfile:
				spamwriter = csv.writer(csvfile, delimiter=',')
				spamwriter.writerow([-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, tt])
		print time.time()-start
		
"""End analysis"""




#Change directory and make temp & output directories
dir		=	'/scratch/PI/gdaily/Jeff/GIS_Layers/GBIF/Plants/Outputs'
try:
	os.makedirs(dir + '/summary/')
except:
	a = 1

#Determine number of ouptut files

os.chdir(dir)
analyzeResults(10000)



