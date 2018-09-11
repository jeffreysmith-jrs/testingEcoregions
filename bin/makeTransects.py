"""

This goes through a cleaned GBIF dataset (see cleaningGBIF.py) and makes transects of species communities
The inputs are at 1 map ecoregions and any cleaned CSV from GBIF

"""

#Written by Jeffrey R. Smith
#Last updated May 24, 2018


#Import libraries
import os, sys, random, csv, time, shutil, fnmatch, scipy
import pandas as pd
import numpy as np
from pandas import DataFrame




"""Define your dependent funcitons"""

#Define where you are going to sample
def makeAOI(t):
	
	#Set the constraints for how long a transect has to be to be included (note some of these pixels may not have occs)
	minPix	=	750
	maxPix	=	3000
	
	#Make AOIs until you have one of sufficient length
	sufficientLength	=	0
	while sufficientLength < 1:
		
		#Make array to track pixels in your AOI
		samplingArray	=	np.zeros(shape=(1,2))
		
		#Make list to track biomes you've used
		biomesUsed		=	[]
		
		#Use the b variable to make sure that you're meeting various conditions as you generate transects, othrwise you can start over
		b	=	0
		while b < 1:
			
			
			#Make sure you strat on a land pixel
			findStart = 0
			while findStart	< 1:	
			
				#Pick a random pixel
				xPixel	=	random.randint(0,ncols-2)
				yPixel	=	random.randint(0,nrows-2)
			
				
				#If you're on a land pixel you can break out of the first while loop
				if biome[yPixel,xPixel] > 0 and biome[yPixel,xPixel] < 20:
					b = 1
					samplingArray[0,0]	= xPixel
					samplingArray[0,1]	= yPixel
					biomesUsed			= [int(ecoreg[yPixel,xPixel])]
					
					findStart	=	1
					
			#Loop through pixels until you reach the number of max pixels	
			k		=	0
			exit 	=	0
			
			
			while k < maxPix:
				d	=	0
				
				#Define the possible moves for your AOI
				moves	=	[(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1), (-2, 0), (2, 0), (0, -2), (0, 2), (-2, -1), (2, -1), (-1, -2), (-1, 2), (-2, 1), (2, 1), (1, -2), (1, 2)]
				
				#Loop through until you find a valid pixel
				while d < 1:
				

					#Choose move at random
					dx, dy = random.choice(moves)
					xPixel	= 	int(samplingArray[k-1,0] + dx)
					yPixel	=	int(samplingArray[k-1,1] + dy)			
						
					#Track which moves you've used (in case you need to remove it from array)
					moved	=	dx * 10 + dy
						
					#Check to see if the pixel you just moved into was already in the array
					e = 0
					c = 0 
					for c in range(len(samplingArray)):
						if xPixel == samplingArray[c,0] and yPixel == samplingArray[c,1]:
							e += 1
							moves = [i for i in moves if (i[0] * 10 + i[1]) != moved]
						c += 1
											
					#If pixel is off the side of the earth have it loop around to the other side of the world, fix the y coordinate part
					if xPixel == ncols - 1:
						xPixel	= 	0
						
					elif xPixel == ncols:
						xPixel	= 	1
					
					if yPixel == nrows - 1:
						yPixel	=	0	
						
					elif yPixel == nrows:
						yPixel	=	1		
						
					#If that pixel wasn't in the array of pixels you've already used you can keep going
					if e == 0: 
						
						#Check to see if you're in the same biome that you were in last time, if so continue sampling as normal
						xPixel = int(xPixel)
						yPixel = int(yPixel)


						if ecoreg[int(yPixel),int(xPixel)] == ecoreg[int(samplingArray[k-1,1]),int(samplingArray[k-1,0])]:
							
							#Exit the d loop which means you found a viable pixel
							d = 1
							break
						
						#If you've entered a different biome check to see if you've been in that biome before
						else:
							
							#If it is not in that biome, great start sampling that biome 
							if int(ecoreg[yPixel,xPixel]) not in biomesUsed and int(ecoreg[yPixel,xPixel]) != -9999 and int(ecoreg[yPixel,xPixel]) != 0:
								
								#Exit the d loop which means you found a viable pixel
								biomesUsed	=	biomesUsed + [int(ecoreg[yPixel,xPixel])]
								d = 1
								break
								
							#If you've already been to that biome you can't go back so break out of the loop and remove that move from the options	 							
							else:
								moves = [i for i in moves if (i[0] * 10 + i[1]) != moved]
								d = 0
					
					
					#If you've used all the moves and still have not found a viable pixel then you're done trying to make your transect, we'll now check if its long enough
					if len(moves) == 0:
						exit = 1
						d = 1
						k = maxPix
						break
					
				#If you've been told to exit because there are no more good pixels go onto checking to see if you're long enough	
				if exit == 1:
					k = maxPix
					break	
				
				#Make sure its not water, it it is then your transect is over and you have to check if you're long enough
				elif biome[yPixel,xPixel] <= 0 or biome[yPixel,xPixel] > 20:
					k = maxPix
					break	
				
				#If you have viable pixels that aren't water, great, you can add that pixel to the sampling array
				else:
					samplingArray = 	np.append(samplingArray, [[xPixel, yPixel]], axis = 0)
					
		#Check to see if your biome is sufficiently long if so you're done	
		if len(samplingArray) >= minPix:
			sufficientLength = 1
		
		#If not start the whole thing over again
		else:
			k	=	0
			

	#Make up your array that will be you AOI
	pixels		= 	np.zeros(shape=(len(samplingArray),4))
	npix		=	len(pixels)


	#Fill your transect array
	for c in range(npix):
		
		#From the random walk AOI we just amade
		xpixel		=	samplingArray[c,0]
		ypixel		=	samplingArray[c,1]
		
		#Lower x bound, upper x bound, lower y bound, upper y bound, x pixel, y pixel, biome value, ecoregion value, pixel number
		pixels[c,0]	=	c
		pixels[c,1]	=	xpixel
		pixels[c,2]	=	ypixel
		pixels[c,3]	= 	ecoreg[int(ypixel), int(xpixel)]
		c += 1
	
	#Save your AOI
	return(pixels)


#Find the occs in each pixel in your AOI		
def extractSpeciesPerPixel(pixelBySpeciesPivot, pixels, taxa):		
	
	try:
		os.remove(pixelBySpeciesPivot)
	except:
		aaa = 1

	#Set up transect AOI dataframe
	transect	=	pd.DataFrame(data = pixels, columns=['pixNumb', 'xpixA', 'ypixA', 'ecorA']).astype(int)
	
	
	#Sample that tarnsect for species of focal taxa
	df	=	pd.merge(left = transect, right = taxa, how='left', left_on=['xpixA', 'ypixA'], right_on=['xpix', 'ypix'])
		
	#Sum up the number of each species in each pixel and do a bit of dataframe management to reset column and index levels
	reshape		=	pd.pivot_table(df,index=['pixNumb'],values=["Count"],columns=["Species"], aggfunc = np.sum, margins = True, margins_name = 'total')
	
	reshape.columns = reshape.columns.droplevel()
	reshape		=	reshape.reset_index()


	#Rejoin to transet to make sure that you aren't dropping pixels that didn't have any point occurences (which you lost in the last step)
	reshape2	=	pd.merge(left = transect, right = reshape, how = 'left', left_on = ['pixNumb'], right_on=['pixNumb'])
	
	
	#Save to CSV
	reshape2.to_csv(pixelBySpeciesPivot, sep = ',', mode = 'w', header = True)	


"""End Functions section"""
	





"""Define global variables"""
#Start timer
startTime = time.time()

#Import basemaps (ecoregion asciis) - Always use the 2001 biome baseamp (the values work better with the code)
ascii 					= 	"/scratch/PI/gdaily/Jeff/GIS_Layers/Ecoregion/biom_2001_3410.asc"
asciiE 					= 	"/scratch/PI/gdaily/Jeff/GIS_Layers/Ecoregion/ecor_2001_3410.asc"
biome 					= 	np.genfromtxt(ascii,  dtype = int, skip_header = 6)
ecorega					=	np.genfromtxt(asciiE,  dtype = int, skip_header = 6)

#Get rid of bad data values
ecorega[ecorega < 9999] = 0
biome[biome > 9999] = 0


#Compare to the 2001 map and remove pixels that are not present on both maps
ascii2017 					= 	"/scratch/PI/gdaily/Jeff/GIS_Layers/Ecoregion/ecor_2017_3410.asc"
ecoreg2017					=	np.genfromtxt(ascii2017,  dtype = int, skip_header = 6)
ecorega[ecoreg2001 > 9999] = 0

#Save your final output map of ecoregions 
ecoreg = np.copy(ecorega)


#Process the headers of ascii 
asciiParameters				=	pd.read_csv(ascii, nrows = 5, header = None, delimiter = " ")
asciiParameters['Value'] 	= 	asciiParameters.sum(axis = 1, skipna = True)


#Save ascii properties as varaibles for use later
ncols 			= 	int(asciiParameters['Value'][0])
nrows 			= 	int(asciiParameters['Value'][1])
xorig	 		= 	asciiParameters['Value'][2]
yorig	 		= 	asciiParameters['Value'][3]
cellsize		= 	asciiParameters['Value'][4]

print nrows, ncols

#Set root directory
rootDirectory	=	'/scratch/PI/gdaily/Jeff/GIS_Layers/GBIF/'

#Dictate directories where you have your cleaned CSVs
taxaDirectories	=	['Plants', 'Mammals', 'Arthropods', 'Reptiles', 'Amphibians', 'Birds', 'Fungi']


#Make output directories and read-in each taxa's data
for k in range(len(taxaDirectories)):
	outputDirectory		=	str(str(rootDirectory) + taxaDirectories[k] + '/Outputs/')
	try:
		os.makedirs(outputDirectory)
	except:
		print 'already exists'

#Define locations of your directories
plant_dir	=	str(str(rootDirectory) + '/Plants/')
bug_dir		=	str(str(rootDirectory) + '/Arthropods/')
mam_dir		=	str(str(rootDirectory) + '/Mammals/')
bird_dir	=	str(str(rootDirectory) + '/Birds/')
fungi_dir	=	str(str(rootDirectory) + '/Fungi/')
herp_dir	=	str(str(rootDirectory) + '/Herps/')
reptile_dir	=	str(str(rootDirectory) + '/Reptiles/')
amphib_dir	=	str(str(rootDirectory) + '/Amphibians/')

	 
#Read in each of your GBIF occurence datasets 		
plants		=	pd.read_csv('/scratch/PI/gdaily/Jeff/GIS_Layers/GBIF/Plants/cleaned0.csv', 		sep=',', header = 0, dtype={'xpix': np.int, 'ypix':np.int, 'Species':str})
bugs		=	pd.read_csv('/scratch/PI/gdaily/Jeff/GIS_Layers/GBIF/Arthropods/cleaned0.csv', 		sep=',', header = 0, dtype={'xpix': np.int, 'ypix':np.int, 'Species':str})
mammals		=	pd.read_csv('/scratch/PI/gdaily/Jeff/GIS_Layers/GBIF/Mammals/cleaned0.csv', 		sep=',', header = 0, dtype={'xpix': np.int, 'ypix':np.int, 'Species':str})
fungi		=	pd.read_csv('/scratch/PI/gdaily/Jeff/GIS_Layers/GBIF/Fungi/cleaned0.csv', 		sep=',', header = 0, dtype={'xpix': np.int, 'ypix':np.int, 'Species':str})
birds		=	pd.read_csv('/scratch/PI/gdaily/Jeff/GIS_Layers/GBIF/Birds/cleaned0.csv', 		sep=',', header = 0, dtype={'xpix': np.int, 'ypix':np.int, 'Species':str})
reptiles 		=	pd.read_csv('/scratch/PI/gdaily/Jeff/GIS_Layers/GBIF/Reptiles/cleaned0.csv', 	sep=',', header = 0, dtype={'xpix': np.int, 'ypix':np.int, 'Species':str})
amphibians 	=	pd.read_csv('/scratch/PI/gdaily/Jeff/GIS_Layers/GBIF/Amphibians/cleaned0.csv', 		sep=',', header = 0, dtype={'xpix': np.int, 'ypix':np.int, 'Species':str})



#Make a column of ones to sum up for future analysis 
plants['Count'] = 1
bugs['Count'] = 1
mammals['Count'] = 1
fungi['Count'] = 1
birds['Count'] = 1
reptiles['Count'] = 1
amphibians['Count'] = 1





#Transects to skip (for running multiple codes at once on Sherlock)
transectsToSkip = 0

"""End Global variables"""

#Do a certain number of runs (change what is in the range)
for l in range(10):
	
	#Add this so you don't overwrite exisiting output files you may have (change the number to what you want to kep)
	t	=	l + transectsToSkip
	
	#Run AOI code that you will use for all of the taxa you are considering
	pixels	=	makeAOI(t)
	print 'Made AOI after', (time.time() - startTime)/60, 'minutes of total program time'

	#Select the species on that transect for each taxa
	def taxaSpecific(dir, pixels, taxa):
		try:
			#Navigate to your correct directory
			os.chdir(rootDirectory)
			os.chdir(dir)
			
			#Run your transect
			pixelBySpeciesPivot 	= 'Outputs/' + str(t) + 'b.csv'
			extractSpeciesPerPixel(pixelBySpeciesPivot, pixels, taxa)	
			
		#Put in an out in case there are no occcurences on your transect (i.e. Sahara, Siberia) prevents pandas pivot table from failing
		except:
			print 'Did not find any point occurences in this transect'
	#Run for each taxa
	taxaSpecific(plant_dir, 	pixels,	 	plants)	
	taxaSpecific(mam_dir,  		pixels,		mammals)
	taxaSpecific(bug_dir,  		pixels, 	bugs)
	taxaSpecific(bird_dir, 		pixels, 	birds)
	taxaSpecific(fungi_dir, 	pixels, 	fungi)
	taxaSpecific(reptile_dir, 	pixels, 	reptiles)
	taxaSpecific(amphib_dir, 	pixels, 	amphibians)
	
	
	
	#Bump up counter
	t += 1
	print 'Run', t, 'completed after', (time.time() - startTime)/60, 'minutes of total program time'
