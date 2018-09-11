#Import necessary libraries and functions
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from math import radians, sin, cos, tan, pi, sqrt
import pyproj

"""Start User Inputs"""

#What type of data you would like to use
obserationTypes	=	['HUMAN_OBSERVATION',
					'PRESERVED_SPECIMEN', 
					'SPECIMEN', 
					'OBSERVATION', 
					'LIVING_SPECIMEN', 
					'MATERIAL_SAMPLE', 
					'MACHINE_OBSERVATION']

#Set year cutoff
minYear	=	1950

#Set precision estimate cutoff
precision = 10000

#Set cutoff for number of ecoregion represented for a dataset to be included
ecoregionMin	=	1

#File path to ascii location
ascii	=	"/scratch/PI/gdaily/Jeff/GIS_Layers/Ecoregion/ecor_2017_3410.asc"


"""End user inputs"""
#Process the headers of ascii 
asciiParameters				=	pd.read_csv(ascii, nrows = 5, header = None, delimiter = " ")
asciiParameters['Value'] 	= 	asciiParameters.sum(axis = 1, skipna = True)


#Save ascii properties as varaibles for use later
ncols 			= 	int(asciiParameters['Value'][0])
nrows 			= 	int(asciiParameters['Value'][1])
xorig	 		= 	asciiParameters['Value'][2]
yorig	 		= 	asciiParameters['Value'][3]
cellsize		= 	asciiParameters['Value'][4]


#Read in ascii minus the headers to a numpy array
ecor	=	np.genfromtxt(ascii,  dtype = int, skip_header = 5)

#Define which projection systems you might need
wgs84	=	pyproj.Proj("+init=EPSG:4326") 
p3410	=	pyproj.Proj("+init=EPSG:3410")

										
	
def cleanGBIF(directory, toSave):
	
	#Read in raw GBIF data
	inputFile	=	str(directory + 'raw.csv')
	full = pd.read_csv(inputFile, header = 0, sep = '\t',
						usecols 	= 	['species',
										'decimallatitude',
										'decimallongitude',
										'coordinateuncertaintyinmeters',
										'year',
										'basisofrecord',
										'datasetkey'],
										
						dtype		=	 {'species': str, 
										'decimallatitude': str,
										'decimallongitude': str,
										'coordinateuncertaintyinmeters': str,
										'year':  str,
										'basisofrecord': str,
										'datasetkey': str
										},
						error_bad_lines = False)
	#Change data types here (this can't be done in the step above because of bad data-entry (trust us we tried it works for some GBIF data, but not all)
	full['decimallatitude']						=	pd.to_numeric(full['decimallatitude'], errors='coerce')
	full['decimallongitude']					=	pd.to_numeric(full['decimallongitude'], errors='coerce')
	full['year']								=	pd.to_numeric(full['year'], errors='coerce')
	full['coordinateuncertaintyinmeters']		=	pd.to_numeric(full['coordinateuncertaintyinmeters'], errors='coerce')


	
	#Keep only those with valid records for species, lat, and long
	full 	= 	full.dropna(subset=['species', 'decimallatitude','decimallongitude']) 
	
	#Drop records that are too inpercise
	#full 	= 	full.loc[full['coordinateuncertaintyinmeters'] <= precision]
	
	#Drop points with 0.0 lat and/or long
	full 	= 	full.loc[full['decimallatitude']  != 0.0]
	full 	= 	full.loc[full['decimallongitude'] != 0.0]
	
	#Drop records that are too old
	full 	= 	full.loc[full['year'] >= minYear]
	
	#Make sure it is in your selected observation types
	full 	= 	full[full['basisofrecord'].isin(obserationTypes)]


	#Convert lat, long to x,y
	toExtract	=	full.as_matrix(columns = ['decimallongitude', 'decimallatitude'])
	full['x'], full['y']		=	p3410(toExtract[:,0], toExtract[:,1])
		
	#Calculate which pixel of the ecoregion ASCII that puts you into
	full['xpix']	=	((full['x'] - xorig) / cellsize).astype(int)
	full['ypix']	=	(nrows  - ((full['y']	 - yorig).astype(int) / cellsize)).astype(int)
	
	#Find out what ecoregion each pixel is in and drop those that have a value of 0
	full['ecor']	=	ecor[full['ypix'], full['xpix']]
	full			=	full.loc[full['ecor'] != 0]

	
	#Pivot table to calculate how many unique ecoregions each dataset represents
	full['ones']	=	1
	screening		=	pd.pivot_table(full, columns = ['ecor'], index=['datasetkey'], values = ['ones'], aggfunc=np.sum)
	screening		=	screening.fillna(0)
	screening[screening > 1]		=	1
	
	#Sum up the number of unique ecoregions in each dataset
	toMatch			=	screening.sum(axis = 1)
	
	#Drop those which have too few ecoregions represented and then do some data management
	toMatch			=	toMatch[toMatch >= ecoregionMin].rename('ecoregionsRepresented')
	toMatch			=	toMatch.reset_index()

	#Drop all records that aren't in the datasets you've decided to keep
	full			=	pd.merge(full, toMatch, on = 'datasetkey', how = 'inner')

	#Define what you want to save
	toKeep				=	full[['xpix', 'ypix']].astype(int)
	toKeep['Species']	=	full['species']
	
	#Save some data
	pathToSave		=	str(directory + toSave)
	toKeep.to_csv(pathToSave, sep = ',', index = False)

cleanGBIF('/scratch/PI/gdaily/Jeff/GIS_Layers/GBIF/Birds/', 'cleaned_no_precision.csv')

