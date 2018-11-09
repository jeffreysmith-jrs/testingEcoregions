

#Import libraries
import os,sys,random,csv,sys,time,shutil,fnmatch,scipy,sklearn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sin, cos, sqrt, log, radians, fabs
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA
import seaborn as sns
from scipy.stats import mstats
from scipy.optimize import curve_fit
from matplotlib import cm
from matplotlib.gridspec import GridSpec
import ecopy as ep
from collections import Counter
import statsmodels.formula.api as smf
import matplotlib.image as image

def makeGraph(year, minOccPixels, minPoints, minSpecies, minEcoregions, allGraphs, outfile):
	#Import summary files
	plants	=	np.genfromtxt(str('E:\\GIS_layers\\GBIF\\Plants\\' + year + '.csv'),  	dtype = float, delimiter = ',', usecols = (0,1,2,3,4,5,6,7,8))
	mammals	=	np.genfromtxt(str('E:\\GIS_layers\\GBIF\\Mammals\\' + year + '.csv'),  	dtype = float, delimiter = ',', usecols = (0,1,2,3,4,5,6,7,8))
	herps	=	np.genfromtxt(str('E:\\GIS_layers\\GBIF\\Herps\\' + year + '.csv'),  		dtype = float, delimiter = ',', usecols = (0,1,2,3,4,5,6,7,8))
	bugs	=	np.genfromtxt(str('E:\\GIS_layers\\GBIF\\Arthropods\\' + year + '.csv'), 	dtype = float, delimiter = ',', usecols = (0,1,2,3,4,5,6,7,8))
	fungi	=	np.genfromtxt(str('E:\\GIS_layers\\GBIF\\Fungi\\' + year + '.csv'),  		dtype = float, delimiter = ',', usecols = (0,1,2,3,4,5,6,7,8))
	birds	=	np.genfromtxt(str('E:\\GIS_layers\\GBIF\\Birds\\' + year + '.csv'),  		dtype = float, delimiter = ',', usecols = (0,1,2,3,4,5,6,7,8))
	reptiles	=	np.genfromtxt(str('E:\\GIS_layers\\GBIF\\Reptiles\\' + year + '.csv'),  		dtype = float, delimiter = ',', usecols = (0,1,2,3,4,5,6,7,8))
	amphibians	=	np.genfromtxt(str('E:\\GIS_layers\\GBIF\\Amphibians\\' + year + '.csv'),  		dtype = float, delimiter = ',', usecols = (0,1,2,3,4,5,6,7,8))



	#Subset based ont he minimum number of occupied pixels
	plants	=	plants[plants[:,8] 		>=minOccPixels]
	mammals	=	mammals[mammals[:,8] 	>=minOccPixels]
	bugs	=	bugs[bugs[:,8] 			>=minOccPixels]
	fungi	=	fungi[fungi[:,8] 		>=minOccPixels]
	herps	=	herps[herps[:,8] 		>=minOccPixels]
	birds	=	birds[birds[:,8] 		>=minOccPixels]
	reptiles	=	reptiles[reptiles[:,8] 		>=minOccPixels]
	amphibians	=	amphibians[amphibians[:,8] 		>=minOccPixels]


	#Summarize by the minimum number of unique species
	plants	=	plants[plants[:,7]	 	>=minSpecies]
	mammals	=	mammals[mammals[:,7] 	>=minSpecies]
	bugs	=	bugs[bugs[:,7] 			>=minSpecies]
	fungi	=	fungi[fungi[:,7] 		>=minSpecies]
	herps	=	herps[herps[:,7] 		>=minSpecies]
	birds	=	birds[birds[:,7] 		>=minSpecies]
	reptiles	=	reptiles[reptiles[:,8] 		>=minOccPixels]
	amphibians	=	amphibians[amphibians[:,8] 		>=minOccPixels]


	#Summarize by the minimum number of ecoregions on your transect
	plants	=	plants[plants[:,4] 		>=minEcoregions]
	mammals	=	mammals[mammals[:,4] 	>=minEcoregions]
	bugs	=	bugs[bugs[:,4] 			>=minEcoregions]
	fungi	=	fungi[fungi[:,4] 		>=minEcoregions]
	herps	=	herps[herps[:,4] 		>=minEcoregions]
	birds	=	birds[birds[:,4] 		>=minEcoregions]
	reptiles	=	reptiles[reptiles[:,8] 		>=minOccPixels]
	amphibians	=	amphibians[amphibians[:,8] 		>=minOccPixels]


	#Summarize by the minimum number of total points
	plants	=	plants[plants[:,5] 		>=minPoints]
	mammals	=	mammals[mammals[:,5]	>=minPoints]
	bugs	=	bugs[bugs[:,5] 			>=minPoints]
	fungi	=	fungi[fungi[:,5] 		>=minPoints]
	herps	=	herps[herps[:,5] 		>=minPoints]
	birds	=	birds[birds[:,5] 		>=minPoints]
	reptiles	=	reptiles[reptiles[:,8] 		>=minOccPixels]
	amphibians	=	amphibians[amphibians[:,8] 		>=minOccPixels]


	if allGraphs == 'yes':
		print float(len(plants[plants[:,0] <= 0.05]))	/float(len(plants)), 'plants'
		print float(len(bugs[bugs[:,0] <= 0.05]))		/float(len(bugs)), 'bugs'
		print float(len(mammals[mammals[:,0] <= 0.05]))	/float(len(mammals)), 'mammals'
		print float(len(herps[herps[:,0] <= 0.05]))		/float(len(herps)), 'herps'
		print float(len(birds[birds[:,0] <= 0.05]))		/float(len(birds)), 'birds'
		print float(len(fungi[fungi[:,0] <= 0.05]))		/float(len(fungi)), 'fungi'
		print float(len(reptiles[reptiles[:,0] <= 0.05]))		/float(len(reptiles)), 'reptiles'
		print float(len(amphibians[amphibians[:,0] <= 0.05]))		/float(len(amphibians)), 'amphibians'



		tests = ['sac']
		for z in range(len(tests)):
			taxa = [plants, bugs, herps, mammals, birds, fungi, reptiles, amphibians]
			taxaS = ['plants', 'bugs', 'herps', 'mammals', 'birds', 'fungi', 'reptiles', 'amphibians']
			
			
			for x in range(len(taxa)):
				statist, p = scipy.stats.kstest(taxa[x][:,z], 'norm')
				print tests[z], taxaS[x], 'normal', statist, p
					
					
				
				statist, p = scipy.stats.kstest(taxa[x][:,z], 'uniform')
				print tests[z], taxaS[x], 'uniform', statist, p
				
				for y in range(len(taxa)):
					statist, p = scipy.stats.ks_2samp(taxa[x][:,z], taxa[y][:,z])
					
					print tests[z], taxaS[x], taxaS[y], statist, p
					y += 1
				x += 1
				y = 0




	#Define plot style
	plt.style.use(u'seaborn-white')
	fig1 	= 	plt.figure()
	gs	= GridSpec(2,4)
	
	if allGraphs == 'yes':
		fig1.set_size_inches(16,8)
		ax		=	fig1.add_subplot(gs[:,1:3])
		
	else:
		fig1.set_size_inches(10,8)
		ax		=	fig1.add_subplot(gs[:,:])

	#Plot transect specific p-values from Species accumulation curves
	
	data	=	[fungi[:,0], bugs[:,0], birds[:,0],reptiles[:,0],  plants[:,0],  mammals[:,0], amphibians[:,0]]


	plt.violinplot(data, showmeans=False, showextrema=False, showmedians=True, vert = False)
	ax.set_yticks([7,6,5,4,3,2,1])


	ax.text(0.5, 5.25, 'Plants ' + '(n = ' + str(len(plants)) + ')', fontsize = 16, horizontalalignment='center')
	ax.text(0.5, 2.25, 'Arthropods ' + '(n = ' + str(len(bugs)) + ')', fontsize = 16, horizontalalignment='center')
	ax.text(0.5, 4.25, 'Reptiles ' + '(n = ' + str(len(reptiles)) + ')', fontsize = 16, horizontalalignment='center')
	ax.text(0.5, 7.25, 'Amphibians ' + '(n = ' + str(len(amphibians)) + ')', fontsize = 16, horizontalalignment='center')
	ax.text(0.5, 6.25, 'Mammals ' + '(n = ' + str(len(mammals)) + ')', fontsize = 16, horizontalalignment='center')
	ax.text(0.5, 3.25, 'Birds ' + '(n = ' + str(len(birds)) + ')', fontsize = 16, horizontalalignment='center')
	ax.text(0.5, 1.25, 'Fungi ' + '(n = ' + str(len(fungi)) + ')', fontsize = 16, horizontalalignment='center')

	#ax.set_yticklabels(label, rotation = 45, ha = 'right')
	ax.set(xlim=(0, 1))
	ax.set(ylim=(0.5,7.5))
	ax.set_title('Species-Accumulation Curves', weight = 'bold', fontsize = 22)
	ax.set_xlabel(	u"\u2190" + '  More supportive of Sharp-Transition Hypothesis', fontsize = 16)
	ax.set_yticklabels([])

	im = image.imread('C:/Users/Jeffrey/Documents/Academic/Stanford/Steg/Drafts/v9/taxa/plant.tif')
	ax.imshow(im, aspect='auto', extent=(0.9, 0.95, 5, 5.5))

	im = image.imread('C:/Users/Jeffrey/Documents/Academic/Stanford/Steg/Drafts/v9/taxa/bug.png')
	ax.imshow(im, aspect='auto', extent=(0.89, 0.96, 2, 2.5))

	im = image.imread('C:/Users/Jeffrey/Documents/Academic/Stanford/Steg/Drafts/v9/taxa/amphibian.tif')
	ax.imshow(im, aspect='auto', extent=(0.87, 0.97, 7, 7.5))

	im = image.imread('C:/Users/Jeffrey/Documents/Academic/Stanford/Steg/Drafts/v9/taxa/mammal.tif')
	ax.imshow(im, aspect='auto', extent=(0.87, 0.97, 6, 6.5))

	im = image.imread('C:/Users/Jeffrey/Documents/Academic/Stanford/Steg/Drafts/v9/taxa/bird.tif')
	ax.imshow(im, aspect='auto', extent=(0.9, 0.95, 3, 3.5))

	im = image.imread('C:/Users/Jeffrey/Documents/Academic/Stanford/Steg/Drafts/v9/taxa/reptile.png')
	ax.imshow(im, aspect='auto', extent=(0.87, 0.97, 4, 4.5))

	im = image.imread('C:/Users/Jeffrey/Documents/Academic/Stanford/Steg/Drafts/v9/taxa/fungi.png')
	ax.imshow(im, aspect='auto', extent=(0.9, 0.95, 1, 1.5))

	if allGraphs == 'no':
		plt.tight_layout()
		plt.savefig(outfile)
	else:
		ax.text(0.05, 7.6, 'C.', fontsize = 16, horizontalalignment='center')

		ax.text(0.85, 7.2, 'A', fontsize = 16, horizontalalignment='center')
		ax.text(0.85, 6.2, 'A', fontsize = 16, horizontalalignment='center')
		ax.text(0.85, 5.2, 'B', fontsize = 16, horizontalalignment='center')
		ax.text(0.85, 4.2, 'BC', fontsize = 16, horizontalalignment='center')
		ax.text(0.85, 3.2, 'CD', fontsize = 16, horizontalalignment='center')
		ax.text(0.85, 2.2, 'E', fontsize = 16, horizontalalignment='center')
		ax.text(0.85, 1.2, 'F', fontsize = 16, horizontalalignment='center')



		def singleTrialGraph(pixelBySpeciesPivot):
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
				if np.sum(data[x,6:-1]) > maxOccs:
					
					#Extract that line to play with
					subArray	=	np.reshape(data[x,6:-1], len(data[x,6:-1]))
					
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
					data[x,6:-1]	=	replacementRow	



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
			
			#Make randomized bell curves
			trials		=	25000
			aucValue	=	np.zeros(shape=(trials, 2))
			gridValue	=	np.zeros(shape=(trials, 2))
			regression	=	np.zeros(shape=(trials, 2))

			lowest	=	100000000000000
			highest	=	0
			holdingOld	=	np.copy(holding)	
			#Loop through for each permutation
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
				aucValue[w,1]	=	np.sum((numSpecies - modelSpecies)**2)				


			
				
				if aucValue[w,1] > highest:
					highest		=	aucValue[w,1]
					highestA	=	np.copy(holding)
					
				
				if aucValue[w,1] < lowest:
					lowest		=	aucValue[w,1]
					lowestA 	=	np.copy(holding)
					
				
				#Randomly shuffle by ecoregion								
				np.random.shuffle(holding)
				
				x	=	0
				y	=	0
				
				#Recalculate the number of species you would've found in each of the 
				for y in range(numUniques):
					
					x += holding[y,0]
					holding[y,1]	=	numSpecies[x]
				holding[-1,1]	=	maxSpp	

			randomB	=	np.copy(holding)
			return data, numSpecies, holdingOld, aucValue, npix, highestA, lowestA, randomB
			
			
			

		data, numSpecies, holdingOld, aucValue, npix, highestA, lowestA, randomB	=	singleTrialGraph('E:\\GIS_layers\\GBIF\\Plants\\Outputs\\7450a.csv')				

					
		ax3		=	fig1.add_subplot(gs[0,0])


		ax3.set_title('Sharp-Transition' + '\n' +  'Hypothesis', weight = 'bold', fontsize = 18)
		ax3.set_ylabel('Cumulative number' + '\n' +  'of species', fontsize = 16, color = 'k')
		ax3.set_xlabel('Area sampled', fontsize = 16, color = 'k')
		plt.plot(data[:,0],numSpecies, color = 'red')
		plt.plot([np.sum(holdingOld[:-1, 0]), npix], [(holdingOld[-1, 1]), (holdingOld[-1, 1])], color='k', linestyle='-', linewidth=2)

		ax3.set_yticklabels([])
		ax3.set_xticklabels([])
		x = 0
		for x in range(len(holdingOld)):
			if x == 0:
				plt.plot([0, holdingOld[0,0]], [holdingOld[0,1], holdingOld[0,1]], color='k', linestyle='-', linewidth=2)
			else:
				plt.plot([np.sum(holdingOld[:x-1, 0]), np.sum(holdingOld[:x, 0])], [(holdingOld[x-1, 1]), (holdingOld[x-1, 1])], color='k', linestyle='-', linewidth=2)
				plt.plot([np.sum(holdingOld[:x, 0]), np.sum(holdingOld[:x, 0])], [(holdingOld[x-1, 1]), (holdingOld[x, 1])], color='k', linestyle='-', linewidth=2)

		ymin, ymax = ax3.get_ylim()
		xmin, xmax = ax3.get_xlim()
		ax3.text(xmin, ymax + 3, 'A.', fontsize = 16, horizontalalignment='left')
		holding		=	np.copy(randomB)
		x = 0
		plt.plot([np.sum(holding[:-1, 0]), npix], [(holding[-1, 1]), (holding[-1, 1])], color='k', linestyle=':', linewidth=2)
		for x in range(len(holding)):
			if x == 0:
				plt.plot([0, holding[0,0]], [holding[0,1], holding[0,1]], color='g', linestyle=':', linewidth=2)
			else:
				plt.plot([np.sum(holding[:x-1, 0]), np.sum(holding[:x, 0])], [(holding[x-1, 1]), (holding[x-1, 1])], color='k', linestyle=':', linewidth=2)
				plt.plot([np.sum(holding[:x, 0]), np.sum(holding[:x, 0])], [(holding[x-1, 1]), (holding[x, 1])], color='k', linestyle=':', linewidth=2)





		ax4		=	fig1.add_subplot(gs[1,0])
		ax4.set_ylabel('Frequency', fontsize = 16)
		ax4.set_xlabel('Residual', fontsize = 16)
		ax4.set_title('B.', loc='left', fontsize = 16)
		sns.kdeplot(aucValue[:,1], color = 'pink', kernel = 'gau', shade = True)
		plt.axvline(x = aucValue[0,1], linestyle='-', color = 'k')
		plt.axvline(x = aucValue[-1,1], linestyle=':', color = 'k')
		ax4.set_yticklabels([])
		ax4.set_xticklabels([])

		data, numSpecies, holdingOld, aucValue, npix, highestA, lowestA, randomB	=	singleTrialGraph('E:\\GIS_layers\\GBIF\\Plants\\Outputs\\9660b.csv')				
		holding		=	np.copy(lowestA)
					
		ax5		=	fig1.add_subplot(gs[0,1])


		ax5.set_title('Gradual-Transition' + '\n' +  'Hypothesis', weight = 'bold', fontsize = 18)
		ax5.set_ylabel('Cumulative number' + '\n' +  'of species', fontsize = 16, color = 'k')
		ax5.set_xlabel('Area sampled', fontsize = 16, color = 'k')
		plt.plot(data[:,0],numSpecies, color = 'red')
		ax5.set_yticklabels([])
		ax5.set_xticklabels([])

		plt.plot([np.sum(holdingOld[:-1, 0]), npix], [(holdingOld[-1, 1]), (holdingOld[-1, 1])], color='k', linestyle='-', linewidth=2)
		x = 0
		for x in range(len(holdingOld)):
			if x == 0:
				plt.plot([0, holdingOld[0,0]], [holdingOld[0,1], holdingOld[0,1]], color='k', linestyle='-', linewidth=2)
			else:
				plt.plot([np.sum(holdingOld[:x-1, 0]), np.sum(holdingOld[:x, 0])], [(holdingOld[x-1, 1]), (holdingOld[x-1, 1])], color='k', linestyle='-', linewidth=2)
				plt.plot([np.sum(holdingOld[:x, 0]), np.sum(holdingOld[:x, 0])], [(holdingOld[x-1, 1]), (holdingOld[x, 1])], color='k', linestyle='-', linewidth=2)

		ymin, ymax = ax5.get_ylim()
		xmin, xmax = ax5.get_xlim()
		ax5.text(xmin, ymax + 3, 'D.', fontsize = 16, horizontalalignment='left')
		holding		=	np.copy(highestA)
		x = 0
		plt.plot([np.sum(holding[:-1, 0]), npix], [(holding[-1, 1]), (holding[-1, 1])], color='k', linestyle=':', linewidth=2)
		for x in range(len(holding)):
			if x == 0:
				plt.plot([0, holding[0,0]], [holding[0,1], holding[0,1]], color='g', linestyle='-', linewidth=2)
			else:
				plt.plot([np.sum(holding[:x-1, 0]), np.sum(holding[:x, 0])], [(holding[x-1, 1]), (holding[x-1, 1])], color='k', linestyle=':', linewidth=2)
				plt.plot([np.sum(holding[:x, 0]), np.sum(holding[:x, 0])], [(holding[x-1, 1]), (holding[x, 1])], color='k', linestyle=':', linewidth=2)



			

		ax6		=	fig1.add_subplot(gs[1,1])
		ax6.set_ylabel('Frequency', fontsize = 16)
		ax6.set_xlabel('Residual', fontsize = 16)
		ax6.set_title('E.', loc='left', fontsize = 16)
		sns.kdeplot(aucValue[:,1], color = 'pink', kernel = 'gau', shade = True)
		plt.axvline(x = aucValue[0,1], linestyle='-', color = 'k')
		plt.axvline(x = aucValue[-1,1], linestyle=':', color = 'k')
		ax6.set_yticklabels([])
		ax6.set_xticklabels([])
		
		plt.tight_layout()
		plt.savefig(outfile)
	plt.close()

makeGraph('summary2017', 75, 250, 40, 5, 'no', 'C:\\Users\\Jeffrey\\Documents\\Academic\\Stanford\\Steg\\Drafts\\Fig3.png')
#makeGraph('summary2001', 75, 250, 40, 5, 'no', 'C:\\Users\\Jeffrey\\Documents\\Academic\\Stanford\\Steg\\Drafts\\v14_revision2\\Figs\\FigS1.png')
#makeGraph('summarySub', 75, 250, 40, 5, 'no', 'C:\\Users\\Jeffrey\\Documents\\Academic\\Stanford\\Steg\\Drafts\\v14_revision2\\Figs\\FigS3.png')

#makeGraph('summary2017', 50, 100, 20, 4, 'no', 'C:\\Users\\Jeffrey\\Documents\\Academic\\Stanford\\Steg\\Drafts\\v14_revision2\\Figs\\FigS8.png')
#makeGraph('summary2017', 100, 500, 60, 8, 'no', 'C:\\Users\\Jeffrey\\Documents\\Academic\\Stanford\\Steg\\Drafts\\v14_revision2\\Figs\\FigS10.png')
