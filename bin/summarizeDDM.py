


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
from collections import Counter
from matplotlib import cm
from matplotlib.gridspec import GridSpec
import ecopy as ep
from collections import Counter
import statsmodels.formula.api as smf
import matplotlib.image as image
from astropy.convolution import convolve
from astropy.convolution.kernels import Gaussian2DKernel


def makeGraph(year, minOccPixels, minPoints, minSpecies, minEcoregions, toTest, allGraphs, outfile):
	#Import summary files

	def readData(directory, toTest, numFiles):
		

		filePath	=	str('E:\\GIS_layers\\GBIF\\' + str(directory) + str(year) +'_0.csv')
		data		=	pd.read_csv(filePath,  	delimiter = ',', usecols = (0,1,2,3,4,5,6,8,9,10,11), names = ['identity', 'numSpecies', 'numPoints', 'numPixels', 'occPixels', 'comps', 'numEcoregions', 'log', 'sqrt', 'linear', 'uncorrected'])
		
		if numFiles > 1:
			for x in range(numFiles-1):
				y			=	x + 1
				filePath	=	str('E:\\GIS_layers\\GBIF\\' + str(directory) + '\\ddm_summary_2017_jaccard_' + str(y) + '.csv')
				data		=	data.append(pd.read_csv(filePath,  	delimiter = ',', usecols = (0,1,2,3,4,5,6,8,9,10,11), names = ['identity', 'numSpecies', 'numPoints', 'numPixels', 'occPixels', 'comps', 'numEcoregions', 'identity', 'log', 'sqrt', 'linear', 'uncorrected']))
		
		data	=	data.drop_duplicates(subset=['identity'])
		
		


		ND_value		=	-9999
		minComps		=	0

		data	=	data[data['identity'] >= 0]

		data	=	data[data['occPixels'] 		>=minOccPixels]
		data	=	data[data['numSpecies']	 	>=minSpecies]
		data	=	data[data['numEcoregions'] 	>=minEcoregions]
		data	=	data[data['numPoints']		>=minPoints]
		data	=	data[data['comps'] 			>=minComps]
		

		if allGraphs == 'yes':
			print float(len(data[data[toTest] >= 0.95]))	/float(len(data)), directory
		return data.copy()
		
	toTest	=	'log'
	plants		=	readData('Plants', toTest, 4)
	bugs		=	readData('Arthropods', toTest, 2)
	mammals		=	readData('Mammals', toTest, 2)
	herps		=	readData('Herps', toTest, 2)
	birds		=	readData('Birds', toTest, 4)
	fungi		=	readData('Fungi', toTest, 2)
	reptiles	=	readData('Reptiles', toTest, 2)
	amphibians	=	readData('Amphibians', toTest, 2)


	taxa 	= 	[plants, bugs, herps, mammals, birds, fungi, reptiles, amphibians]
	taxaS	=	['plants', 'bugs', 'herps', 'mammals', 'birds', 'fungi', 'reptiles', 'amphibians']

	if allGraphs == 'yes':
		for x in range(len(taxa)):
			statist, p = scipy.stats.kstest(taxa[x][toTest], 'uniform')
			print 'ddm', str(taxaS[x]), 'uniform', statist, p
				
			for y in range(len(taxa)):
				statist, p = scipy.stats.ks_2samp(taxa[x][toTest], taxa[y][toTest])

				print str(taxaS[x]), str(taxaS[y]), statist, p
				y += 1
			x += 1
			y = 0
				

	#Define plot style
	plt.style.use(u'seaborn-white')
	fig2 	= 	plt.figure()
	gs	= GridSpec(2,4)
	
	if allGraphs == 'yes':
		fig2.set_size_inches(16,8)
		ax		=	fig2.add_subplot(gs[:,1:3])
		
	else:
		fig2.set_size_inches(10,8)
		ax		=	fig2.add_subplot(gs[:,:])


	#Plot transect specific p-values from Species accumulation curves
	data	=	[1-fungi[toTest], 1-bugs[toTest], 1-reptiles[toTest], 1-amphibians[toTest],  1-plants[toTest], 1-birds[toTest], 1-mammals[toTest]]


	plt.violinplot(data, showmeans=False, showextrema=False, showmedians=True, vert = False)
	ax.set_yticks([7,6,5,4,3,2,1])


	ax.text(0.5, 5.25, 'Plants ' + '(n = ' + str(len(plants)) + ')', fontsize = 16, horizontalalignment='center')
	ax.text(0.5, 2.25, 'Arthropods ' + '(n = ' + str(len(bugs)) + ')', fontsize = 16, horizontalalignment='center')
	ax.text(0.5, 3.25, 'Reptiles ' + '(n = ' + str(len(reptiles)) + ')', fontsize = 16, horizontalalignment='center')
	ax.text(0.5, 4.25, 'Amphibians ' + '(n = ' + str(len(amphibians)) + ')', fontsize = 16, horizontalalignment='center')
	ax.text(0.5, 7.25, 'Mammals ' + '(n = ' + str(len(mammals)) + ')', fontsize = 16, horizontalalignment='center')
	ax.text(0.5, 6.25, 'Birds ' + '(n = ' + str(len(birds)) + ')', fontsize = 16, horizontalalignment='center')
	ax.text(0.5, 1.25, 'Fungi ' + '(n = ' + str(len(fungi)) + ')', fontsize = 16, horizontalalignment='center')

	#ax.set_yticklabels(label, rotation = 45, ha = 'right')
	ax.set(xlim=(0, 1))
	ax.set(ylim=(0.5,7.5))
	ax.set_title('Distance-Similarity Matrices', weight = 'bold', fontsize = 22)
	ax.set_xlabel(	u"\u2190" + '  More supportive of Sharp-Transition Hypothesis', fontsize = 16)
	ax.set_yticklabels([])

	im = image.imread('C:/Users/Jeffrey/Documents/Academic/Stanford/Steg/Drafts/v9/taxa/plant.tif')
	ax.imshow(im, aspect='auto', extent=(0.9, 0.95,5, 5.5))

	im = image.imread('C:/Users/Jeffrey/Documents/Academic/Stanford/Steg/Drafts/v9/taxa/bug.png')
	ax.imshow(im, aspect='auto', extent=(0.89, 0.96, 2, 2.5))

	im = image.imread('C:/Users/Jeffrey/Documents/Academic/Stanford/Steg/Drafts/v9/taxa/amphibian.tif')
	ax.imshow(im, aspect='auto', extent=(0.87, 0.97, 4, 4.5))

	im = image.imread('C:/Users/Jeffrey/Documents/Academic/Stanford/Steg/Drafts/v9/taxa/mammal.tif')
	ax.imshow(im, aspect='auto', extent=(0.87, 0.97, 7, 7.5))

	im = image.imread('C:/Users/Jeffrey/Documents/Academic/Stanford/Steg/Drafts/v9/taxa/bird.tif')
	ax.imshow(im, aspect='auto', extent=(0.9, 0.95, 6, 6.5))

	im = image.imread('C:/Users/Jeffrey/Documents/Academic/Stanford/Steg/Drafts/v9/taxa/reptile.png')
	ax.imshow(im, aspect='auto', extent=(0.87, 0.97, 3, 3.5))

	im = image.imread('C:/Users/Jeffrey/Documents/Academic/Stanford/Steg/Drafts/v9/taxa/fungi.png')
	ax.imshow(im, aspect='auto', extent=(0.9, 0.95, 1, 1.5))

	if allGraphs == 'no':
		plt.tight_layout()
		plt.savefig(outfile)
	
	else:
		ax.text(0.05, 7.6, 'C.', fontsize = 16, horizontalalignment='center')

		ax.text(0.85, 7.2, 'A', fontsize = 16, horizontalalignment='center')
		ax.text(0.85, 6.2, 'A', fontsize = 16, horizontalalignment='center')
		ax.text(0.85, 5.2, 'A', fontsize = 16, horizontalalignment='center')
		ax.text(0.85, 4.2, 'B', fontsize = 16, horizontalalignment='center')
		ax.text(0.85, 3.2, 'B', fontsize = 16, horizontalalignment='center')
		ax.text(0.85, 2.2, 'C', fontsize = 16, horizontalalignment='center')
		ax.text(0.85, 1.2, 'D', fontsize = 16, horizontalalignment='center')




		def distanceMatrix(pixelBySpeciesPivot):	
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

				return data

			def setUpHoldingArray(data):
				#Set up for randomly shuffling while mantaining ecorgion bin width
				unique		=	list(np.unique(data[:,4]))
				numUniques	=	len(unique)
				holding		=	np.empty(len(unique))
				npix		=	len(data)
				
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
				holdingOld	=	np.copy(holding)
			
				return holding, unique, numUniques, npix, holdingOld
				
				
			def calculatePredictedGrid(data, holding):
				#Calculate real and stattistical differences between pixels
				X_true	=	data[:,6:-1]
				statDistances	=	scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(X_true, metric='jaccard', p=2, w=None, V=None, VI=None)).astype('float')
				statDistances	=	np.subtract(1, statDistances)
				
				
				
				realDistances	=	scipy.spatial.distance.cdist(data[:,2:3], data[:,2:3], metric='euclidean', p=2, V=None, VI=None, w=None).astype('float')


				
				x = 0
				for x in range(len(X_true)):
					statDistances[x,x] = np.nan
					realDistances[x,x] = np.nan
					if np.nansum(X_true[x,:]) == 0:
						statDistances[x,:] = np.nan
						statDistances[:,x] = np.nan
						
				
				
				realDistances[statDistances == 0] = 'nan'
				statDistances[statDistances == 0] = 'nan'
				statDistances	=	np.log(statDistances)
				statDistances 	= 	np.divide(np.subtract(statDistances, np.nanmean(statDistances)), np.nanstd(statDistances))	


				
				realDistances	=	np.sqrt(realDistances)
				realDistances 	= 	np.multiply(np.divide(np.subtract(realDistances, np.nanmean(realDistances)), np.nanstd(realDistances)), -1)		
					

				
				flatStatA 	= 	np.ndarray.flatten(statDistances)
				flatStat 	= 	flatStatA[~np.isnan(flatStatA)]	

				
				distFlat	= 	np.ndarray.flatten(realDistances)	
				distFlat 	= 	distFlat[~np.isnan(flatStatA)]

				
				
				
				
				d = {'distance': distFlat, 'jaccard': flatStat}
				df = pd.DataFrame(data=d)
				


				nick	=	smf.ols(formula = 'jaccard~distance', data = df)
				res		=	nick.fit()

				intercept		=	res.params[0]
				distValue		=	res.params[1]

				
				
				predictedGrid	=	np.add(np.multiply(realDistances, distValue), intercept)
				observedMinusPredicted	=	np.subtract(statDistances,  predictedGrid)
				


				graphingGrid	=	np.zeros_like(observedMinusPredicted)
				x_left	=	0
				y_left	=	0
				x_right	=	0
				y_right	=	0
				
				for i in range(len(holding)):
					x_left	=	int(np.sum(holding[:i]))
					x_right	=	int(np.sum(holding[:i+1]))
					
					for j in range(len(holding)):
						y_left	=	int(np.sum(holding[:j]))
						y_right	=	int(np.sum(holding[:j+1]))
						
						z 	=	0
						
						try:
							z	=	np.nansum(observedMinusPredicted[x_left:x_right, y_left:y_right])
						except:
							z 	= 	0
						
						if np.isfinite(z) == True:			
							graphingGrid[x_left:x_right, y_left:y_right]	=	np.nanmean(observedMinusPredicted[x_left:x_right, y_left:y_right])
						
						j += 1
					j = 0
					i += 1
				

				
				return observedMinusPredicted, X_true, graphingGrid
			

			def doPermuatations(observedMinusPredicted, holding, trials): 
				gridValues	=	np.zeros(trials)
				highest		=	-100000
				lowest		=	100000
				
				
				for i in range(trials):

					#Do modularity calculation
					a = 0
					b = 0
					j = 0
					
					
					
					for j in range(len(unique)):
						a 	=	int(np.sum(holding[:j]))
						b	=	int(np.sum(holding[:j+1]))
						z 	=	0
						
						
						
						try:
							z	=	np.nansum(observedMinusPredicted[a:b,a:b])
						except:
							continue
						
						
						
						if np.isfinite(z) == True:
							gridValues[i]	+= z
						
						

						
						a	+=	b
					
					if gridValues[i] >= highest:
						highest		=	gridValues[i]
						highestA	=	np.copy(holding)
					
					elif gridValues[i] <= lowest:
						lowest		=	gridValues[i]
						lowestA		=	np.copy(holding)
					
					np.random.shuffle(holding)
					 
				randomB	=	np.copy(holding)
				return gridValues, highestA, lowestA, randomB
				
				

			#Read in species pivot table data
			dataO 	=	np.genfromtxt(pixelBySpeciesPivot, dtype= 'float' , delimiter=',', skip_header=1)

			#Change non-occurences to zeros (important for futrue analysis)
			data	=	np.nan_to_num(dataO)

			#Filter out extra occurances and pixels with too few occurences
			data	=	removeMinOccs(data, 0)
			data	=	reduceMaxOccs(data, 50)

			#Figure out where ecoregion boundaries fall
			holding, unique, numUniques, npix, holdingOld	=	setUpHoldingArray(data)
			
			
			observedMinusPredicted, X_true, graphingGrid	=	calculatePredictedGrid(data, holding)


			trials 	=	5000
			gridValues, highestA, lowestA, randomB	=	doPermuatations(observedMinusPredicted, holding, trials)

			
			return gridValues, observedMinusPredicted, holdingOld, graphingGrid, highestA, lowestA, randomB
			
		gridValues, observedMinusPredicted, holdingOld, graphingGrid, highestA, lowestA, randomB		=		distanceMatrix('E:\\GIS_layers\\GBIF\\Plants\\Outputs\\9660b.csv')




		ax7		=	fig2.add_subplot(gs[0,0])
		ax7.set_title('Sharp-Transition' + '\n' + 'Hypothesis', weight = 'bold', fontsize = 18)
		ax7.set_xlabel('Pixel Number', fontsize = 16, color = 'k')
		ax7.set_ylabel('Pixel Number', fontsize = 16, color = 'k')
		ax7.text(1, -30, 'A.', fontsize = 16, horizontalalignment='left')

		setMin	=	(max(np.nanmin(observedMinusPredicted) * -1, np.nanmax(observedMinusPredicted)) * -1)/2
		setMax	=	(max(np.nanmin(observedMinusPredicted) * -1, np.nanmax(observedMinusPredicted)))/2



		plt.imshow(convolve(observedMinusPredicted, Gaussian2DKernel(stddev=8)), interpolation='none', cmap = 'bwr', vmin = setMin, vmax = setMax)



		toPlot	=	np.zeros_like(holdingOld)
		for i in range(len(holdingOld)):
			toPlot[i]	=	np.sum(holdingOld[:i])
			
		for xc in range(len(holdingOld)):
			plt.axvline(x=toPlot[xc], linestyle='-', color = 'k')
			plt.axhline(y=toPlot[xc], linestyle='-', color = 'k')	
			
		toPlot	=	np.zeros_like(randomB)
		for i in range(len(randomB)):
			toPlot[i]	=	np.sum(randomB[:i])
			
		for xc in range(len(holdingOld)):
			plt.axvline(x=toPlot[xc], linestyle=':', color = 'k')
			plt.axhline(y=toPlot[xc], linestyle=':', color = 'k')	
				

				

		ax8		=	fig2.add_subplot(gs[1,0])
		ax8.set_ylabel('Frequency', fontsize = 16)
		ax8.set_xlabel('Modularity explained', fontsize = 16)
		ax8.set_title('B.', loc='left', fontsize = 16)
		sns.kdeplot(gridValues, color = 'pink', kernel = 'gau', shade = True)
		plt.axvline(x = gridValues[0], linestyle='-', color = 'k')
		plt.axvline(x = gridValues[-1], linestyle=':', color = 'k')
		ax8.set_yticklabels([])
		ax8.set_xticklabels([])	
				

				
				
				
				
				
		gridValues, observedMinusPredicted, holdingOld, graphingGrid, highestA, lowestA, randomB		=		distanceMatrix('E:\\GIS_layers\\GBIF\\Plants\\Outputs\\1321a.csv')
				



		ax9		=	fig2.add_subplot(gs[0,1])
		ax9.set_title('Gradual-Transition' + '\n' + 'Hypothesis', weight = 'bold', fontsize = 18)
		ax9.set_xlabel('Pixel Number', fontsize = 16, color = 'k')
		ax9.set_ylabel('Pixel Number', fontsize = 16, color = 'k')
		ax9.text(1, -30, 'D.', fontsize = 16, horizontalalignment='left')


		setMin	=	(max(np.nanmin(observedMinusPredicted) * -1, np.nanmax(observedMinusPredicted)) * -1)/2
		setMax	=	(max(np.nanmin(observedMinusPredicted) * -1, np.nanmax(observedMinusPredicted)))/2


		plt.imshow(convolve(observedMinusPredicted, Gaussian2DKernel(stddev=8)), interpolation='none', cmap = 'bwr', vmin = setMin, vmax = setMax)



		toPlot	=	np.zeros_like(holdingOld)
		for i in range(len(holdingOld)):
			toPlot[i]	=	np.sum(holdingOld[:i])
			
		for xc in range(len(holdingOld)):
			plt.axvline(x=toPlot[xc], linestyle='-', color = 'k')
			plt.axhline(y=toPlot[xc], linestyle='-', color = 'k')	
			
		toPlot	=	np.zeros_like(randomB)
		for i in range(len(randomB)):
			toPlot[i]	=	np.sum(randomB[:i])
			
		for xc in range(len(holdingOld)):
			plt.axvline(x=toPlot[xc], linestyle=':', color = 'k')
			plt.axhline(y=toPlot[xc], linestyle=':', color = 'k')	
				

			


		ax6		=	fig2.add_subplot(gs[1,1])
		ax6.set_ylabel('Frequency', fontsize = 16)
		ax6.set_xlabel('Modularity explained', fontsize = 16)
		ax6.set_title('E.', loc='left', fontsize = 16)
		sns.kdeplot(gridValues, color = 'pink', kernel = 'gau', shade = True)
		plt.axvline(x = gridValues[0], linestyle='-', color = 'k')
		plt.axvline(x = gridValues[-1], linestyle=':', color = 'k')
		ax6.set_yticklabels([])
		ax6.set_xticklabels([])	


			


				
		#Show the plot to fix formatting before saving (click on tight layout)
		plt.savefig(outfile)
	
	plt.close()
makeGraph('\\ddm_summary_2017_jaccard', 75, 250, 40, 5, 'log', 'no', 'C:\\Users\\Jeffrey\\Documents\\Academic\\Stanford\\Steg\\Drafts\\v14_revision2\\Figs\\Fig3.png')

makeGraph('\\ddm_summary_subset_jaccard', 75, 250, 40, 5, 'log', 'no', 'C:\\Users\\Jeffrey\\Documents\\Academic\\Stanford\\Steg\\Drafts\\v14_revision2\\Figs\\FigS2.png')

makeGraph('\\ddm_summary_2001_jaccard', 75, 250, 40, 5, 'log', 'no', 'C:\\Users\\Jeffrey\\Documents\\Academic\\Stanford\\Steg\\Drafts\\v14_revision2\\Figs\\FigS4.png')




makeGraph('\\ddm_summary_2017_jaccard', 75, 250, 40, 5, 'sqrt', 'no', 'C:\\Users\\Jeffrey\\Documents\\Academic\\Stanford\\Steg\\Drafts\\v14_revision2\\Figs\\FigS5.png')
makeGraph('\\ddm_summary_2017_jaccard', 75, 250, 40, 5, 'linear', 'no', 'C:\\Users\\Jeffrey\\Documents\\Academic\\Stanford\\Steg\\Drafts\\v14_revision2\\Figs\\FigS6.png')


makeGraph('\\ddm_summary_2017_chao', 75, 250, 40, 5, 'log', 'no', 'C:\\Users\\Jeffrey\\Documents\\Academic\\Stanford\\Steg\\Drafts\\v14_revision2\\Figs\\FigS7.png')


makeGraph('\\ddm_summary_2017_jaccard', 50, 100, 20, 4, 'log', 'no', 'C:\\Users\\Jeffrey\\Documents\\Academic\\Stanford\\Steg\\Drafts\\v14_revision2\\Figs\\FigS9.png')
makeGraph('\\ddm_summary_2017_jaccard', 100, 500, 60, 8, 'log', 'no', 'C:\\Users\\Jeffrey\\Documents\\Academic\\Stanford\\Steg\\Drafts\\v14_revision2\\Figs\\FigS11.png')
