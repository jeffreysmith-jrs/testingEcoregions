# testingEcoregions
Code associated with the paper 'A global test of ecoregions' by Smith et al.

----------------
General instructions:

In order to run the code for Smith et al. (2018) you should first run codes 1 and 2. Then, to make the species-accumulation curve figures you should run codes 3a and 4a. For the distance-similarity matrix figures you should run codes 3b and 4b. The code in each case is dependent on the previous steps (e.g. code 2 will only work after running code 1).

----------------
File descriptions:

- 1. cleaningGBIF.py - This code takes raw input files from www.gbif.org and processes them into a csv containing the x and y pixel values of a specified basemap ascii and the species name.

- 2. makeTransects.py - This code takes the cleaned GBIF file from the above code and generates random transects around the globe following the methods laid out in Smith et al. (2018)

- 3a. sacTest.py - This code performs the species area accumulation curve test described in Smith et al. (2018), relying on the inputs generated from the makeTransects.py code

- 3b. ddmTest.py - This code performs the distance-similarity matrix test described in Smith et al. (2018), relying on the inputs generated from the makeTransects.py code

- 4a. summarizeSAC.py - This code summarizes the results of the sacTest.py code, allowing for the user to modify data quality control 

- 4b. summarizeDDM.py - This code summarizes the results of the ddmTest.py code, allowing for the user to modify data quality control 

----------------
Data download links:

Amphibians - http://api.gbif.org/v1/occurrence/download/request/0011811-180412121330197.zip

Arthropods - http://api.gbif.org/v1/occurrence/download/request/0011802-180412121330197.zip

Fungi - http://api.gbif.org/v1/occurrence/download/request/0011799-180412121330197.zip

Plants - http://api.gbif.org/v1/occurrence/download/request/0011797-180412121330197.zip

Reptiles - http://api.gbif.org/v1/occurrence/download/request/0011815-180412121330197.zip

Mammals - http://api.gbif.org/v1/occurrence/download/request/0011818-180412121330197.zip

Birds - http://api.gbif.org/v1/occurrence/download/request/0011822-180412121330197.zip


