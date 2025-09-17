# Read Me: Bird-scented nests as a mechanism for olfactory homing in a burrow nesting seabird

## **Overview**

This repository contains data and code for the analyses in the manuscript “Bird-scented nests as a mechanism for olfactory homing in a burrow nesting seabird”. 
The analyses are organized into an R Project. This Read Me file describes the required software and the organization of the R Project and the associated files.

Author: Sarah L. Jennings

Other contributors: Susan E. Ebeler, Gail L. Patricelli

Dataset title: Data and code for the article "Bird-scented nests as a mechanism for olfactory homing in a burrow nesting seabird"

Persistent Identifier: https://doi.org/10.5281/zenodo.16878922

Date created: 03/28/2025

Dataset citation: Jennings, S.L., G.L. Patricelli, S.E. Ebeler. 2025. Data and code for "Bird-scented nests as a mechanism for olfactory homing in a burrow-nesting seabird". Zenodo. https://doi.org/10.5281/zenodo.16878922

Corresponding publication: Jennings, S.L., G.L. Patricelli, S.E. Ebeler. 2025. Bird-scented nests as a mechanism for olfactory homing in a burrow-nesting seabird. The American Naturalist

## **Correspondence**

Please direct questions to:
 
Name: Sarah L. Jennings, PhD

Affliations: California State Polytechnic University, San Luis Obispo

ORCID: 0000-0003-2995-8813

Email: sjenni02@calpoly.edu

## **Metadata**

Dates and Locations: fieldwork completed in July, 2015 on Bon Portage Island, Nova Scotia, Canada. Laboratory and data analyses conducted at the University of California, Davis, USA.

Methodological information: provided in the corresponding manuscript

## **Software** 

R programming language version 4.3.2 (2023-10-31)

RStudio IDE version 2024.12.0+467

R packages with version number:
broom (1.0.6),
colorspace (2.1-1),
ecodist (2.3.1),
geodist (0.0.8),
here (1.0.1),
patchwork (1.2.0),
perm (1.0-0.4),
RColorBrewer (1.1-3),
tidyverse (2.0.0),
vegan (2.6-4),
writexl (1.4.2)

**Other software:** 

PRIMER 7 with PERMANOVA+ (https://www.primer-e.com/software)
This is a paid software program. It offers a wide range of methods for multivariate data and offers several advantages to the options that are currently available in free, open source analysis packages in R. One of the key differences is that the methods used in PRIMER to partition the variation in multivariate data are better able to handle unbalanced or complex designs compared with R, which is why we have opted to use it for this study. For more information see: https://learninghub.primer-e.com/books/should-i-use-primer-or-r

In this R Project, the raw data files are manipulated and subsetted to create the files that were imported and analyzed within PRIMER. We also provide copies of the saved outputs of the analyses performed in PRIMER to increase transparency. A user who has access to the PRIMER software would be able to reproduce or confirm our findings with these files.

## **Description of Folders and Files**

This repository contains an R project with various folders that are organized and named for their contents. A description of each folder and the files contained within each is provided below.

***Overview:***
* Number of Folders: 8
* Number of Files: 52
* File formats: .cvs, .xlsx, .R, .rtf, .png

***Workflow***

The R scripts are numbered in Steps (1 through 6), which is noted as part of their file name. They should be run in numerical order starting with Soil_Step1. R and ending with Feathers_Step6.R

-------------------------
### **Scripts Folder**

***Description:*** contains R scripts to process the data, to generate some of findings in the manuscript, to prepare data files for import into the PRIMER software program, and to create figures. The scripts are sequential and are labeled accordingly (Step 1 through Step 6).

***Contents:*** 6 files. All .R

1)	Soil_Step1.R

    *File description:* organize and process chemical data from soil samples. Examine the repeatability across triplicate samples
  	

3)	Soil_Step2.R

    *File description:* perform qualitative and quantitative comparisons of the chemical profiles of the different types of soil samples
  	

5)	Soil_Step3.R

    *File description:* examine whether soil sample chemical profiles vary across the landscape. Create multiple files for export to PRIMER and generate figures from the results of the PRIMER analyses
  	

4)	Feathers_Step4.R

    *File description:* organize and process chemical data from feather samples. Examine the repeatability across triplicate samples
  	

6)	Feathers_Step5.R

    *File description:* examine the similarities and differences in the chemical profiles of birds that are part of a breeding pair that shares a burrow
  	

8)	Feathers_Step6.R

    *File description:* compare chemical profiles of burrows (soil samples) and their avian occupants (feather samples) to determine whether there is evidence to support a transfer of chemical fingerprints from the bird to the burrow and vice versa

-------------------------
### **Data Folder**

***Description:*** contains original datafiles used for the analysis, stored as .csv

***Contents:*** 6 files. All .csv

1)	soilpeakarea.csv

    *File description:* peak areas of compounds that were measured and detected in the soil samples

    *Columns:*

     - File: unique identifier that contains multiple pieces of information. Bkgd indicates the sample came from the background (forest floor). Alternatively, Burr indicates the sample came from inside a burrow. The numbers that follow Bkgd or Burr (e.g., BurrXXX) denote the unique ID number associated with the marked burrow where the sample was collected. The information after the first underscore is the Site. The study area was divided into 3 seperate areas of marked storm-petrel burrows that are referred to as Sites (e.g., Site1, Site2, Site3). The final number after the second underscore is the replicate number (e.g., _1, _2, _3); as each sample was analyzed in triplicate, this number is either 1, 2 or 3.   


  	 - RT: the time in minutes when the compound exited the GC (gas chromatograph) column and was detected by the MS (mass spectrometer). Each retention time identifies a unique chemical.


  	  - Area: peak area of each compound integrated from the chromatogram. The area of the peak reflects the abundance of the compound.
        

3)	soilsampleinfo.csv

    *File description:* additional identifying and measured information about each soil sample

    *Columns:*

     - File: unique identifier that contains multiple pieces of information. Bkgd indicates the sample came from the background (forest floor). Alternatively, Burr indicates the sample came from inside a burrow. The numbers that follow Bkgd or Burr (e.g., BurrXXX) denote the unique ID number associated with the marked burrow where the sample was collected. The information after the first underscore is the Site. The study area was divided into 3 seperate areas of marked storm-petrel burrows that are referred to as Sites (e.g., Site1, Site2, Site3). The final number after the second underscore is the replicate number (e.g., _1, _2, _3); as each sample was analyzed in triplicate, this number is either 1, 2 or 3.
       
     - Type: Bkgd or Burr. Identifies if the sample was taken inside the burrow (Burr) or from the forest floor outside the burrow (Bkgd).
  
       
     - Location: combines the columns Type (Bkgd or Burr) and Burrow (unique ID number for each sampled burrow)
  
       
     - Burrow: the unique ID number associated with the sampled burrow.
  
       
     - Site: Site1, Site2 or Site3. Marks the site within the larger study area where the sample was collected. See corresponding manuscript for a map.
  
       
     - Feathers: Y or N. Denotes whether there is a feather sample(s) from occupant(s) of this burrow that will be used to examine the degree of overlap between the burrow and the bird.
  
       
     - Occupied: Y or N. Denotes whether the burrow was actively used by nesting birds during the study year.
  
       
     - Pair: Y or N. Indicates whether there are feather samples from both indivduals that occupied this burrow (aka a breeding pair).
  
       
     - IS_Area: the peak area of the internal standard compound that was added to each sample (for soil samples this compound was 0.5mL of 25 mg/L d8-naphthalene in 100% ethanol).
  
       
     - Wet_Soil_Mass: mass of the sample that was used for chemical analysis in grams. Scent compounds were extracted from the headspace above soil samples that were thawed from frozen but still moist/wet. This measurement reflects the mass of the sample immediately prior to starting the extraction process.
  
       
     - WaterPercentage: the percentage of water in soil at each location. This value was determined by drying approximately 2 grams of soil in an oven for 24 hours, and comparing the initial (wet) mass with the resulting dry mass to obtain the percentage of the original sample that was water and had evaporated during the drying process.

4)	featherpeakarea.csv

    *File description:* peak areas of compounds that were measured and detected in the feather samples
  	
    *Columns:*
  	
    - File: unique identifier that contains multiple pieces of information. The first 12 numbers (e.g., XXXX-XXXXX) are taken from the bird's metal ID band. The number that comes after the underscore (e.g., _1, _2, _3) is the replicate number; as each sample was analyzed in triplicate, this number is either 1, 2 or 3.
  
      
    - Band_Number: 12-digit unique ID on the metal band placed on each bird.
  
      
    - RT: the time in minutes when the compound exited the GC (gas chromatograph) column and was detected by the MS (mass spectrometer). Each retention time identifies a unique chemical.
  
      
     - Area: peak area of each compound integrated from the chromatogram. The area of the peak reflects the abundance of the compound.


5)	feathersampleinfo.csv
   
    *File description:* additional identifying and measured information about each feather sample
  	
    *Columns:*
  	
  	- File: unique identifier that contains multiple pieces of information. The first 12 numbers (e.g., XXXX-XXXXX) are taken from the bird's metal ID band. The number that comes after the underscore (e.g., _1, _2, _3) is the replicate number; as each sample was analyzed in triplicate, this number is either 1, 2 or 3.

    - Band_Number: 12-digit unique ID on the metal band placed on each bird.
  
      
    - Burrow: the unique ID number of the burrow each bird occupied during the sampling year.
  
      
    - Pair: Y or N. Denotes whether this individual is part of a breeding pair where feather samples were analyzed for both individuals.
  

    - BurrowOverlap: Y or N. Indicates whether the individual came from a burrow that had soil collected and analyzed.
  
      
    - IS_Area: the peak area of the internal standard compound that was added to each sample (for soil samples this compound was 0.5mL of 10 mg/L d8-naphthalene in 100% ethanol).
  
      
    - Sample_Mass: mass of the sample that was used for chemical analysis in grams.


5)	Filename: burrow_coordinates.csv

    *File description:* geographic coordinates associated with each burrow (sampling location)

    *Columns:*
  	
     - Burrow: the unique ID number associated with the burrow.
  
       
     - Site: Site1, Site2 or Site3. Denotes the site within the larger study area where the burrow is located.
  
       
     - longitude: geographic longitude coordinate for the burrow's location in decimal degrees.
  
       
     - latitude: geographic latitude coordinate for the burrow's location in decimal degrees.

  
     - Sampling_Day: denotes the day of sampling for each burrow as 1 through 4.
       
  
6)	Filename: feathercompoundlist.csv

    *File description:* list of 155 compounds detected in feather samples

    *Columns:*

      - FeatherRT: the time in minutes when the compound exited the GC (gas chromatograph) column and was detected by the MS (mass spectrometer) in the feather samples. Each retention time identifies a unique chemical. While the same instrument and analysis program were used for both soil and feather samples, we needed to trim a small amount of the front end of the GC column as part of regular maintenance that occurred in between analyzing the two types of samples. This resulted in slight differences retention times between the two sample types such that the same chemical has a slightly earlier retention time in the soil data compared with the feather data. This column along with the SoilRT column help align the two types of data.
  
        
      - Compound_Name: the name of the chemical compound obtained through matching the retention time and mass spectra of the detected compound with mass spectral reference libraries.
  
        
      - Type: Bird, Plant or  Contaminant. Compounds were classified into three groups based on the most likely source of the chemical. These groups were used for analyses.
  
        
      - In_Burrow: Y or N. Indicates whether the compound was also detected in the burrow soil sample, and is therefore common to both sample types.
        
        
      - SoilRT: the time in minutes when the compound exited the GC (gas chromatograph) column and was detected by the MS (mass spectrometer) in the soil samples. Each retention time identifies a unique chemical. While the same instrument and analysis program were used for both soil and feather samples, we needed to trim a small amount of the front end of the GC column as part of regular maintenance that occurred in between analyzing the two types of samples. This resulted in slight differences retention times between the two sample types such that the same chemical has a slightly earlier retention time in the soil data compared with the feather data. This column along with the FeatherRT column help align the two types of data.

-------------------------
### **PRIMER Folder**

**Subfolder 1: Imported Files**

***Description:*** data frames created in R and exported as .xlxs that were imported into PRIMER, a multivariate software program, where certain analyses were performed

***Contents:*** 10 files. All .xlsx

1) AveStdPA_All.xlsx
 
      *File description:* file for import into PRIMER. Contains abundances (standardized and averaged peak areas) of all compounds in all soil samples (burrow and background).
   
2) AveStdPA_bkgd.xlsx

     *File description:* file for import into PRIMER. Contains abundances (standardized and averaged peak areas) of all compounds in background soil samples.
   
3) AveStdPA_burr.xlsx

    *File description:* file for import into PRIMER. Contains abundances (standardized and averaged peak areas) of all compounds in burrow soil samples.
   
4) burrowdistancematrix.xlsx

    *File description:* matrix of pairwise distances between sampled burrows. Distances are in meters.

5) Occburrows_elevbirdcomps.xlsx
 
   *File description:* the coordinates (scores) of each soil sample on the CAP1 and CAP2 axes that were created by the CAP model for Class. Used to create an ordination plot in R to visualize differences between samples.

6) Pairs_AllBirdComps.xlsx
 
   *File description:* file for import into PRIMER. Contains abundances (standardized and averaged peak areas) of all bird-derived compounds on feather samples of mated pairs.

7) Pairs_AllComps.xlsx

   *File description:* file for import into PRIMER. Contains abundances (standardized and averaged peak areas) of all compounds on feather samples of mated pairs.

8) Pairs_Contaminants.xlsx

   *File description:* file for import into PRIMER. Contains abundances (standardized and averaged peak areas) of synthetic contaminant compounds on feather samples of mated pairs.

9) Pairs_ElevBirdComps.xlsx

    *File description:* file for import into PRIMER. Contains abundances (standardized and averaged peak areas) of 6 elevated storm-petrel compounds on feather samples of mated pairs.

10) Pairs_PlantComps.xlsx

    *File description:* file for import into PRIMER. Contains abundances (standardized and averaged peak areas) of plant-derived compounds on feather samples of mated pairs.

--------------------------------------

**Subfolder 2: Exported Results**

***Description:*** the results of each model performed in PRIMER was exported and saved as .rtf to provide transparency about these findings beyond the values reported in the manuscript. Some additional results were also exported from the PRIMER as .xlsx files to allow for creation of figures in R. This included the scores (aka the x and y coordinates) of the samples identified by multivariate oridination methods in PRIMER, and the correlations between odor chemicals and the axes identified by multivariate models.

***Contents:*** 15 files. Saved as either .csv or .rtf

1) BVSTEP Bkgd Compounds.rft
 
      *File description:* output from PRIMER of BVSTEP model using background soil samples
   
2) BVSTEP Burrow Compounds.rtf

     *File description:* output from PRIMER of BVSTEP model using burrow soil samples
   
3) CAP for Site.rtf

    *File description:* output from PRIMER of Canonical Analysis of Principal Coordinates (CAP) model testing for differences between Sites of soil samples
   
4) CAP for Class.rft

    *File description:* output from PRIMER of Canonical Analysis of Principal Coordinates (CAP) model testing for differences between Classes of soil samples

5) CAP Scores Class.csv
 
   *File description:* the coordinates (scores) of each soil sample on the CAP1 and CAP2 axes that were created by the CAP model for Class. Used to create an ordination plot in R to visualize differences between samples.

6) CAP Scores Site.csv
 
   *File description:* the coordinates (scores) of each sample on the CAP1 and CAP2 axes that were created by the CAP model for Site. Used to create an ordination plot in R to visualize differences between samples.

7) Class CAP Axes Compound Correlations.csv

   *File description:* correlation scores between the soil chemicals and the CAP model axes (CAP1 and CAP2) created by the CAP model for Class

8) Site CAP Axes Compound Correlations.csv

   *File description:* correlation scores between the soil chemicals and the CAP model axes (CAP1 and CAP2) created by the CAP model for Site

9) Permanova Occupied Burrows Elevated Bird Compounds.rtf

    *File description:* output from PRIMER of PERMANOVA model testing for differences between soil chemical profiles of occupied burrows using the 6 chemicals derived from storm-petrels that were elevated in the occupied burrows

10) Permanova Pairs All Bird Compounds.rtf

    *File description:* output from PRIMER of PERMANOVA model testing for differences between feather chemical profiles of mated pairs using all bird-derived compounds

11) Permanova Pairs All Compounds.rtf

    *File description:* output from PRIMER of PERMANOVA model testing for differences between feather chemical profiles of mated pairs using all compounds

12) Permanova Pairs Contaminants.rtf

    *File description:* output from PRIMER of PERMANOVA model testing for differences between feather chemical profiles of mated pairs using compounds identified as synthetic contaminants

13) Permanova Pairs Elevated Bird Compounds.rtf

    *File description:* output from PRIMER of PERMANOVA model testing for differences between feather chemical profiles of mated pairs using 6 storm-petrel compounds that were elevated inside occupied burrows

14) Permanova Pairs Plant Compounds.rtf

    *File description:* output from PRIMER of PERMANOVA model testing for differences between feather chemical profiles of mated pairs using compounds identified as plant-derived

15) Permanova Site and Class.rtf

    *File description:* output from PRIMER of PERMANOVA model testing for the importance of Site and Sample Class in the soil samples chemical profiles


-------------------------
None of the other folders and files are required to reproduce our analysis. All of their contents can all be generated using the provided R scripts and data files. However, we have provided a brief description of what can be found in each of these folders.

### **Output Folder**

***Description:*** contains various files stored as .rds that contain organized data frames that were generated within the R scripts that are being moved between R scripts and to be used in subsequent analysis steps.

***Contents:*** 4 files. All .rds

### **Models Folder**

***Description:*** contains various models produced by the R scripts stored as .rds files. These files can be imported into R and viewed to see the results of our analyses.

***Contents:*** 4 files. All .rds

### **Figures Folder**

***Description:*** contains versions of the figures in the manuscript that were generated by the scripts in the R Project

***Contents:*** 7 files. All .png
