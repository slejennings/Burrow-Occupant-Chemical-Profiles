# GitHub Repository Read Me: Burrow & Occupant Chemical Profiles

## **Overview**

This repository contains data and code for the analyses in the manuscript “XXXXX”. The analyses are organized into an R Project. This Read Me file describes the required software and the organization of the R Project and the associated files.

## **Correspondence**

Please direct questions about the data, analysis, and results to:

XXXXX

## **Software** 

R programming language version 4.3.2 (2023-10-31)

RStudio IDE version 2024.12.0+467

R packages with version number:
broom (1.0.6),
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
In this R Project, the raw data files are manipulated and subseted to create the files that were imported and analyzed within PRIMER. We also provide copies of the saved outputs of the analyses performed in PRIMER to increase transparency. A user who has access to the PRIMER software would be able to reproduce or confirm our findings with these files.

## **Description of Folders and Files**

This repository contains an R project with various folders that are organized and named for their contents. A description of each folder and the files contained within each is provided below.

### **Scripts Folder**

***Description:*** contains R scripts to process the data, to generate some of findings in the manuscript, to prepare data files for import into the PRIMER software program, and to create figures. The scripts are sequential and are labeled accordingly (Step1 through Step 6).

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

    *File description:* compare chemical profiles of burrows (soil samples) and their avian occupants (feather samples) to determine whether there is evidence to suport a transfer of chemical fingerprints from the bird to the burrow and vice versa
  	
### **Data Folder**

***Description:*** contains original datafiles used for the analysis, stored as .csv

***Contents:*** 6 files. All .csv

1)	soilpeakarea.csv

    *File description:* peak areas of compounds that were measured and detected in the soil samples

    *Columns:*

     - File: unique identifier that contains multiple pieces of information. Bkgd indicates the sample came from the background (forest floor). Alternatively, Burr indicates the sample came from inside a burrow. The numbers that follow Bkgd or Burr (e.g., BurrXXX) denote the unique ID number associated with the marked burrow where the sample was collected. The information after the first underscore is the Site. The study area was divided into 3 seperate areas of marked storm-petrel burrows that are referred to as Sites (e.g., Site1, Site2, Site3). The final number after the second underscore is the replicate number (e.g., _1, _2, _3); as each sample was analyzed in triplicate, this number is either 1, 2 or 3.   


  	 - RT: the retetion time when the compound exited the GC (gas chromatograph) column and was detected by the MS (mass spectrometer). Each retention time identifies a unique chemical.


  	  - Area: peak area of each compound integrated from the chromatogram. The area of the peak reflects the abundance of the compound.
        

3)	soilsampleinfo.csv

    *File description:* additional identifying and measured information about each soil sample

    *Columns:*

     - File: unique identifier that contains multiple pieces of information. Bkgd indicates the sample came from the background (forest floor). Alternatively, Burr indicates the sample came from inside a burrow. The numbers that follow Bkgd or Burr (e.g., BurrXXX) denote the unique ID number associated with the marked burrow where the sample was collected. The information after the first underscore is the Site. The study area was divided into 3 seperate areas of marked storm-petrel burrows that are referred to as Sites (e.g., Site1, Site2, Site3). The final number after the second underscore is the replicate number (e.g., _1, _2, _3); as each sample was analyzed in triplicate, this number is either 1, 2 or 3.
  
       
     - Name:
  
       
     - Type: Bkgd or Burr. Identifies if the sample was taken inside the burrow (Burr) or from the forest floor outside the burrow (Bkgd).
  
       
     - Location:
  
       
     - Burrow: the unique ID number associated with the sampled burrow.
  
       
     - Site: Site1, Site2 or Site3. Marks the site within the larger study area where the sample was collected.
  
       
     - Feathers: Y or N. Denotes whether there is a feather sample(s) from occupant(s) of this burrow that will be used to examine the degree of overlap between the burrow and the bird.
  
       
     - Occupied: Y or N. Denotes whether the burrow was actively used by nesting birds during the study year.
  
       
     - Pair: Y or N. Indicates whether there are feather samples from both indivduals that occupied this burrow (aka a breeding pair).
  
       
     - IS_Area: the peak area of the internal standard compounds that was added to each sample (0.5mL of 25 mg/L d8-naphthalene in 100% ethanol).
  
       
     - Wet_Soil_Mass: mass of the sample that was used for chemical analysis in grams. Scent compounds were extracted from the headspace above soil samples that were thawed from frozen but still moist/wet. This measurement reflects the mass of the sample immediately prior to starting the extraction process.
  
       
     - WaterPercentage: the percentage of water in soil at each location. This value was determined by drying approximately 2 grams of soil in an oven for 24 hours, and comparing the initial (wet) mass with the resulting dry mass to obtain the percentage of the original sample that was water and had evaporated during the drying process.

4)	featherpeakarea.csv

    *File description:* peak areas of compounds that were measured and detected in the feather samples
  	
    *Columns:*
  	
    - File: unique identifier that contains multiple pieces of information. The first 12 numbers (e.g., XXXX-XXXXX) are taken from the bird's metal ID band. The number that comes after the underscore (e.g., _1, _2, _3) is the replicate number; as each sample was analyzed in triplicate, this number is either 1, 2 or 3.
  
      
    - Band_Number: 12-digit unique ID on the metal band placed on each bird.
  
      
    - RT: the retetion time when the compound exited the GC (gas chromatograph) column and was detected by the MS (mass spectrometer). Each retention time identifies a unique chemical.
  
      
     - Area: peak area of each compound integrated from the chromatogram. The area of the peak reflects the abundance of the compound.

5)	feathersampleinfo.csv
   
    *File description:* additional identifying and measured information about each feather sample
  	
    *Columns:*
  	
  	- File: unique identifier that contains multiple pieces of information. The first 12 numbers (e.g., XXXX-XXXXX) are taken from the bird's metal ID band. The number that comes after the underscore (e.g., _1, _2, _3) is the replicate number; as each sample was analyzed in triplicate, this number is either 1, 2 or 3.

    - Band_Number: 12-digit unique ID on the metal band placed on each bird.
  
      
    - Burrow: the unique ID number of the burrow each bird occupied during the sampling year.
  
      
    - Pair: Y or N. Denotes whether this individual is part of a breeding pair where feather samples were analyzed for both individuals.
  

    - BurrowOverlap: Y or N. Indicates whether the individual came from a burrow that had soil collected and analyzed.
  
      
    - IS_Area: the peak area of the internal standard compounds that was added to each sample (0.5mL of 10 mg/L d8-naphthalene in 100% ethanol).
  
      
    - Sample_Mass: mass of the sample that was used for chemical analysis in grams.

5)	Filename: burrow_coordinates.csv

    *File description:* geographic coordinates associated with each burrow (sampling location)

    *Columns:*
  	
     - Burrow: the unique ID number associated with the burrow.
  
       
     - Site: Site1, Site2 or Site3. Denotes the site within the larger study area where the burrow is located.
  
       
     - longitude: geographic longitude coordinate for the burrow's location in decimal degrees.
  
       
     - latitude: geographic latitude coordinate for the burrow's location in decimal degrees.

### **PRIMER Folder**

**Subfolder 1: Imported Files**

***Description:*** data frames created in R and exported as .xlxs that were imported into PRIMER, a multivariate software program, where certain analyses were performed

***Contents:*** 10 files. All .xlsx


**Subfolder 2: Exported Results**

***Description:*** the results of each model performed in PRIMER was exported and saved as .rtf to provide transparency about these findings beyond the values reported in the manuscript. Some additional results were also exported from the PRIMER as .xlsx files to allow for creation of figures in R. This included the scores (aka the x and y coordinates) of the samples identified by multivariate oridination methods in PRIMER, and the correlations between odor chemicals and the axes identified by multivariate models.

***Contents:*** 10 files. Saved as either .xlsx or .rtf


-------------------------
None of the other folders and files are required to reproduce our analysis. All of their contents can all be generated using the provided R scripts and data files. However, we have provided a brief description of what can be found in each of these folders.

### **Output Folder**

***Description:*** contains various files stored as .rds that contain organized data frames that were generated within the R scripts that are being moved between R scripts and to be used in subsequent analysis steps.

***Contents:*** 3 files. All .rds

### **Models Folder**

***Description:*** contains various models produced by the R scripts stored as .rds files. These files can be imported into R and viewed to see the results of our analyses.

***Contents:*** 5 files. All .rds

### **Tables and Figures Folder**

***Description:*** contains versions of the tables and figures in the manuscript that were generated by the scripts in the R Project

