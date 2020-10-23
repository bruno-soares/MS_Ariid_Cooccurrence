# Data analysis for the manuscript: "Environmental conditions promote local segregation, but functional distinctiveness allow aggregation of catfishes in the Amazonian estuary"

## Authors
Bruno E. Soares, Naraiana L. Benone, Ronaldo B. Barthem, Alexandre P. Marceniuk and Luciano F. A. Montag

## Published in
Preprint available in EvoEvoRxiv: https://doi.org/10.32942/osf.io/nmyxt
Submitted to Estuarine, Coastal and Shelf Science

## To cite this data
Data and coding* available in Zenodo: https://zenodo.org/record/4118631

*But see/cite D'Amen et al. (2018) to the coding of C-score analyses, null modelling and Blois framework.

## Description
Data and script to analyze the local co-occurrence patterns of ariids in the Amazonian estuary. These results are discussed in the manuscript "Environmental conditions promote local segregation of catfishes in the Amazonian estuary, but ecomorphological differences may allow aggregation", in preparation.

## How to use this directory
In R, you may download all the data contained in this repository by using the download.file() function, then use the function unzip(). Example:

`
download.file(url="https://github.com/bruno-soares/MS_Ariid_Cooccurrence/archive/master.zip", destfile = "MS_Ariid_Cooccurrence.zip")
`

`
unzip(zipfile = "MS_Ariid_Cooccurrence.zip")
`

R scripts are numbered in the order that the analyses are performed, from 01 to 06. All the data necessary to run analyses are availagle in /data and /results, so one can run any part of the script independently from running prior analysis. Make sure you install the necessary R packages before running the analysis.

## Data description
Data encompasses the ecomorphological description of nine marine catfishes of the family Ariidae (Siluriformes) and their occurrence in the inner estuary of the Amazon River mouth. Data is available in three files:

1. ecomorphological data: seven ecomorphological traits from nine species, including the relative depth of the caudal peduncle, the relative area of the eyes, the relative length of the maxillary barbel, the relative length of the mentonian barbel, the relative width of the mouth, maximum size and shape of the teeth.

2. trawling data: information on the occurrence of ariids in the inner Amazonian estuary, including sampling date and season, latitude and longitude, depth at hauling (prof), maximum salinity in trawled area (sal), average temperature in trawled area (temp).

3. avg. salinity: information on the average salinity in trawled area (sal_avg) as complementary information to trawling data.

## Contact information
Please, feel free to contact me in my personal e-mail: soares.e.bruno@gmail.com.

## References
D'Amen, M; Mod, H. K.; Gotelli, N. J.; Guisan, A. 2018. Disentangling biotic interactions, environmental filters, and dispersal limitation as drivers of species co‐occurrence. Ecography, 41 (8): 1233-1244. doi: https://doi.org/10.1111/ecog.03148
