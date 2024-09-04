# Two-ColorDiFC 
Processing and plotting DiFC data for:
Williams et.al. (2024) "Two-color diffuse in vivo flow cytometer," **Journal of Biomedical Optics** 29(6), 065003.

## Table of Contents
* [General Info](#general-info)
* [Technologies](#technologies)
* [Processing Data](#processing-data)
* [Preprocessed Data](#preprocessed-data)
* [Data analysis and plotting](#data-analysis-and-plotting)

## General Info
This files in this github repository were used to process, analyze and plot data used in the 2024 Journal of Biomedical Optics article "Two-color diffuse in vivo flow cytometer". The scripts and functions used for creating the figures in the article are included in this repository.

The data files exceed the GitHub repository size limit, so they are located in a Pennsieve repository: DOI: 10.26275/q2w0-keol.

## Technologies
Matlab 2020b

## Processing Data
To process the raw data download:
* ProcessAllData.m
* Two_Color_DiFC_Processing.m
* Pennsieve repository folder "Two-ColorDiFC_RawData"
* Folder "proccodes"
* remove4ChannelCoincPeaks.m
* match2ColorPeaks.m

## Preprocessed Data
To skip processing the raw data files, download the folder:
* Pennsieve repository folder "Two-ColorDiFC_Data"

Any downloaded folders should be located in your working directory and/or paths to the folders should be manually changed at the top of each script or function (where noted). If you are not using a windows computer, you must also manually change these paths.

## Data analysis and plotting
The following scripts and functions are used for creating part or all of a respective figure from Williams et.al. 2024.
* Fig2_SNR.m
* Fig3_Fig6_2FP.m
  * load2ColorData.m
* Fig4_MM_Mouse.m
  * rasterplot.m
  * rasterplot_2colors.m
* Fig5_movAvg_PCC.m
  * Count_CTCs_per_interval.m

