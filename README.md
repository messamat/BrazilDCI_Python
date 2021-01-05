# Safeguarding migratory fish via strategic planning of future small hydropower in Brazil - Python scripts
Thiago B. A. Couto, Mathis L. Messager & Julian D. Olden

Purpose: We quantified the trade-offs between hydroelectric generation capacity and impacts on river connectivity for thousands of current and
projected-future dams across Brazil.

## Set-up
All GIS analyses in this study require an ESRI ArcGIS license including the Spatial and Network Analyst extensions.
We used the Python Arcpy module associated with ArcGIS 10.5 in Python 2.7 with 64-bit background processing.

These scripts are annotated but could be challenging to follow. Please contact Mathis L. Messager for comments and clarifications. 

Files needed to run this analysis are available by downloading the study's [figshare permanent repository](https://figshare.com/s/5ba67b7f58ccc812ae70).
The /data folder in the figshare repository contains raw data not currently available in the same form as was used in this study and 
the directory structure enables users to reproduce our study using the scripts herein.

Additional sources of data to download prior to run this analysis include (folder to insert them in: link to data source):
- */data/ibge/BR/BRUFE250GC_SIR.shp : from the ibge ftp: ftp://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2018/Brasil/BR 
- */data/HydroATLAS/RiverATLAS_v10.gdb/RiverATLAS_v10: from [HydroSHEDS website](https://hydrosheds.org/page/hydroatlas) (Linke et al. 2019)  
- */data/HydroATLAS/BasinATLAS_v10.gdb/BasinATLAS_lev0{4-8}_v10: from [HydroSHEDS website](https://hydrosheds.org/page/hydroatlas) (Linke et al. 2019)  

## Workflow
*/src/python/dci_format.py: 
	Contains the python script used in spatially pre-processing (GIS analysis) the river network and hydroelectric dam datasets
	Requires downloading the data sources cited above (keeping the same directory structure).
  
 ## Dendritic Connectivity Index analysis
 Subsequent data analysis and figure production was performed in R 4.0. The corresponding source codes are available at [messamat/BrazilDCI_R](https://github.com/messamat/BrazilDCI_R)
