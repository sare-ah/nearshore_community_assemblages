#########################################################################
# Prep data for community assemblages analysis from BHM dive survey data
# 
# Objective: Organize data to for community assemblages analysis
#
# Notes:      This script assumes that observations have already been 
#             extracted from BHM db (using 32-bit R) using 
#             '1.ExtractDataFromMSAccess.R' to create SpeciesObs.csv
#
# Author:     Sarah Davies
#             Sarah.Davies@dfo-mpo.gc.ca
#             250-756-7124
# Date:       October 1st, 2019
#########################################################################
 
# start fresh
rm(list=ls())

# Check which version of R is being used and reset if necessary
Sys.getenv("R_ARCH")   

# 1. Update Observations
# A. Recode species observations to fix typos & add Invert/Algae designation to species code
# Input = SpeciesObs.csv; fields = HKey,Quadrat,SpNum,SpType,Species
# Output = SpeciesObs_updated.csv; fields = HKey,Quadrat,SpNum,SpType,Species,Species_Code
source('C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Scripts/2.CreateUpdatedSppObs.R')

# B1. Build species by site matrix
# Input = SpeciesObs_updated.csv, Quadrat.csv, Headers.csv
# Output = AllSpeciesByDepthCategory.csv; fields = TransDepth,A_AA,...,I_ZS *** Need to add HKey here? ***
source('C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Scripts/3c.A&I_BuildMtrx.R')

# C. Recode the substrate category from 11 to 4
# Input = Quadrat.csv; fields = HKey,Quadrat,GaugeDepth,CorDepthM,DepthCat,Time,Substrate1,...,DriftSpecies
# Output = Quadrat_RMSM.csv; fields = 	HKey,Quadrat,GaugeDepth,...,Sub.cat,SubstrateCat.Nme,BType1,BType2,BType.dscrptn
source('C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Scripts/ConvertSubstrate2RMSM.R')

# D. Summarise dominant substrates and calculate mean depth and change in elevation for each depth category
# Input = Quadrat.csv OR Quadrat_RMSM.csv; fields = HKey,Quadrat,GaugeDepth,...,BType.dscrptn
# Output = DepthCat_CalcSub.csv OR DepthCat_CalcSub_RMSM.csv; fields = HKey,Quadrat,MeanDepth.m,Slope,...,BType.dscrptn
source('C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Scripts/SummariseSubstrate.R')

# 3. Spatialize points *** Code partially complete ***
# Input = Quadrat.csv, Headers.csv
# Output = SpatialPts.shp; fields = HKey,LatStart,LatEnd,LonStart,LonEnd
source('C:/Users/daviessa/Documents/R/PROJECTS_OTHERS/SpatializeDiveTransects/Prep4SpatializeDiveTransects.R')

# 4a. Export bathy derivatives for observations (ExportRastersAsTiffs.py - Python script)
# and attach to species observations
# source('F:/R/MY_PROJECTS/RCPmod/Scripts/1.ExtractAttachEnviro.R')

# Extract environmental variables

# 4b. Add fetch values to each header key (transect)
# Input = Transect.coords.csv; fields = HKey,LatStart,LatEnd,LonStart,LonEnd
# Output = Trans.coords.fetch.csv; fields = HKey,LatStart,LatEnd,LonStart,LonEnd,Sum_Fetch
source('C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Scripts/AddFetch2Transect.R')

# 5. Remove correlated environmental variables
# source('F:/R/MY_PROJECTS/RCPmod/Scripts/2.RemoveCorrelatedEnviros.R')

# # B2. Remove rare species and sites 
# # Input = AllSpeciesByDepthCategory.csv
# # Output = remRareByDepthCat.csv
# source('C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Scripts/RemoveRareSpp.R')

# 6. Combining all data together into final matrices
# Input = remRareByDepthCategory.csv,DepthCat_CalcSub.csv,DepthCat_CalcSub_RMSM.csv,Trans.coords.fetch.csv, SpatialPts.shp
# Output = ../CommunityAssemblages/Data/byDepthCat.csv
source('C:/Users/daviessa/Documents/R/PROJECTS_MY/DiveSurveys_DataPrep/Scripts/BuildSppEnviroMtrx.R')




# Data Exploration
# Figures and summaries of species data
# Do I need another script to summarise environmental data?
#source()

# Cluster Tendency
# Various tests for clustering tendency
#source('ClusterTendency.R')


# Cluster analysis
#source('F:/R/MY_PROJECTS/Cluster_Analysis/Scripts/clusteranalysis/1.ClusterAnalysis.R')

# ################################
# # Import rasters for analysis  #
# ################################
# setwd("T:/Environmental_Layers/Rasters for SDMs/Nearshore/")
# ls <- c("HG/Rasters/bathy.tif", "NCC/Rasters/bathy.tif", "QCS/Rasters/bathy.tif", "SoG/Rasters/bathy.tif", "WCVI/Rasters/bathy.tif")
# 
# rasfiles<-list.files(getwd(), pattern = "(*.)tif$",recursive=T)
# #ras<-stack(rasfiles[c(1:3, 5:30)]) #Skip Btype
# ras <- raster::stack(rasfiles) 
# 
# ##################################################
# # Bring in spatialized points with species data 
# ##################################################
# setwd("C:/Users/daviessa/Documents/R/PROJECTS_OTHERS/SpatializeDiveTransects/Haida_Gwaii/Output")
# sp<-readOGR("HG_bhm_pts.shp")
# sp<-spTransform(sp, CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs"))
# proj4string(ras)==proj4string(sp)
# 
# ########################################################
# # Extract enviromental data from rasters to points
# ########################################################
# spEnv = data.frame(sp,raster::extract(ras, sp)) ##add .parameters to each record
# summary(is.na(spEnv)) 
# 
# spEnv_complete<-spEnv[complete.cases(spEnv[!names(spEnv) %in% names(sp)]),] #remove cases where environmental data aren't available
# summary(is.na(spEnv_complete)) 
# write.csv(spEnv_complete,"C:/Users/daviessa/Documents/R/PROJECTS_OTHERS/Cluster_Analysis/Models/2019.04_HG/SpeciesByDepthSpatialEnviro.csv", row.names=F)
# head(spEnv_complete,3)
 
