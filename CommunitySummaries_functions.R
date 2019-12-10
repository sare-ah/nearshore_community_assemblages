#######################################################################################################################
# Benthic Habitat Mapping Dive Survey Data Exploration
# 
# Objective:  Functions necessary for CommunitySummaries.R script
# 
# Author:     Sarah Davies
#             Sarah.Davies@dfo-mpo.gc.ca
#             250-756-7124
# Date:       October 9, 2019
######################################################################################################################


# # Determine frequency of occurence for each species
# spp.freq.Original <- function(speciesObs, speciesLut, samplingUnits){
#   # Build a dataframe with species counts
#   sppNmes <- colnames(spp)
#   cnts <- colSums(spp)
#   sppCnts <- data.frame(sppNmes, cnts)
#   # Count the number of occurences for each species
#   colnames(sppCnts) <- c("Species_Code","Counts")
#   sppCnts$Species_Code <- as.character(sppCnts$Species_Code)
#   #head(sppCnts)
#   # Join each foo table to look-up table to add Hart codes
#   sppSciNme <- full_join( sppCnts, luTbl, by="Species_Code" )
#   # Remove extra fields 
#   sppSciNme <- dplyr::select( sppSciNme, Species_Code, Name, Counts )
#   # Determine the frequency of occurence throughout all dive sites
#   sppSciNme$Counts[is.na(sppSciNme$Counts)] <- 0   # Set any NA values to 0
#   sppSciNme$Freq <- round( sppSciNme$Counts/nUnits, digits=3 )  # Number of sampling units = nUnits
#   # Arrange in correct order 
#   #df <- dplyr::select( sppSciNme, Name, Species_Code, Counts,Freq)
#   #df <- dplyr::select( Species_Code, Counts, Freq)
#   df <- df[order(-df$Freq),]
#   df[is.na(df)] <- " "
#   return(df)
# }

# Determine frequency of occurence for each species
spp.freq <- function(speciesObs){
  # Build a dataframe with species counts
  sppNmes <- colnames(spp)
  cnts <- colSums(spp)
  sppCnts <- data.frame(sppNmes, cnts)
  # Count the number of occurences for each species
  colnames(sppCnts) <- c("Species_Code","Counts")
  sppCnts$Species_Code <- as.character(sppCnts$Species_Code)
  # Determine the frequency of occurence throughout all dive sites
  sppCnts$Counts[is.na(sppCnts$Counts)] <- 0   # Set any NA values to 0
  sppCnts$Freq <- round( sppCnts$Counts/nrow(spp), digits=3 )  # Number of sampling units = nUnits
  # Arrange in correct order 
  sppCnts <- dplyr::select(sppCnts, Species_Code, Counts, Freq)
  sppCnts <- sppCnts[order(-sppCnts$Freq),]
  sppCnts[is.na(sppCnts)] <- " "
  # Remove irriating region prefix
  sppCnts$Species_Code = gsub("^.*\\.","", sppCnts$Species_Code)
  return(sppCnts)
}



# 
# my.list <- list(name = c("A", "B", "C"), age = c(20, 30, 40, 50))
# lengths(my.list)
# 
# dat <- function(df, arrange_var) {
#   quo_arrange_var <- enquo(arrange_var)
#   depths %>%
#     melt(id="Name",value.name="Depths") %>%
#     arrange(-((!!quo_arrange_var)==2), !!quo_arrange_var)
# }
# 
# dat <- dat(dcat, Species_Code)
# 
# arrange_var <- dcat$Species_Code
# quo_arrange_var <- enquo(arrange_var)
# depths %>%
#   melt(id="Name",value.name="Depths") %>%
#   arrange(-((!!quo_arrange_var)==2), !!quo_arrange_var)
# 
# 
#  
# depthsA <- dcat %>%
#   group_by(Species_Code) %>%
#   mutate(minD=min(DepthCat)) %>%
#   mutate(maxD=max(DepthCat)) %>%
#   select(Species_Code,minD,maxD) %>%
#   distinct(Species_Code,minD,maxD)
# 
# 
# # 
# #3B.Plot
# # dat <- function(df, arrange_var) {
# #   quo_arrange_var <- enquo(arrange_var)
# #   depths %>%
# #   melt(id="Name",value.name="Depth") %>%
# #   arrange(-((!!quo_arrange_var)==2), !!quo_arrange_var)
# # }
# 
# 
