# Set working directory
# setwd("/Volumes/ifs/DCEG/Branches/LTG/Chanock/Lea/HiCProConfig")
# setwd("/Volumes/ifs/DCEG/Branches/LTG/Chanock/Tim/HiC_Mai")

#######################################################################################################
# Purpose:  Take a region of interest (ROI) and expand it by finding the genomic coordinates
# upstream and downstream of that fragment that includes 4 restriction enzyme cutting sites.
# Send output to a date frame that includes 1, 2, 3, & 4 restriction enzyme sites.
# See schemaitc diagram below.
# -- NOTE:  On the downstream side of the ROI, the whole restriction fragment is not included in the
# -- output coordinates but only its 'start' coordinate.  However, on the upstream side of the ROI,
# -- the whole restriction fragment is included in the output.
#######################################################################################################

############################################################
# Schematic diagram of regions of interest (ROI) and
# restrction enzyme cut sites.
############################################################

# '|' defines boundaries of ROI's and restriction fragments
# '^' indicates restriction enzyme cut site

# |_^___|_^___|_^___|_^___|__________________________|_^___|_^___|_^___|_^___|
# |..4..|..3..|..2..|..1..|...........ROI............|..1..|..2..|..3..|..4..|
# |.......up stream.......|....region of interest....|......down stream......|


############################
# Two Data Files Required
############################
# Data File 1:  User defined list of genomic regions of interest (ROI) in 'bed' format (chr# start end)
#                                                                                      (e.g. chr6 1234567 3456789)
# Data File 2:  Genomic coordinate list of restriction enzyme cut fragments in 'bed' format, same as above.  (In this first case,
# the restriction fragments were from DpnII) and was generated using the "digest.genome.py" utility in the "hic-pro" tool
# (https://nservant.github.io/HiC-Pro/).  A copy of this file can be found at /Volumes/ifs/DCEG/Branches/LTG/Chanock/Lea/HiCProConfig.


########################
# Read in Data Files
########################
# Read Data File 1
roi_list <- read.delim("/Volumes/ifs/DCEG/Branches/LTG/Chanock/Tim/HiC_Mai/data/Mai_Regions_1.txt", header = FALSE, sep = "\t", as.is=T)
# Read Data File 2
enz_fragment <- read.delim("/Volumes/ifs/DCEG/Branches/LTG/Chanock/Lea/HiCProConfig/hg19_dpnII.bed", header = FALSE, sep = "\t", as.is=T)


#################################################################################################################################
# Begin FUNCTION
# Function for looping through Data File 2, the 'bed' file of DpnII genome-wide restriction fragments using input from region
# of interest (ROI).  Input into this function is chromosome number, genomic start coordinate, genomic stop coordinate
# for the ROI.
#################################################################################################################################
roi_expand <- function(chr_roi, start_roi, end_roi) {
  # begin for1 loop, upstream of ROI
  for(i in 1:nrow(enz_fragment)) {
    # 1st 'if' statement to evaluate chromosome number
    if(chr_roi == enz_fragment$V1[i]) {
      # 2nd 'if' statement to evaluate genomic start coordinate
      if(start_roi < enz_fragment$V2[i]) {
        x <- data.frame(chr = chr_roi, start.0 = start_roi, start.1 = enz_fragment$V2[i-1], start.2 = enz_fragment$V2[i-2], start.3 = enz_fragment$V2[i-3], start.4 = enz_fragment$V2[i-4], stringsAsFactors = FALSE)
        break
      # end if2
      }
  # end if1
    }
  # end for1 loop
  }

  # begin for2 loop, downstream of ROI
  for(j in i:nrow(enz_fragment)) {
    if(end_roi < enz_fragment$V3[j]) {
      y <- data.frame(chr.a = chr_roi, end.0 = end_roi, end.1 = enz_fragment$V2[j+1], end.2 = enz_fragment$V2[j+2], end.3 = enz_fragment$V2[j+3], end.4 = enz_fragment$V2[j+4], stringsAsFactors = FALSE)
      break
    # end if1
    }

  # end for2 loop
  }
  # combine the two data frames; 'x' is from upstream of ROI; 'y' is from downstream of ROI
  z <- cbind(x,y)
  return(z)
}
########################
# End function
########################


#######################################################################################
# Batch call; loop through Data File 1 (roi_list) and call the function (roi_expand)
#######################################################################################

# create empty data frame
df_expand <- data.frame()
# 'for' loop through ROI list; Data File 1
  for(k in 1:nrow(roi_list)) {
    chr_roi <- roi_list$V1[k]
   start_roi <- roi_list$V2[k]
    end_roi <- roi_list$V3[k]

    # Function call
   df_tmp <- roi_expand(chr_roi, start_roi, end_roi)
    # Append to the bottom of the data frame 'df_expand'
    df_expand <- rbind(df_expand, df_tmp)
  }

# Rearrange data frame returned from function 'roi_expand'
# delete column 7, "chr.a", second listing of the chromosome number, not needed anymore
df_expand[7] <- NULL
# reorder data frame by column index
df_expand <- df_expand[,c(1, 2, 7, 3, 8, 4, 9, 5, 10, 6, 11)]
head(df_expand)

# write to csv file
write.csv(df_expand, file="/Volumes/ifs/DCEG/Branches/LTG/Chanock/Tim/HiC_Mai/output/hg19_dpnII_final_list.csv")


#################################################################################################
# Single Function Call which uses only the first line of roi_list (primarily used for testing)
#################################################################################################
# Initial input variables
chr_roi <- roi_list$V1[1]
start_roi <- roi_list$V2[1]
end_roi <- roi_list$V3[1]

# Single Function call
df_expand <- roi_expand(chr_roi, start_roi, end_roi)
# Rearrange data frame
# delete column 7, "chr.a", second listing of the chromosome number, not needed anymore
df_expand[7] <- NULL
# reorder data frame by column index
df_expand <- df_expand[,c(1, 2, 7, 3, 8, 4, 9, 5, 10, 6, 11)]

# write to csv file
# write.csv(df_expand, file="/Volumes/ifs/DCEG/Branches/LTG/Chanock/Tim/HiC_Mai/output/dpnII_test7.csv")


#####################################
# Record Session Info
#####################################
writeLines(capture.output(sessionInfo()), "/Volumes/ifs/DCEG/Branches/LTG/Chanock/Tim/HiC_Mai/sessionInfo/HiC_Mai_regions_v0.1.5_sessionInfo.txt")


