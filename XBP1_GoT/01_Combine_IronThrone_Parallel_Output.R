library(dplyr)
library(tidyr)
library(parallel)

#Read in options from Rscript command to set working directory, pcr ratio threshold, levenshtein distance for UMI collapse, and pcr duplicate threshold ####
options <- commandArgs(trailingOnly = TRUE)
wd <- options[1]
pcr_thresh <- 0.5 # PCR ratio threshold for assigning WT or MUT after collapsing UMIs. 
ld <- 1 #Levenshetein distance for declaring two UMIs "similar enough to merge"
dupcut <- 2 #minimum number of WT+MUT reads needed for keeping a UMI

setwd(wd)

#Run for individual samples
HP23166_0h_S1 <- read.table("HP23166_0h_S1.summTable.txt",header=T,sep="\t")


#Collapse UMIs #####
raw_GoT_table <- HP23166_0h_S1
pcr_ratio_thresh <- pcr_thresh


#Function to collapse UMIs for a single cell barcode
list_collapse <- function(single_got_row){
  UMIs <- unlist(strsplit(single_got_row[,"UMI"], ";")) #Separate all UMIs
  num.WT.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.WT.in.dups"], ";"))) #Separate number of WT PCR duplicates for each UMI
  num.MUT.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.MUT.in.dups"], ";"))) #Separate number of MUT PCR duplicates for each UMI
  num.amb.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.amb.in.dups"], ";"))) #Separate number of ambiguous PCR duplicates for each UMI
  call.in.dups <- unlist(strsplit(single_got_row[,"call.in.dups"], ";")) #Separate WT/MUT/Amb call for each UMI
  
  
  match_list <- lapply(UMIs, FUN =  function(x) agrep(x, UMIs, max.distance = ld)) #For each UMI, use the agrep function to identify all similar UMIs within the Levenshtein distance threshold
  num_of_matches <- unlist(lapply(match_list, FUN = function(x) length(x))) #Count the number of similar UMIs for each UMI
  
  #Continue to collapse while there are similar UMIs remaining in the set
  while (sum(num_of_matches) > length(num_of_matches)){
    to_collapse <- which(num_of_matches == max(num_of_matches))[1] #Identify the UMI with the most number of similar UMIs in the set
    matches_t0 <- numeric()
    matches_t1 <- match_list[[to_collapse]] #Identify the group of UMIs to be collapsed
    while (length(matches_t1) > length(matches_t0)){
      to_add <- unique(unlist(match_list[c(matches_t1)]))
      matches_t0 <- matches_t1
      matches_t1 <- to_add
    }
    
    WT_dups <- num.WT.in.dups[matches_t1] #Identify the number of WT PCR duplicates for each UMI to be collapsed
    MUT_dups <- num.MUT.in.dups[matches_t1] #Identify the number of MUT PCR duplicates for each UMI to be collapsed
    AMB_dups <- num.amb.in.dups[matches_t1] #Identify the number of Amb PCR duplicates for each UMI to be collapsed
    num.WT.in.dups[to_collapse] <- sum(WT_dups) #Add the WT PCR duplicates for all UMIs to be collapsed
    num.MUT.in.dups[to_collapse] <- sum(MUT_dups) #Add the MUT PCR duplicates for all UMIs to be collapsed
    num.amb.in.dups[to_collapse] <- sum(AMB_dups) #Add the Amb PCR duplicates for all UMIs to be collapsed
    if (sum(WT_dups) + sum(MUT_dups) == 0){
      call.in.dups[to_collapse] <- "AMB" #If no WT or MUT supporting PCR duplicates, the UMI is ambiguous
    } else {
      pcr_ratio <- (max(c(sum(WT_dups), sum(MUT_dups)))[1])/(sum(WT_dups) + sum(MUT_dups)) #Calculate the ratio of supporting PCR duplicates in favor of the majority genotype (WT or MUT)
      if(pcr_ratio > pcr_ratio_thresh){
        call.in.dups[to_collapse] <- ifelse(sum(WT_dups) > sum(MUT_dups), "WT", "MUT") #If the ratio passes the specified threshold, call the UMI in favor of the majority genotype
      } else {
        call.in.dups[to_collapse] <- "AMB" #If the ratio is below the specified threshold, the UMI is ambiguous
      }
    }
    
    #Remove the UMIs and corresponding data for all but the UMI to which the group was collapsed
    matches_rm <- matches_t1[matches_t1 != to_collapse] 
    UMIs <- UMIs[-matches_rm]
    num.WT.in.dups <- num.WT.in.dups[-matches_rm]
    num.MUT.in.dups <- num.MUT.in.dups[-matches_rm]
    num.amb.in.dups <- num.amb.in.dups[-matches_rm]
    call.in.dups <- call.in.dups[-matches_rm]
    
    #Recalculate the number of matches within the UMI set for the single cell barcode
    match_list <- lapply(UMIs, FUN =  function(x) agrep(x, UMIs, max.distance = ld))
    num_of_matches <- unlist(lapply(match_list, FUN = function(x) length(x)))
  }
  #Reformat the data for remaining UMIs following collapsing back into the format of the original IronThone data frame
  single_got_row[,"UMI"] <- paste0(UMIs, collapse = ";")
  single_got_row[,"num.WT.in.dups"] <- paste0(num.WT.in.dups, collapse = ";")
  single_got_row[,"num.MUT.in.dups"] <- paste0(num.MUT.in.dups, collapse = ";")
  single_got_row[,"num.amb.in.dups"] <- paste0(num.amb.in.dups, collapse = ";")
  single_got_row[,"call.in.dups"] <- paste0(call.in.dups, collapse = ";")
  single_got_row[,"WT.calls"] <- sum(call.in.dups == "WT")
  single_got_row[,"MUT.calls"] <- sum(call.in.dups == "MUT")
  single_got_row[,"amb.calls"] <- sum(call.in.dups == "AMB")
  
  #Keep only those UMIs with a total number of supporting PCR duplicates above the pre-specified threshold
  dup_thresh <- dupcut
  
  UMIs <- unlist(strsplit(single_got_row[,"UMI"], ";"))
  num.WT.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.WT.in.dups"], ";")))
  num.MUT.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.MUT.in.dups"], ";")))
  num.amb.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.amb.in.dups"], ";")))
  call.in.dups <- unlist(strsplit(single_got_row[,"call.in.dups"], ";"))
  sum_reads <- num.WT.in.dups + num.MUT.in.dups
  threshold_filter <- sum_reads >= dup_thresh
  
  single_got_row[,"UMI"] <- paste0(UMIs[threshold_filter], collapse = ";")
  single_got_row[,"num.WT.in.dups"] <- paste0(num.WT.in.dups[threshold_filter], collapse = ";")
  single_got_row[,"num.MUT.in.dups"] <- paste0(num.MUT.in.dups[threshold_filter] , collapse = ";")
  single_got_row[,"num.amb.in.dups"] <- paste0(num.amb.in.dups[threshold_filter], collapse = ";")
  single_got_row[,"call.in.dups"] <- paste0(call.in.dups[threshold_filter], collapse = ";")
  single_got_row[,"WT.calls"] <- sum(call.in.dups[threshold_filter] == "WT")
  single_got_row[,"MUT.calls"] <- sum(call.in.dups[threshold_filter] == "MUT")
  single_got_row[,"amb.calls"] <- sum(call.in.dups[threshold_filter] == "AMB")
  
  return(single_got_row)
}

#Run the single-barcode UMI collapsing function on the entire dataset
UMI_collapse <- function(GoT_table){
  GoT_table_to_collapse <- data.frame(BC = GoT_table[,"BC"],
                                      UMI = GoT_table[,"UMI"],
                                      num.WT.in.dups = GoT_table[,"num.WT.in.dups"],
                                      num.MUT.in.dups = GoT_table[,"num.MUT.in.dups"],
                                      num.amb.in.dups = GoT_table[,"num.amb.in.dups"],
                                      call.in.dups = GoT_table[,"call.in.dups"],
                                      WT.calls = "",
                                      MUT.calls = "",
                                      amb.calls = "",
                                      stringsAsFactors = FALSE)
  #Convert the GoT data frame into a list of single-row data frames where each entry is the data for a single cell barcode
  GoT_list <- split(GoT_table_to_collapse, seq(nrow(GoT_table_to_collapse)))
  collapsed_GoT_list <- mclapply(GoT_list, FUN = list_collapse)
  #Convert collapsed list back to a data frame
  collapsed_GoT_table <- do.call("rbind", collapsed_GoT_list)
  return(collapsed_GoT_table)
}

max_collapsed_GoT_table_1 <- UMI_collapse(raw_GoT_table)

#Save output
write.table(max_collapsed_GoT_table_1, file = paste0(wd,"/collapsed/HP23166_0h_S1_myGoT.summTable.concat.umi_collapsed.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
            
            
            