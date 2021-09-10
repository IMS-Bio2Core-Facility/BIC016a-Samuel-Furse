rm(list = ls())
options(warn=-1)

#### PRELIMINARIES ############################################################################################# 

#*Uploads the needed libraries --------------------------------------------------------------------------------

require(ggplot2)

require(data.table)

require(plotly)

require(DT)

require(R2HTML)

require(stringr)

#** Sets the number of significant digits for the output --------------------------
sig_dig = 4

#** Sets the working directory ---------------------------------------------------------------------------------

#Gets the default wd
default_wd <- getwd()

setwd("C:/Users/Furse/.../300. Data outputs/")# <--- insert here the path to the working (output) directory

new_wd <- getwd()

#Sets the input directory
inputdir <-"C:/Users/Furse/.../200. Data sheets/"# <--- insert here the path to the input directory


#### DEFINES FUNCTIONS ############################################################################################# 

count_zeroes <- function(x){length(which(x==0))}
#considered_mode <- "\\+ve"
considered_mode <- "\\+ve"
considered_generation <- "PWD" 
considered_model_1 <-"lean"
considered_model_2 <- "obese"
percentage_of_zeroes <- 33 # Sets a threshold to exclude rows from the analysis: the rows containing more than [percentage_of_zeroes/100] zeroes will be excluded from the analysis
#### DATA UPLOAD ###############################################################################################################
# - Uploads the .csv files containing the information related to the F1A, PW, -ve datasets.
# - Separates the metadata from the main data
# - Picks only the considered_model_1 and considered_model_2 sets

files_names_originals <- list.files(inputdir)
files_names_originals <- files_names_originals[which(str_length(files_names_originals)==17)] 

modes <- unique(substr(files_names_originals, start=1, stop=3))
tissues <- unique(substr(files_names_originals, start=6, stop=8))
generations <- unique(substr(files_names_originals, start=11, stop=13)) 

# PRODUCES TISSUE-SPECIFIC MATRICES
# For each tissue, produces two matrices (one for each considered_model) 
# in which the columns are the samples and the rows are the lipids 

for(j in 1: length(tissues)){
  
  tissue <- tissues[j]
  aa <- files_names_originals[grep(files_names_originals, pattern=considered_mode)]
  aa <- aa[grep(aa, pattern=considered_generation)]
  aa<- aa[grep(aa, pattern=tissue)]
  
  if(length(aa)!=0){
    
    bb<- read.csv(paste0(inputdir, aa), stringsAsFactors = F)
    assign(paste0(tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_metadata"), bb)
    
    cc<- read.csv(paste0(inputdir, aa), stringsAsFactors = F , skip=10 ) 
    cc <- cc[!is.na(cc$m.z),]
    assign(paste0(tissue, "_", considered_generation, "_", substr(considered_mode, 2,4)), cc)
    
    cc_1 <- cc[,grep(as.vector(bb[4, ]),pattern=considered_model_1)]
    cc_2 <- cc[,grep(as.vector(bb[4, ]),pattern=considered_model_2)]
    rownames(cc_1) <- cc$Lipid.variable
    rownames(cc_2) <- cc$Lipid.variable
    
    cc_1_zeroes <- apply(cc_1[, c(3:ncol(cc_1))], 1, count_zeroes)
    cc_2_zeroes <- apply(cc_2[, c(3:ncol(cc_2))], 1, count_zeroes)
    
    cc_1_nozeroes <- cc_1[-which(cc_1_zeroes > ncol(cc_1)*percentage_of_zeroes/100) ,] # Excludes the rows in which the number of zeroes is greater than 20% 
    cc_2_nozeroes <- cc_2[-which(cc_2_zeroes > ncol(cc_2)*percentage_of_zeroes/100) ,]
    
    assign(paste0(tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1), cc_1)
    assign(paste0(tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_2), cc_2)
    
    assign(paste0(tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_nozeroes"), cc_1_nozeroes)
    assign(paste0(tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_2, "_nozeroes"), cc_2_nozeroes)
  }
}

#### A LIPIDS ###############################################################################################################

# Finds the A-type lipids for both the considered models -------------------------------------------------------------------------------------------------- 

# ** Considered_model_1 --------------------------------------------------------------------------------------

# Finds the A-lipids

# Creates a list (called all_tissues) in which each element is 
# a list of lipids coming from the row names of
# each Tissue-specific Matrix
all_tissues <-list()
for(k in 1:length(tissues)){
  tissue <- tissues[k]
  aa <- files_names_originals[grep(files_names_originals, pattern=considered_mode)]
  aa <- aa[grep(aa, pattern=considered_generation)]
  aa<- aa[grep(aa, pattern=tissue)]
  
  if(length(aa)!=0){
    print(tissue)
    yy <- rownames(get(paste0(tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_nozeroes") ) ) 
    assign(paste0("xx_",k), yy)
    all_tissues[[k]] <- yy
    
  }
  
}
names(all_tissues) <- tissues

if(length(which(all_tissues=="NULL")) !=0){
  all_tissues<-all_tissues[-which(all_tissues=="NULL")]
}

# intersects all the elements of the all_tissues list to find the A-lipids  
gg_1<-Reduce(intersect, all_tissues)
assign(paste0 ("A_Lipids_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1), gg_1)
write.csv(gg_1, file=paste0 ("A_Lipids_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, ".csv"))

# Counts the A-lipids in each classes
classes <- unique(substr(gg_1, start=1, stop =5))
classes_counts_model_1 <- matrix(ncol=1, nrow=length(classes))
rownames(classes_counts_model_1) <- classes 
colnames(classes_counts_model_1) <- paste0("A_lipids_",considered_model_1)
for(i in 1:length(classes)){
  
  classes_counts_model_1[i] <- length(grep(gg_1, pattern=paste0("^",classes[i])))
}

# ** Considered_model_2 --------------------------------------------------------------------------------------

# Finds the A-lipids

# Creates a list (called all_tissues) in which each element is 
# a list of lipids coming from the row names of
# each Tissue-specific Matrix
all_tissues <-list()
for(k in 1:length(tissues)){
  tissue <- tissues[k]
  aa <- files_names_originals[grep(files_names_originals, pattern=considered_mode)]
  aa <- aa[grep(aa, pattern=considered_generation)]
  aa<- aa[grep(aa, pattern=tissue)]
  
  if(length(aa)!=0){
    print(tissue)
    yy <- rownames(get(paste0(tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_2, "_nozeroes") ) ) 
    assign(paste0("xx_",k), yy)
    all_tissues[[k]] <- yy
    
  }
  
}
names(all_tissues) <- tissues
if(length(which(all_tissues=="NULL")) !=0){
  all_tissues<-all_tissues[-which(all_tissues=="NULL")]
}

# intersects all the elements of the all_tissues list to find the A-lipids  
gg_2<-Reduce(intersect, all_tissues)
assign(paste0 ("A_Lipids_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_2), gg_2)
write.csv(gg_2, file=paste0 ("A_Lipids_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_2, ".csv"))


# Counts the A-lipids in each class
classes <- unique(substr(gg_2, start=1, stop =5))
classes_counts_model_2 <- matrix(ncol=1, nrow=length(classes))
rownames(classes_counts_model_2) <- classes 
colnames(classes_counts_model_2) <- paste0("A_lipids_",considered_model_2)
for(i in 1:length(classes)){
  
  classes_counts_model_2[i] <- length(grep(gg_2, pattern= paste0("^",classes[i])))
} 

# Creates the A_lipids_classes_counts_tot dataframe, where the vectors classes_counts_model_1  and classes_counts_model_2 are merged
classes_counts_model_1 <- as.data.frame(classes_counts_model_1)
classes_counts_model_2 <- as.data.frame(classes_counts_model_2)
A_lipids_classes_counts_tot <- merge(classes_counts_model_1, classes_counts_model_2, by.x="row.names", by.y="row.names", all=T)
A_lipids_classes_counts_tot[is.na(A_lipids_classes_counts_tot)]<-0
rownames(A_lipids_classes_counts_tot) <- A_lipids_classes_counts_tot[,1]
A_lipids_classes_counts_tot <- A_lipids_classes_counts_tot[,-1]
#
All_glicerids <- A_lipids_classes_counts_tot[grep(rownames(A_lipids_classes_counts_tot), pattern="DGX|MGX|TGX"),]
A_lipids_classes_counts_tot[nrow(A_lipids_classes_counts_tot)+1,] <- colSums(All_glicerids)
rownames(A_lipids_classes_counts_tot)[nrow(A_lipids_classes_counts_tot)] <- "Glyc"

assign(paste0("A_lipids_classes_counts_", considered_generation, "_",  substr(considered_mode, 2,4), "_", considered_model_1, "_", considered_model_2), A_lipids_classes_counts_tot)
write.csv(A_lipids_classes_counts_tot, file=paste0("A_lipids_classes_counts_", considered_generation, "_",  substr(considered_mode, 2,4), "_", considered_model_1, "_", considered_model_2,".csv"))

# Computes the Jaccard distances between models -------------------------------------------------------------------------------------------------- 

# This part is not needed for the computation Jaccard distances but only for printing the  A_lipids_all_matrix - START
# Merges the vectors containing the A lipids "substituted" for each model into the A_lipids_matrix 
ss <- union(gg_1, gg_2)
ss<- sort(ss)
#
A_lipids_matrix<- cbind(ss, rep(0, length(ss)), rep(0, length(ss)) )
A_lipids_matrix[,2][which(A_lipids_matrix[,1] %in% gg_1)] <- A_lipids_matrix[,1][which(A_lipids_matrix[,1] %in% gg_1)]
A_lipids_matrix[,3][which(A_lipids_matrix[,1] %in% gg_2)] <- A_lipids_matrix[,1][which(A_lipids_matrix[,1] %in% gg_2)]
#
colnames(A_lipids_matrix) <- c("rownames",considered_model_1, considered_model_2)
rownames(A_lipids_matrix) <- A_lipids_matrix[,1]
A_lipids_matrix <- A_lipids_matrix[,-1]
assign(paste0 ("A_Lipids_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2), A_lipids_matrix)
write.csv(file= paste0 ("A_Lipids_all_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2,".csv"), A_lipids_matrix)
# This part is not needed for the computation Jaccard distances but only for printing the  A_lipids_all_matrix - END


# Unifies the classes MG, DG and TG under the class Glyc, by changing the row names of the A-lipids 
gg_1_substituted <- gsub("MGXX|DGXX|TGXX", "Glyc", gg_1)
gg_2_substituted <- gsub("MGXX|DGXX|TGXX", "Glyc", gg_2)

# Merges the vectors containing the A lipids "substituted" for each model into the A_lipids_matrix 
ss <- union(gg_1_substituted, gg_2_substituted)
ss<- sort(ss)
#
A_lipids_matrix<- cbind(ss, rep(0, length(ss)), rep(0, length(ss)) )
A_lipids_matrix[,2][which(A_lipids_matrix[,1] %in% gg_1_substituted)] <- A_lipids_matrix[,1][which(A_lipids_matrix[,1] %in% gg_1_substituted)]
A_lipids_matrix[,3][which(A_lipids_matrix[,1] %in% gg_2_substituted)] <- A_lipids_matrix[,1][which(A_lipids_matrix[,1] %in% gg_2_substituted)]
#
colnames(A_lipids_matrix) <- c("rownames",considered_model_1, considered_model_2)
rownames(A_lipids_matrix) <- A_lipids_matrix[,1]
A_lipids_matrix <- A_lipids_matrix[,-1]
assign(paste0 ("A_Lipids_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2), A_lipids_matrix)
write.csv(file= paste0 ("A_Lipids_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2,".csv"), A_lipids_matrix)

# Computes the Jaccard distances
require(jaccard)

if(considered_mode =="\\+ve"){
  classes <- rownames(A_lipids_classes_counts_tot)[-grep(rownames(A_lipids_classes_counts_tot), pattern="DG|TGX|MGX")]
  classes <- c(classes, "Glyc")
  Jaccard_distances <- matrix(ncol=2, nrow = length(classes))
  colnames(Jaccard_distances) <- c("Distance", "Pvalue")
  rownames(Jaccard_distances) <- classes
} else {
  classes <- rownames(A_lipids_classes_counts_tot)
  Jaccard_distances <- matrix(ncol=2, nrow = length(classes))
  colnames(Jaccard_distances) <- c("Distance", "Pvalue")
  rownames(Jaccard_distances) <- classes
}

Global_jaccard_matrix <- A_lipids_matrix
Global_jaccard_matrix[which(Global_jaccard_matrix!=0)] = 1

for(i in 1:length(classes)){  
  zz<- as.matrix(Global_jaccard_matrix[grep(rownames(Global_jaccard_matrix), pattern=classes[i]),])
  if(length(grep(rownames(Global_jaccard_matrix), pattern=classes[i]))==1 ){
    uu<- jaccard(as.numeric(zz[1]), as.numeric(zz[2]))
    vv <-jaccard.test(as.numeric(zz[1]), as.numeric(zz[2]), method = "exact")
    Jaccard_distances[i,c(1,2)] <- c(uu,vv$pvalue)
  }
  else if (length(grep(rownames(Global_jaccard_matrix), pattern=classes[i]))==0) {
    Jaccard_distances[i,c(1,2)] <- c("NA","NA")
  }
  else{
    uu<- jaccard(as.numeric(zz[,1]), as.numeric(zz[,2]))
    vv <-jaccard.test(as.numeric(zz[,1]), as.numeric(zz[,2]), method = "exact")
    Jaccard_distances[i,c(1,2)] <- c(uu,vv$pvalue)
  }
  #print(zz)
}

assign(paste0 ("Jaccard_distances_A_Lipids_", "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2), Jaccard_distances)
write.csv(file=paste0 ("Jaccard_distances_A_Lipids_", "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2,".csv"), Jaccard_distances)

#### U LIPIDS ###############################################################################################################

# Finds the U-lipids for each tissue

# ** Considered_model_1 --------------------------------------------------------------------------------------
# Creates a list (called all_tissues_1) in which each element is 
# a list of lipids coming from the row names of
# each Tissue-specific and model-specific Matrix

all_tissues_1 <-list()
for(k in 1:length(tissues)){
  tissue <- tissues[k]
  aa <- files_names_originals[grep(files_names_originals, pattern=considered_mode)]
  aa <- aa[grep(aa, pattern=considered_generation)]
  aa <- aa[grep(aa, pattern=tissue)]
  
  if(length(aa)!=0){
    yy <- rownames(get(paste0(tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_nozeroes") ) ) 
    assign(paste0("xx_",k), yy)
    all_tissues_1[[k]] <- yy
    
  }
  
}
names(all_tissues_1) <- tissues
if(length(which(all_tissues_1=="NULL")) !=0){
  all_tissues<-all_tissues_1[-which(all_tissues_1=="NULL")]
  all_tissues_1<-all_tissues_1[-which(all_tissues_1=="NULL")]
}

# ** Considered_model_2 --------------------------------------------------------------------------------------
# Creates a list (called all_tissues_2) in which each element is 
# a list of lipids coming from the row names of
# each Tissue-specific and model-specific Matrix

all_tissues_2 <-list()
for(k in 1:length(tissues)){
  tissue <- tissues[k]
  aa <- files_names_originals[grep(files_names_originals, pattern=considered_mode)]
  aa <- aa[grep(aa, pattern=considered_generation)]
  aa <- aa[grep(aa, pattern=tissue)]
  
  if(length(aa)!=0){
    yy <- rownames(get(paste0(tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_2, "_nozeroes") ) ) 
    assign(paste0("xx_",k), yy)
    all_tissues_2[[k]] <- yy
    
  }
  
}
names(all_tissues_2) <- tissues
if(length(which(all_tissues_2=="NULL")) !=0){
  all_tissues<-all_tissues_2[-which(all_tissues_2=="NULL")]
  all_tissues_2<-all_tissues_2[-which(all_tissues_2=="NULL")]
}

# ** Tissue-specific U-lipids --------------------------------------------------------------------------------------

# **** Produces the Tissue-specific U-lipids matrices --------------------------------------------------------------------------------------

for(o in 1:length(names(all_tissues_1))){
  Tissue <- names(all_tissues)[o]
  
  considered_tissue_1 <- all_tissues_1[which(names(all_tissues_1)==Tissue)]
  other_tissues_1 <- all_tissues_1[-which(names(all_tissues_1)==Tissue)]
  
  ss_1 <- Reduce(union, other_tissues_1) # Lists the lipids present in all the tissues but the considered one
  tt_1 <- considered_tissue_1[[1]] # Lists the lipids present in the considered tissue
  uu_1<- setdiff(tt_1,ss_1) # Lists the lipids that are in the considered tissue but not in all the others
  
  considered_tissue_2 <- all_tissues_2[which(names(all_tissues_2)==Tissue)]
  other_tissues_2 <- all_tissues_2[-which(names(all_tissues_2)==Tissue)]
  
  ss_2 <- Reduce(union, other_tissues_2)
  tt_2 <- considered_tissue_2[[1]]
  uu_2<- setdiff(tt_2,ss_2)
  
  vv <- union(uu_1, uu_2)
  vv<- sort(vv)
  #
  U_lipids_matrix<- cbind(vv, rep(0, length(vv)), rep(0, length(vv)) )
  U_lipids_matrix[,2][which(U_lipids_matrix[,1] %in% uu_1)] <- U_lipids_matrix[,1][which(U_lipids_matrix[,1] %in% uu_1)]
  U_lipids_matrix[,3][which(U_lipids_matrix[,1] %in% uu_2)] <- U_lipids_matrix[,1][which(U_lipids_matrix[,1] %in% uu_2)]
  #
  colnames(U_lipids_matrix) <- c("rownames",considered_model_1, considered_model_2)
  rownames(U_lipids_matrix) <- U_lipids_matrix[,1]
  U_lipids_matrix <- U_lipids_matrix[,-1]
  assign(paste0 ("U_Lipids_", "all_", Tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2), U_lipids_matrix)
  write.csv(file= paste0 ("U_Lipids_", "all_", Tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2,".csv"), U_lipids_matrix)
  
  # **** Counts the U-lipids in each class --------------------------------------------------------------------------------------
  
  classes_1 <- unique(substr(uu_1, start=1, stop =5))
  classes_counts_model_1 <- matrix(ncol=1, nrow=length(classes_1))
  rownames(classes_counts_model_1) <- classes_1 
  colnames(classes_counts_model_1) <- paste0("U_lipids_",considered_model_1)
  
  if(length(uu_1) !=0){
    for(i in 1:length(classes_1)){
      
      classes_counts_model_1[i] <- length(grep(uu_1, pattern=classes_1[i]))
    } 
  }
  
  classes_2 <- unique(substr(uu_2, start=1, stop =5))
  classes_counts_model_2 <- matrix(ncol=1, nrow=length(classes_2))
  rownames(classes_counts_model_2) <- classes_2 
  colnames(classes_counts_model_2) <- paste0("U_lipids_",considered_model_2)
  
  if(length(uu_2) !=0){
    for(i in 1:length(classes_2)){
      
      classes_counts_model_2[i] <- length(grep(uu_2, pattern=classes_2[i]))
    } 
  }
  # Creates the U_lipids_classes_counts_tot dataframe, where the vectors classes_counts_model_1  and classes_counts_model_2 are merged
  classes_counts_model_1 <- as.data.frame(classes_counts_model_1)
  classes_counts_model_2 <- as.data.frame(classes_counts_model_2)
  U_lipids_classes_counts_tot <- merge(classes_counts_model_1, classes_counts_model_2, by.x="row.names", by.y="row.names", all=T)
  U_lipids_classes_counts_tot[is.na(U_lipids_classes_counts_tot)]<-0
  rownames(U_lipids_classes_counts_tot) <- U_lipids_classes_counts_tot[,1]
  U_lipids_classes_counts_tot <- U_lipids_classes_counts_tot[,-1]
  #
  All_glicerids <- U_lipids_classes_counts_tot[grep(rownames(U_lipids_classes_counts_tot), pattern="DGX|MGX|TGX"),]
  U_lipids_classes_counts_tot[nrow(U_lipids_classes_counts_tot)+1,] <- colSums(All_glicerids)
  rownames(U_lipids_classes_counts_tot)[nrow(U_lipids_classes_counts_tot)] <- "Glyc"
  
  assign(paste0("U_lipids_classes_counts_", Tissue, "_", considered_generation, "_",  substr(considered_mode, 2,4), "_", considered_model_1, "_", considered_model_2), U_lipids_classes_counts_tot)
  write.csv(U_lipids_classes_counts_tot, file=paste0("U_lipids_classes_counts_", Tissue, "_", considered_generation, "_",  substr(considered_mode, 2,4), "_", considered_model_1, "_", considered_model_2,".csv"))
  
  
  # Computes the Jaccard distances between models -------------------------------------------------------------------------------------------------- 
  
  # Unifies the classes MG, DG and TG under the class Glyc, by changing the row names of the A-lipids 
  uu_1_substituted <- gsub("MGXX|DGXX|TGXX", "Glyc", uu_1)
  uu_2_substituted <- gsub("MGXX|DGXX|TGXX", "Glyc", uu_2)
  
  # Merges the vectors containing the A lipids "substituted" for each model into the A_lipids_matrix 
  ss <- union(uu_1_substituted, uu_2_substituted)
  ss<- sort(ss)
  #
  U_lipids_matrix<- cbind(ss, rep(0, length(ss)), rep(0, length(ss)) )
  U_lipids_matrix[,2][which(U_lipids_matrix[,1] %in% uu_1_substituted)] <- U_lipids_matrix[,1][which(U_lipids_matrix[,1] %in% uu_1_substituted)]
  U_lipids_matrix[,3][which(U_lipids_matrix[,1] %in% uu_2_substituted)] <- U_lipids_matrix[,1][which(U_lipids_matrix[,1] %in% uu_2_substituted)]
  #
  colnames(U_lipids_matrix) <- c("rownames",considered_model_1, considered_model_2)
  rownames(U_lipids_matrix) <- U_lipids_matrix[,1]
  U_lipids_matrix <- U_lipids_matrix[,-1]
  assign(paste0 ("U_Lipids_", Tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2), U_lipids_matrix)
  write.csv(file= paste0 ("U_Lipids_", Tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2,".csv"), U_lipids_matrix)
  
  if(length(uu_2)!=0&length(uu_1)!=0){
    # Computes the Jaccard distances
    require(jaccard)
    
    classes <- rownames(U_lipids_classes_counts_tot)
    Jaccard_distances <- matrix(ncol=2, nrow = length(classes))
    colnames(Jaccard_distances) <- c("Distance", "Pvalue")
    rownames(Jaccard_distances) <- classes
    
    Global_jaccard_matrix <- U_lipids_matrix
    Global_jaccard_matrix[which(Global_jaccard_matrix!=0)] = 1
    
    for(i in 1:length(classes)){  
      zz<- as.matrix(Global_jaccard_matrix[grep(rownames(Global_jaccard_matrix), pattern=classes[i]),])
      if(length(grep(rownames(Global_jaccard_matrix), pattern=classes[i]))==1 ){
        uu<- jaccard(as.numeric(zz[1]), as.numeric(zz[2]))
        vv <-jaccard.test(as.numeric(zz[1]), as.numeric(zz[2]), method = "exact")
        Jaccard_distances[i,c(1,2)] <- c(uu,vv$pvalue)
      }
      else if (length(grep(rownames(Global_jaccard_matrix), pattern=classes[i]))==0) {
        Jaccard_distances[i,c(1,2)] <- c("NA","NA")
      }
      else{
        uu<- jaccard(as.numeric(zz[,1]), as.numeric(zz[,2]))
        vv <-jaccard.test(as.numeric(zz[,1]), as.numeric(zz[,2]), method = "exact")
        Jaccard_distances[i,c(1,2)] <- c(uu,vv$pvalue)
      }
      #print(zz)
    }
    
    assign(paste0 ("Jaccard_distances_U_Lipids_", Tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2), Jaccard_distances)
    write.csv(file=paste0 ("Jaccard_distances_U_Lipids_", Tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2,".csv"), Jaccard_distances)
  }
}    

#### B LIPIDS ###############################################################################################################

# ** Considered_model_1 --------------------------------------------------------------------------------------
# Creates a list (called all_tissues_1) in which each element is 
# a list of lipids coming from the row names of
# each Tissue-specific and model-specific Matrix

all_tissues_1 <-list()
for(k in 1:length(tissues)){
  tissue <- tissues[k]
  aa <- files_names_originals[grep(files_names_originals, pattern=considered_mode)]
  aa <- aa[grep(aa, pattern=considered_generation)]
  aa <- aa[grep(aa, pattern=tissue)]
  
  if(length(aa)!=0){
    yy <- rownames(get(paste0(tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_nozeroes") ) ) 
    assign(paste0("xx_",k), yy)
    all_tissues_1[[k]] <- yy
    
  }
  
}
names(all_tissues_1) <- tissues
if(length(which(all_tissues_1=="NULL")) !=0){
  all_tissues<-all_tissues_1[-which(all_tissues_1=="NULL")]
}

# ** Considered_model_2 --------------------------------------------------------------------------------------
# Creates a list (called all_tissues_2) in which each element is 
# a list of lipids coming from the row names of
# each Tissue-specific and model-specific Matrix

all_tissues_2 <-list()
for(k in 1:length(tissues)){
  tissue <- tissues[k]
  aa <- files_names_originals[grep(files_names_originals, pattern=considered_mode)]
  aa <- aa[grep(aa, pattern=considered_generation)]
  aa <- aa[grep(aa, pattern=tissue)]
  
  if(length(aa)!=0){
    yy <- rownames(get(paste0(tissue, "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_2, "_nozeroes") ) ) 
    assign(paste0("xx_",k), yy)
    all_tissues_2[[k]] <- yy
    
  }
  
}
names(all_tissues_2) <- tissues
if(length(which(all_tissues_2=="NULL")) !=0){
  all_tissues<-all_tissues_2[-which(all_tissues_2=="NULL")]
}

#Creates the IDs of the pairwise comparisons
possible_pairs <- combn(names(all_tissues_1), 2)

#Compares the two models for each pair of tissues
for(d in 1:ncol(possible_pairs)){
  
  tissues_to_compare <- possible_pairs[,d]
  
  B_lipids_1 <- (intersect( all_tissues_1[[tissues_to_compare[1]]] , all_tissues_1[[tissues_to_compare[2]]] ))
  B_lipids_2 <- (intersect( all_tissues_2[[tissues_to_compare[1]]] , all_tissues_2[[tissues_to_compare[2]]] ))
  
  B_lipids_tot <- union(B_lipids_1, B_lipids_2)
  
  B_lipids_matrix<- cbind(B_lipids_tot, rep(0, length(B_lipids_tot)), rep(0, length(B_lipids_tot)) )
  B_lipids_matrix[,2][which(B_lipids_matrix[,1] %in% B_lipids_1)] <- B_lipids_matrix[,1][which(B_lipids_matrix[,1] %in% B_lipids_1)]
  B_lipids_matrix[,3][which(B_lipids_matrix[,1] %in% B_lipids_2)] <- B_lipids_matrix[,1][which(B_lipids_matrix[,1] %in% B_lipids_2)]
  #
  colnames(B_lipids_matrix) <- c("rownames",considered_model_1, considered_model_2)
  rownames(B_lipids_matrix) <- B_lipids_matrix[,1]
  B_lipids_matrix <- B_lipids_matrix[,-1]
  assign(paste0 ("B_Lipids_", "all_", tissues_to_compare[1],  "_", tissues_to_compare[2], "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2), B_lipids_matrix)
  write.csv(file= paste0 ("B_Lipids_", "all_", tissues_to_compare[1],  "_", tissues_to_compare[2], "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2,".csv"), B_lipids_matrix)
  
  
  # Counts the B-lipids in each classes --------------------------------------------------------------------------------------
  
  # Considered_model_1 
  
  classes <- unique(substr(B_lipids_1, start=1, stop =5))
  classes_counts_model_1 <- matrix(ncol=1, nrow=length(classes))
  rownames(classes_counts_model_1) <- classes 
  colnames(classes_counts_model_1) <- paste0("B_lipids_",considered_model_1)
  for(i in 1:length(classes)){
    
    classes_counts_model_1[i] <- length(grep(B_lipids_1, pattern=paste0("^",classes[i])))
  }
  
  # Considered_model_2 
  
  classes <- unique(substr(B_lipids_2, start=1, stop =5))
  classes_counts_model_2 <- matrix(ncol=1, nrow=length(classes))
  rownames(classes_counts_model_2) <- classes 
  colnames(classes_counts_model_2) <- paste0("B_lipids_",considered_model_2)
  for(i in 1:length(classes)){
    
    classes_counts_model_2[i] <- length(grep(B_lipids_2, pattern=paste0("^",classes[i])))
  }
  
  #Creates the B_lipids_classes_counts_tot dataframe, where the vectors classes_counts_model_1  and classes_counts_model_2 are merged
  classes_counts_model_1 <- as.data.frame(classes_counts_model_1)
  classes_counts_model_2 <- as.data.frame(classes_counts_model_2)
  B_lipids_classes_counts_tot <- merge(classes_counts_model_1, classes_counts_model_2, by.x="row.names", by.y="row.names", all=T)
  B_lipids_classes_counts_tot[is.na(B_lipids_classes_counts_tot)]<-0
  rownames(B_lipids_classes_counts_tot) <- B_lipids_classes_counts_tot[,1]
  B_lipids_classes_counts_tot <- B_lipids_classes_counts_tot[,-1]
  #
  All_glicerids <- B_lipids_classes_counts_tot[grep(rownames(B_lipids_classes_counts_tot), pattern="DGX|MGX|TGX"),]
  B_lipids_classes_counts_tot[nrow(B_lipids_classes_counts_tot)+1,] <- colSums(All_glicerids)
  rownames(B_lipids_classes_counts_tot)[nrow(B_lipids_classes_counts_tot)] <- "Glyc"
  
  assign(paste0("B_lipids_classes_counts_", tissues_to_compare[1],  "_", tissues_to_compare[2], "_", considered_generation, "_",  substr(considered_mode, 2,4), "_", considered_model_1, "_", considered_model_2), B_lipids_classes_counts_tot)
  write.csv(B_lipids_classes_counts_tot, file=paste0("B_lipids_classes_counts_", tissues_to_compare[1],  "_", tissues_to_compare[2], "_", considered_generation, "_",  substr(considered_mode, 2,4), "_", considered_model_1, "_", considered_model_2,".csv"))
  
  # Computes the Jaccard Distances --------------------------------------------------------------------------------------
  
  # Unifies the classes MG, DG and TG under the class Glyc, by changing the row names of the B-lipids 
  B_lipids_1_substituted <- gsub("MGXX|DGXX|TGXX", "Glyc", B_lipids_1)
  B_lipids_2_substituted <- gsub("MGXX|DGXX|TGXX", "Glyc", B_lipids_2)
  
  # Merges the vectors containing the A lipids "substituted" for each model into the B_lipids_matrix 
  vv <- union(B_lipids_1_substituted, B_lipids_2_substituted)
  vv<- sort(vv)
  #
  B_lipids_matrix<- cbind(vv, rep(0, length(vv)), rep(0, length(vv)) )
  B_lipids_matrix[,2][which(B_lipids_matrix[,1] %in% B_lipids_1_substituted)] <- B_lipids_matrix[,1][which(B_lipids_matrix[,1] %in% B_lipids_1_substituted)]
  B_lipids_matrix[,3][which(B_lipids_matrix[,1] %in% B_lipids_2_substituted)] <- B_lipids_matrix[,1][which(B_lipids_matrix[,1] %in% B_lipids_2_substituted)]
  #
  colnames(B_lipids_matrix) <- c("rownames",considered_model_1, considered_model_2)
  rownames(B_lipids_matrix) <- B_lipids_matrix[,1]
  B_lipids_matrix <- B_lipids_matrix[,-1]
  assign(paste0 ("B_Lipids_", tissues_to_compare[1],  "_", tissues_to_compare[2], "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2), B_lipids_matrix)
  write.csv(file= paste0 ("B_Lipids_", tissues_to_compare[1],  "_", tissues_to_compare[2], "_",considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2,".csv"), B_lipids_matrix)
  
  # Computes the Jaccard distances
  require(jaccard)
  classes <- rownames(B_lipids_classes_counts_tot)
  Jaccard_distances <- matrix(ncol=2, nrow = length(classes))
  colnames(Jaccard_distances) <- c("Distance", "Pvalue")
  rownames(Jaccard_distances) <- classes
  
  Global_jaccard_matrix <- B_lipids_matrix
  Global_jaccard_matrix[which(Global_jaccard_matrix!=0)] = 1
  
  for(i in 1:length(classes)){  
    zz<- as.matrix(Global_jaccard_matrix[grep(rownames(Global_jaccard_matrix), pattern=classes[i]),])
    if(length(grep(rownames(Global_jaccard_matrix), pattern=classes[i]))==1 ){
      uu<- jaccard(as.numeric(zz[1]), as.numeric(zz[2]))
      vv <-jaccard.test(as.numeric(zz[1]), as.numeric(zz[2]), method = "exact")
      Jaccard_distances[i,c(1,2)] <- c(uu,vv$pvalue)
    }
    else if (length(grep(rownames(Global_jaccard_matrix), pattern=classes[i]))==0) {
      Jaccard_distances[i,c(1,2)] <- c("NA","NA")
    }
    else{
      uu<- jaccard(as.numeric(zz[,1]), as.numeric(zz[,2]))
      vv <-jaccard.test(as.numeric(zz[,1]), as.numeric(zz[,2]), method = "exact")
      Jaccard_distances[i,c(1,2)] <- c(uu,vv$pvalue)
    }
    #print(zz)
  }
  
  assign(paste0 ("Jaccard_distances_B_Lipids_", tissues_to_compare[1],  "_", tissues_to_compare[2], "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2), Jaccard_distances)
  write.csv(file=paste0 ("Jaccard_distances_B_Lipids_", tissues_to_compare[1],  "_", tissues_to_compare[2], "_", considered_generation, "_", substr(considered_mode, 2,4), "_",  considered_model_1, "_",  considered_model_2,".csv"), Jaccard_distances)
}

options(warn=0)

