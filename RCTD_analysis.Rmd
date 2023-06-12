#Robust Cell Type Decomposition (RCTD)n analysis
#Installing spacexr from https://github.com/dmcable/spacexr
```{r}
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
```

```{r}
library(spacexr)
library(Matrix)
library(writexl)
library(ggplot2)
library(ggpubr)
library(writexl)
library(openxlsx)
library(tools)
```


```{r}
#Reference file (expression_Aging_mouse_brain_portal_data_updated.txt) was obtained from "Ximerakis, M., Lipnick, S.L., Innes, B.T. et al. Single-cell transcriptomic profiling of the aging mouse brain. Nat Neurosci 22, 1696–1708 (2019). https://doi.org/10.1038/s41593-019-0491-3"

count_L = read.delim("/path/to/expression_Aging_mouse_brain_portal_data_updated.txt")
count_L_1 = count_L[,-1]
rownames(count_L_1) <- count_L[,1] # Move first column to rownames
count_L_1 = as.data.frame(sapply(count_L_1, as.integer))
rownames(count_L_1) <- count_L[,1]

#load in meta_data (barcodes, clusters, and nUMI)
meta_L <- read.delim("/path/to/meta_Aging_mouse_brain_portal_data.txt") 
meta_L = meta_L[-1,]
cell_types_L <- meta_L$cell_type 
names(cell_types_L) <- meta_L$NAME # create cell_types named list
cell_types_L <- as.factor(cell_types_L) # convert to factor data type
nUMI <- meta_L$nUMI; nUMI = as.numeric(nUMI);names(nUMI) <- meta_L$NAME # create nUMI named list

### Create the Reference object
reference <- Reference(count_L_1, cell_types_L, nUMI)

#Path to filenames
#Sample counts files
sample_counts <- list.files("/path/to/inputdirectory", pattern = "*_counts.txt", full.name = TRUE)
sample_countsNum <- seq(sample_counts)

#Sample coords files
sample_coord <- list.files("/path/to/inputdirectory", pattern = "*_coords.csv", full.name = TRUE)
sample_coordNum <- seq(sample_coord)

resultsdir <- "/path/to/outputdirectory/"
dir.create(resultsdir)

#files matching the pattern
sample_counts
sample_coord
```


#Loop through list of files in path
```{r}
for (fileNumber in sample_countsNum) {
  
  #New name and directory for the output file
  newFileName <-  paste(resultsdir, 
                        file_path_sans_ext(basename(sample_counts[fileNumber])), 
                        ".xlsx", sep = "")
  
  newFileName_jpg <-  paste(resultsdir, 
                        file_path_sans_ext(basename(sample_counts[fileNumber])), 
                        "_", sep = "")
  
  
  #Processing counts file
  counts_L = read.delim(sample_counts[fileNumber])
  counts_1 = counts_L[,-c(1:19)]
  counts_1 = as.data.frame(sapply(counts_1, as.integer))
  rownames(counts_1) = counts_L[,1]
  
  #Processing coords file
  coords_L = read.csv(sample_coord[fileNumber], header = FALSE)
  coords_L$V1 = sub("-", ".", coords_L$V1, fixed = TRUE)
  coords_L_1 = coords_L[,c(5,6)]
  rownames(coords_L_1) = coords_L[,1]
  colnames(coords_L_1)[1] = "xcoord"
  colnames(coords_L_1)[2] = "ycoord"
  nUMI <- colSums(counts_1)
  
  ### Create SpatialRNA object
  puck <- SpatialRNA(coords_L_1, counts_1, nUMI)
  
  ## Examine SpatialRNA object (optional)
  print(dim(puck@counts)) # observe Digital Gene Expression matrix
  hist(log(puck@nUMI,2)) # histogram of log_2 nUMI
  
  print(head(puck@coords)) # start of coordinate data.frame
  barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 

  # This list can be restricted if you want to crop the puck e.g. 
  # puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
  # on the plot:
  plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                     title ='plot of nUMI') 
  
  myRCTD <- create.RCTD(puck, reference, max_cores = 1, test_mode = FALSE, gene_cutoff = 0.000125, fc_cutoff = 0.5, gene_cutoff_reg = 2e-04, 
                        fc_cutoff_reg = 0.75, UMI_min = 0, UMI_max = 2e+07, counts_MIN = 10, UMI_min_sigma = 300, CELL_MIN_INSTANCE=5, 
                        cell_type_names = NULL, keep_reference = F, cell_type_profiles = NULL, CONFIDENCE_THRESHOLD = 5, DOUBLET_THRESHOLD = 20)
  
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  barcodes <- rownames(myRCTD@results$weights)
  cell_types <- colnames(myRCTD@results$weights)
  norm_weights <- normalize_weights(myRCTD@results$weights)
  norm_weights_1 <- data.matrix(norm_weights)
  results_df <- data.frame(barcodes,norm_weights_1)
  colnames(results_df)[2:ncol(results_df)] <- cell_types
  
  write.xlsx(results_df, file=newFileName)
  
  cell_type_names <- myRCTD@cell_type_info$info[[2]]
  spatialRNA <- myRCTD@spatialRNA

  for(i in 1:length(cell_types)){
    plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,cell_type_names[i]], title =cell_type_names[i], size=1.5)
    
    ggsave(paste(newFileName_jpg, cell_type_names[i],'weights.jpg', sep=''), height=5, width=5, units='in', dpi=300)
}
}
```





