#A framework for SPAtial Transcriptomic Analysis
#Installing SPATA2 and dependencies
```{r}
install.packages("devtools")

if (!base::requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
  }

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils', 'EBImage'))

install.packages("Seurat")

devtools::install_github(repo = "kueckelj/confuns")
devtools::install_github(repo = "theMILOlab/SPATAData")
devtools::install_github(repo = "theMILOlab/SPATA2")

devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')

install.packages("msigdbr")
install.packages("hdf5r")
BiocManager::install("GSVA")
```

```{r}
# load packages
library(SPATA2)
library(patchwork)
library(tools)
library(readxl)
```

```{r}
#Path to WAM and DAM files
WAM = read_excel("/path/to/neuron paper smartseq set list.xlsx",sheet = 2)
DAM_1 = read_excel("/path/to/neuron paper smartseq set list.xlsx",sheet = 3)
DAM_2 = read_excel("/path/to/neuron paper smartseq set list.xlsx",sheet = 4)

#Path to directories
out_dirs <- list.dirs("/path/to/spata_loop_test", full.names = TRUE, recursive = FALSE)
out_dirs_samp <- list.dirs("/path/to/spata_loop_test", full.names = FALSE, recursive = FALSE)

out_dirs_samp <- gsub('_out', '', out_dirs_samp)

out_dirs_Num <- seq(out_dirs)
out_dirs_samp_Num <- seq(out_dirs_samp)

#Check if numbers and sample names match
out_dirs
out_dirs_samp
```


```{r}
#Path to WAM and DAM files
WAM = read_excel("/path/to/neuron paper smartseq set list.xlsx", sheet = 2)
DAM_1 = read_excel("/path/to/neuron paper smartseq set list.xlsx", sheet = 3)
DAM_2 = read_excel("/path/to/neuron paper smartseq set list.xlsx", sheet = 4)

#Path to directories
out_dirs <- list.dirs("/path/to/inputdirectory", full.names = TRUE, recursive = FALSE)
#used to shorten directory name down to base name
out_dirs_samp <- list.dirs("/path/to/inputdirectory", full.names = FALSE, recursive = FALSE)

out_dirs_samp <- gsub('_out', '', out_dirs_samp)

out_dirs_Num <- seq(out_dirs)
out_dirs_samp_Num <- seq(out_dirs_samp)

#Output directory to save SPATA object
dir_spata <- "/path/to/output/"

#Check if numbers and sample names match
out_dirs
out_dirs_samp
```



```{r}
#Loop through list of files in path
for (fileNumber in out_dirs_Num) {
  
  #New name and directory for the output file
  newFileName <-  paste(dir_spata, 
                        file_path_sans_ext(basename(out_dirs[fileNumber])))
  
  newFileName_pdf <- paste(dir_spata, 
                        file_path_sans_ext(basename(out_dirs[fileNumber])),
                        ".pdf", sep = "")
  

  #initiate spata object
  spata_obj <- initiateSpataObject_10X(directory_10X = c(out_dirs[fileNumber]),
                                            sample_name = c(out_dirs_samp[fileNumber]))
  
  #save spata object
  saveSpataObject(spata_obj, directory_spata = newFileName, combine_with_wd = FALSE, verbose = NULL)
  
  sample <- loadSpataObject(newFileName)
  
  WAM_1 = WAM$`Gene ID`
  DAM_1_1 = DAM_1$`Gene ID`
  DAM_2_1 = DAM_2$`Gene ID`
  DAM2_2 = DAM_2_1[5231:5331]
  
  input_list <- list(Example1 = WAM_1, 
                   Example2 = DAM_1_1,
                   Example3 = DAM2_2 
                   )
  

  p3 <- plotSurfaceAverage(object = sample, 
                   color_by = input_list, 
                   smooth = FALSE, 
                   pt_size = 2,
                   pt_clrsp = "Reds")
  legendTop()

pdf(file = newFileName_pdf, width = 11, height = 8.5)  
print(p3)
dev.off()
  
}
```


