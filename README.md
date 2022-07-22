# MAGMA.SPA
### MAGMA GWAS Genetic Risk Coexpression Module Integration with _Seyfried Pipeline Adaptation_

Requires the following R packages:
 *WGCNA, statmod, doParallel, xlsx, ggplot2, gridBase, grid, gplots, calibrate*

Sample Wrapper for MAGMA.SPA function:

```
#MAGMA-SPA (Seyfried Pipeline Adaptation for MAGMA)
#---------------------------------

# Required parameters, variables, and data must be set as shown above in .GlobalEnv before calling function; currently no defaults are automatic.
##################################
MAGMAinputDir= "E:/5.MAGMAinput/"

MAGMAinputs= c(	"AD_GWAS_ENSEMBLE_averageMinusLogP(Plt0.05).csv",
                "ALS_GWAS_ENSEMBLE_avgMinusLogP(Plt0.05).csv")     #These files exist in MAGMAinputDir

maxP=0.05                 #no genes with a MAGMA summarized p value greater than this will be considered even if in the MAGMA-derived input files.
FDR=0.10		              #FDR or q value (0 < FDR < 1); recommend 0.10, i.e. 10%
barcolors= c("darkslateblue","mediumorchid")  #specify one unique color for each of above MAGMAinputs
                          #common colors: "darkslateblue","mediumorchid", "seagreen3","hotpink","goldenrod","darkorange","darkmagenta", ...
relatednessOrderBar=TRUE  #Plot mean scaled enrichment bar plot in column order (relatedness) of MEs?  If FALSE, they will be plotted in size rank order M1, M2, ...

# Data created during the Seyfried Analysis Pipeline
##################################
NETcolors= net$colors     #module color assignments, vector of length equal to number of rows in cleanDat; should have all colors for modules from 1:minimumSizeRank as printed by WGCNA::labels2colors(1:nModules)
#MEs= MEs                 #Module Eigengenes (or Eigenproteins) with columns of MEs ordered in relatedness order
#cleanDat= cleanDat       #rownames must start with HUMAN gene symbols, separated by any other rowname information using ';' or '|' character

# Other required variables
#################################
outFilePrefix="5"         #Filename prefix; step in the pipeline -- for file sorting by name.
outFileSuffix="AD_ALS_GWAS_ensembleAvg"
parallelThreads=8         #Each permutation analysis is run on a separate thread simultaneously, up to this many threads.
plotOnly=FALSE		        #If plotOnly is TRUE, the variables created by MAGMA.SPA function holding plot data should already exist (xlabels, allBarData).
##################################

# Run the permutation analysis and generate all outputs
source("MAGMA.SPA.R")
MAGMAoutList <- MAGMA.SPA(cleanDat)
# Outputs XLSX, PDF, and list of barplot y values (allBarData), barplot labels (xlabels), and all permutation statistics and gene symbol hits (all_output)


# Rerun function, just to plot previously calculated statistics
allBarData<-MAGMAoutList$allBarData
xlabels<-MAGMAoutList$xlabels
all_output<-MAGMAoutList$all_output
plotOnly=TRUE
MAGMAoutList <- MAGMA.SPA(cleanDat)

