# MAGMA.SPA
### MAGMA GWAS, TWAS, or PWAS Genetic Risk Coexpression Module Integration with _Seyfried Pipeline Adaptation_

Calculates a mean enrichment score for risk in modules or clusters of gene product proteins or trancriptomics, using any comprehensive genome-wide list of genes and their estimated significance of contribution to a trait, e.g. disease risk.  Bootstrap-based permutation calculation is performed as we have published.

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
                "ALS_GWAS_ENSEMBLE_avgMinusLogP(Plt0.05).csv")     #These files must be in MAGMAinputDir

maxP=0.05                 #no genes with a MAGMA summarized p value greater than this will be considered even if in the MAGMA-derived input files.
FDR=0.10                  #FDR or q value (0 < FDR < 1); recommend 0.10, i.e. 10%
barcolors= c("darkslateblue","mediumorchid")  #specify one unique color for each of above MAGMAinputs
                          #common colors: "darkslateblue","mediumorchid", "seagreen3","hotpink","goldenrod","darkorange","darkmagenta", ...
relatednessOrderBar=TRUE  #Plot mean scaled enrichment bar plot in column order (relatedness) of MEs?  If FALSE, they will be plotted in size rank order M1, M2, ...

# Data created during the Seyfried Analysis Pipeline
##################################
NETcolors= net$colors     #module color assignments, vector of length equal to number of rows in cleanDat; should have all colors for modules from 1:minimumSizeRank as printed by WGCNA::labels2colors(1:nModules)
#MEs= MEs                 #Module Eigengenes (or Eigenproteins) with columns of MEs ordered in relatedness order
#cleanDat= cleanDat       #rownames must start with HUMAN gene symbols, separated by any other rowname information using ';' or '|' character

# Other variables
#################################
outFilePrefix="5"         #Filename prefix; step in the pipeline -- for file sorting by name.
outFileSuffix="AD_ALS_GWAS_ensembleAvg"
parallelThreads=8         #Each permutation analysis is run on a separate thread simultaneously, up to this many threads.
calculateMEs=TRUE         #Recalculate MEs and their relatedness order, even if the data already exists.
plotOnly=FALSE            #If plotOnly is TRUE, the variables created by MAGMA.SPA function holding plot data should already exist (xlabels, allBarData).
##################################

# Run the permutation analysis and generate all outputs
source("MAGMA.SPA.R")
MAGMAoutList <- MAGMA.SPA()
# Outputs XLSX, PDF, and list of barplot y values (allBarData), barplot labels (xlabels), and all permutation statistics and gene symbol hits (all_output)


# Rerun function, just to plot previously calculated statistics
allBarData<-MAGMAoutList$allBarData
xlabels<-MAGMAoutList$xlabels
all_output<-MAGMAoutList$all_output
plotOnly=TRUE
MAGMAoutList <- MAGMA.SPA()

```

### <I>Additional Notes</I>

The provided ensemble (multi-GWAS) Alzheimer's disease or ALS risk p values are calculated as the mean -log(p) gene-level risk for all genes reaching nominal significance in any GWAS study considered, following rollup of SNP-level GWAS summary statistics to the gene-level p value using <a href="https://ctg.cncr.nl/software/magma">MAGMA v1.09b</a> command line processing of each GWAS' SNP summary statistics.

MAGMA is only one possible source of genome-wide comprehensive gene-level p values which this bootstrap algorithm can use as input for estimating enrichment of significant risk conferred by multiple gene products represented in a particular systems biology coexpression module or cluster of related gene products.

GWAS studies leveraged to generate the 2-column inputs with human gene SYMBOLS and their associated p value or -log10(p) for disease risk include:

<B>Alzheimer Disease GWAS</B>
1) <a href="https://www.nature.com/articles/s41588-019-0358-2">BW Kunkle et al (2019) AD GWAS</a> (files for this GWAS contain 1834 genes with p<0.05).
2) <a href="https://www.nature.com/articles/ng.2802">J-C Lambert et al (2013) AD GWAS</a> (files for this GWAS contain 1234 genes with p<0.05), and were originally supplemental when used for our publication by <a href="https://www.sciencedirect.com/science/article/pii/S2405471216303702">NT Seyfried et al, Cell Syst (2017)</a>.
3) <a href="https://www.nature.com/articles/s41588-022-01024-z">C Bellenguez et al (2022) AD GWAS</a>, summary statistics available via <a href="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/harmonised/35379992-GCST90027158-MONDO_0004975-Build38.f.tsv.gz">FTP</a> from <a href="https://www.ebi.ac.uk/gwas/studies/GCST90027158">EBI</a>.

<B>Amylotrophic Lateral Sclerosis (ALS) GWAS</B>
1) <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7610013/">A Iacoangelli et al (2020) GWAS</a> (files for this GWAS are labeled by last author "Al Chalibi").
   Related informatics for the paper are on <a href="https://github.com/KHP-Informatics/ALSMetaAnalysis2020">GitHub</a>, and the statistics <a href="http://ec2-user@ec2-35-177-159-177.eu-west-2.compute.amazonaws.com/European_Chinese_MetaAnalysis2020.txt.zip">linked there</a> were downloaded as input for MAGMA.
2) <a href="https://pubmed.ncbi.nlm.nih.gov/34873335/">W van Rheenen et al (2021) GWAS</a>, summary statistics available via <a href="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027164/">FTP</a> from <a href="https://www.ebi.ac.uk/gwas/studies/GCST90027163">EBI</a>.

<I>Preprocessing and linux shell scripts calling MAGMA for each of the above study data are available on request.</I>
