MAGMA.SPA <- function(dummyVar="",env=.GlobalEnv) {

	suppressPackageStartupMessages(require(WGCNA,quietly=TRUE))
	suppressPackageStartupMessages(require(statmod,quietly=TRUE))
	suppressPackageStartupMessages(require(xlsx,quietly=TRUE))
	suppressPackageStartupMessages(require(ggplot2,quietly=TRUE))
	suppressPackageStartupMessages(require(gridBase,quietly=TRUE))
	suppressPackageStartupMessages(require(grid,quietly=TRUE))
	suppressPackageStartupMessages(require(gplots,quietly=TRUE))
	suppressPackageStartupMessages(require(calibrate,quietly=TRUE))

	if (!exists("outFilePrefix")) { outFilePrefix="" } else { if (nchar(outFilePrefix)>0) outFilePrefix=paste0(outFilePrefix,".") }
	if (!exists("outFileSuffix")) { if (exists("FileBaseName")) { outFileSuffix=paste0("-",FileBaseName) } else { outFileSuffix="-unspecified_study" }} else { if (nchar(outFileSuffix)>0) outFileSuffix=paste0("-",outFileSuffix) }

	if (!exists("maxP")) { cat("- maxP variable for maximum p of nominally significant genes to keep for permutation not set. Defaulting to maxP=0.05 ...\n"); maxP=0.05; }
	if (!exists("FDR")) { cat("- FDR variable for q value of significant enrichment of risk in any module not set. Defaulting to FDR=0.10 ...\n"); FDR=0.10; }
	if (!exists("barcolors")) { cat("- barcolors variable not set. Default colors will be used...\n"); barcolors=c("darkslateblue","hotpink","mediumorchid","seagreen3","skyblue","goldenrod","darkorange","darkmagenta","darkred","darkgreen","darkturquoise","saddlebrown","maroon","honeydew","coral","purple","orangered3","lightcoral","cyan","yellow")[1:length(MAGMAinputs)]; }

	if (!exists("relatednessOrderBar") | !is.logical("relatednessOrderBar")) { cat("- relatednessOrderBar not found or not TRUE/FALSE. Ordering modules by relatedness order in MEs.\n"); relatednessOrderBar=TRUE; }
	if (!exists("plotOnly")) { cat("- plotOnly not found. Performing all calculations on provided inputs anew.\n"); plotOnly=FALSE; }
        if (!is.logical(plotOnly)) { cat("- plotOnly not TRUE/FALSE. Performing all calculations on provided inputs anew.\n"); plotOnly=FALSE; }

	if (!plotOnly) {
		if (!exists("cleanDat")) stop("\ncleanDat variable must exist, holding gene product (rows) X sample (columns) data in the form of log2(relative abundance).\n\n")

		if (!exists("NETcolors")) if(exists("net")) { if ("colors" %in% names(net)) { NETcolors=net$colors } else { NETcolors=c() } } else { NETcolors=c() }
		if (!length(NETcolors)==nrow(cleanDat)) { stop("\nNetwork color assignment vector not supplied or not of length in rows of cleanDat; will not be included in output table and data frame.\n\n") }

		MAGMAinputDir=gsub("\\/\\/","/",paste0(MAGMAinputDir,"/"))
		for (file in MAGMAinputs) if (!file.exists(paste0(MAGMAinputDir,file))) stop(paste0("\n",MAGMAinputDir,file," not found. Cannot continue.\n\n"))

		if("MEs" %in% names(net)) MEs=net$MEs
		if(!exists("MEs") | !exists("calculateMEs")) calculateMEs=TRUE
		if(calculateMEs) {
		  cat("- MEs data frame for module eigengenes not found or calculateMEs=TRUE.\n  Attempting to recalculate from cleanDat and net$colors or NETcolors...\n")
		  if(!exists("NETcolors")) if("colors" %in% names(net)) NETcolors=net$colors
		  # if(!exists("cleanDat") | !exists("NETcolors")) stop("Cannot find NETcolors and/or cleanDat for ME calculation.\n")
		  MEs<-tmpMEs<-data.frame()
		  MEList = moduleEigengenes(t(cleanDat), colors = NETcolors)
		  MEs = orderMEs(MEList$eigengenes)
		  colnames(MEs)<-gsub("^ME","",colnames(MEs))  # let's be consistent in case prefix was added, remove it.
		  if("grey" %in% colnames(MEs)) MEs[,"grey"] <- NULL
		}
		colnames(MEs)=gsub("^ME","",colnames(MEs))
		if("grey" %in% colnames(MEs)) MEs[,"grey"] <- NULL

		if(!exists("parallelThreads")) { cat("- parallelThreads variable not set to a number. Attempting to use ",as.numeric(length(MAGMAinputs))," threads for quickest processing...\n"); parallelThreads=as.numeric(length(MAGMAinputs)); }
		if(!is.numeric(as.numeric(parallelThreads))) { cat("- parallelThreads variable not set to a number. Attempting to use ",as.numeric(length(MAGMAinputs))," threads for quickest processing...\n"); parallelThreads=as.numeric(length(MAGMAinputs)); }

		cat(paste0("\nSetting up parallel backend with ",parallelThreads," threads...\n"))
		suppressPackageStartupMessages(require(doParallel, quietly=TRUE))
		# if(exists("clusterLocal")) stopCluster(clusterLocal)
		clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
		registerDoParallel(clusterLocal)

		moduleList=sapply( colnames(MEs),function(x) as.vector(data.frame(do.call("rbind",strsplit(  paste0(data.frame(do.call("rbind",strsplit(rownames(cleanDat),"[|]")))[,1],";")  ,"[;]")))[,1]  )[which(unlist(NETcolors)==x)] )
		nModules=length(names(moduleList))
		for (b in 1:nModules) {
			moduleList[[b]] <- unique(moduleList[[b]][moduleList[[b]] != ""])
			moduleList[[b]] <- unique(moduleList[[b]][moduleList[[b]] != "0"])
		}		
	
		# Order modules by size of modules [as determined using standard ranked colors (no ties)]
		geneList <- list()
		modcolors=unique(unlist(NETcolors))
		modcolors<-modcolors[which(!modcolors=="grey")]
		nModules=length(modcolors)
		modcolors=labels2colors(c(1:nModules))  # INSURE CORRECT RANK ORDER
		for (i in 1:length(modcolors)) 	geneList[[ modcolors[i] ]] <- moduleList[[ modcolors[i] ]]
	
		orderedLabels<- cbind(paste("M",seq(1:nModules),sep=""),labels2colors(c(1:nModules)))
		xlabels.rankOrder <- orderedLabels[,1] 
	
		cat(paste0("Performing boostrap statistics to find mean scaled enrichment scores of significant gene-level risk in ",nModules," modules with ",length(MAGMAinputs)," lists.\n[1 list per each of up to ",parallelThreads," threads at a time]...\n"))
		parallel::clusterExport(cl=clusterLocal, list("maxP","MAGMAinputDir","geneList","xlabels.rankOrder","MAGMAinputs","FDR"), envir=environment())   ## avoid variable not found error during foreach below:
		statOutList <- foreach(thisMAGMAinputFile=as.character(MAGMAinputs)) %dopar% {
	
			## Prepare SNP_data
			SNP_data <- read.csv(paste(MAGMAinputDir,thisMAGMAinputFile,sep=""),stringsAsFactors=FALSE,header=T) #read.table(genePValues, stringsAsFactors = FALSE,header=T)
			SNP_data <- as.data.frame(SNP_data)
			
			#convert csv input's column 2 from p value if not already -log10(p) and filter p<=maxP only
			if(max(SNP_data[,2],na.rm=TRUE)<=1) SNP_data[,2]= -log(SNP_data[,2])
			SNP_data<-SNP_data[which(SNP_data[,2]>= -log10(maxP)),]
			
			nperm <- 10000
			
			# Create function for carrying out permutations for Null distribution
			permute <- function(pVals_SNP,ind){
				pVals_permuted <- sample(pVals_SNP)
				mean(pVals_permuted[ind])
			}
			
			## Permutation p Values can be zero when the significance of the association is very high. The following package (statmod) has a function permp that overcomes this problem and calculates a permutation p value using method described here http://www.statsci.org/webguide/smyth/pubs/permp.pdf
			
			require(statmod)
			commonGenes<-list()
			for (i in 1:nModules){
				pVals_SNP <- as.numeric(SNP_data[c(1:nrow(SNP_data)),2])
				ind <- which(SNP_data[,1] %in% geneList[[i]])
				commonGenes[[i]] <- as.vector(na.omit(intersect(geneList[[i]],SNP_data[,1])))
			
				module_mean <- mean(pVals_SNP[ind])
				module_sd <- sd(pVals_SNP[ind])
				module_sem <- module_sd / sqrt(length(geneList[[i]]))
				permMean <- replicate(nperm,permute(pVals_SNP,ind))
				pVal_module <- sum(abs(permMean) >= abs(module_mean)) / nperm
				numCommon <- length(commonGenes[[i]])
				NES_perm <- (permMean - mean(permMean)) / sd(permMean)
				NES_module <- (module_mean - mean(permMean)) / sd(permMean) 
				NESpVal_module <- sum(abs(NES_perm) >= abs(NES_module)) / nperm
				if (!is.na(NESpVal_module)){
					statmodP <- permp(sum(abs(NES_perm) >= abs(NES_module)),nperm,length(SNP_data[[1]]),length(geneList[[i]]))
				}else{
					statmodP <- "NA"
				}
				if (i == 1){
					mean_allModules <- module_mean
			#		commonGenes_all <- c(commonGenes[[i]])  #,rep(NA,len=50)
					permMean_all <- permMean
					pVal_all <- pVal_module
					numCommon_all <- numCommon
					NES_perm_all <- NES_perm
					NES_module_all <- NES_module
					NESpVal_all <- NESpVal_module
					statmodP_all <- statmodP
				}else{
					mean_allModules <- c(mean_allModules,module_mean)
			#		commonGenes_all <- cbind(commonGenes_all,c(commonGenes[[i]]))  #,rep(NA,len=50)
					permMean_all <- cbind(permMean_all,permMean)
					pVal_all <- cbind(pVal_all,pVal_module)
					numCommon_all <- c(numCommon_all,numCommon)
					NES_perm_all <- cbind(NES_perm_all, NES_perm)
					NES_module_all <- c(NES_module_all, NES_module)
					NESpVal_all <- cbind(NESpVal_all,NESpVal_module)
					statmodP_all <- cbind(statmodP_all,statmodP)
			
				}
				
			}
			
			maxHitListSize<-max(unlist(lapply(commonGenes,length)))
			commonGenes_all1<-lapply(commonGenes,function(x) if (length(x)==maxHitListSize) { sort(x) } else { c(sort(x), rep(NA,maxHitListSize-length(x))) })
			#commonGenes_all<- as.data.frame(matrix(NA,nrow=maxHitListSize,ncol=0))
			#for (i in 1:length(commonGenes_all1)) commonGenes_all<-cbind(commonGenes_all, commonGenes_all1[[i]])
			commonGenes_all <- matrix(unlist(commonGenes_all1),nrow=maxHitListSize,ncol=length(commonGenes_all1),byrow=FALSE)
			
			names(mean_allModules) <- names(geneList)
			names(numCommon_all) <- names(geneList)
			names(NES_module_all) <- names(geneList)
			colnames(commonGenes_all) <- names(geneList)
			colnames(permMean_all) <- names(geneList)
			colnames(pVal_all) <- names(geneList)
			colnames(NES_perm_all) <- names(geneList)
			colnames(NESpVal_all) <- names(geneList)
			colnames(statmodP_all) <- names(geneList)
			
			## Write output to tables and PDF
			pVal_names<-colnames(pVal_all)
			pVal_all<-as.vector(pVal_all)
			names(pVal_all)<-pVal_names
			
			NESpVal_names<-colnames(NESpVal_all)
			NESpVal_all<-as.vector(NESpVal_all)
			names(NESpVal_all)<-NESpVal_names
			
			statmodP_names<-colnames(statmodP_all)
			statmodP_all<-as.vector(statmodP_all)
			names(statmodP_all)<-statmodP_names
			
			all_output <- rbind(mean_allModules, NES_module_all, numCommon_all,
			                    pVal_all, NESpVal_all, statmodP_all,
			                    commonGenes_all)
	#		write.table(all_output,file=paste("./",outFilePrefix,"MAGMA-SPA-",thisMAGMAinputFile,".txt",sep=""), na="", row.names=T,col.names=NA,sep="\t", quote=FALSE)
			
			pdf(paste("./",outFilePrefix,"MAGMA-SPA-",thisMAGMAinputFile,".pdf",sep=""),height=8,width=8)
				par(mfrow=c(4,4))
				par(mar=c(4,4,4,2))
				barplot(mean_allModules,main = "Mean - Enrichment Score", ylab = "Mean",cex.names=0.8, width=0.8,las=2,cex.main=0.95,names.arg=xlabels.rankOrder)
				barplot(NES_module_all,main = "Mean - Scaled Enrichment Score", ylab = "Mean",cex.names=0.8, width=0.8,las=2,cex.main=0.95,names.arg=xlabels.rankOrder)
				
				modColors <- names(geneList)
				for (a in 1:nModules){
					modClr <- xlabels.rankOrder[a]
					NESpValue_module <- NESpVal_all[a]
					if (is.na(NESpValue_module)){
						next
					} else{
						par(mar=c(4,4,4,4))
						permMean_Module <- NES_perm_all[,a]
						hist(permMean_Module,xlab = "Normalized ES (Mean)", ylab = "Frequency", main=modClr, sub = paste("pValue",NESpValue_module),col="grey",freq=T, border="black",cex.main=0.95)
						abline(v=NES_module_all[a],col="red")	
					}
				}
			dev.off()
	
	
			return(list(all_output, NES_module_all,thisMAGMAinputFile))
		}
	
		# re-combine list elements from two outputs, from all MAGMAinputFile runs
		all_output            = do.call(list, lapply(statOutList,function(x){x[[1]]}) )
		NES_module_all        = do.call(list, lapply(statOutList,function(x){x[[2]]}) )
		names(all_output)<-names(NES_module_all) <- MAGMAinputs <- do.call(c,lapply(statOutList,function(x){x[[3]]}))
	
		all_output <- lapply(all_output,function(x) { rownames(x)[7:nrow(x)]<- paste0("geneHit",1:(nrow(x)-6)); x; })
	
		##Export all gene lists to multi-sheet Excel 
		my.xlsx=createWorkbook()
		for (sheetName in MAGMAinputs) addDataFrame(as.data.frame(all_output[[sheetName]]), sheet=createSheet(my.xlsx, sheetName), startColumn=1, row.names=TRUE,col.names=TRUE)
		saveWorkbook(my.xlsx, paste0("./",outFilePrefix,"MAGMA-SPA",outFileSuffix,".xlsx"))
	
	
		#################
		## Summary plot of multiple GWAS MAGMA bootstrap statistics
	
		allBarData<-do.call(rbind,NES_module_all)  #(AD.IGAP.1234=AD1234_NES_all,AD.2019=AD.Kunkle_NES_all,ASD=ASD_NES_all,SCZ=SCZ_NES_all,PD=PD_NES_all,ALS=ALS_NES_all)
		rownames(allBarData)<-MAGMAinputs
		allBarData<-allBarData[,c(1:nModules)]
		allBarData[!is.finite(allBarData)]<-0
	
		xlabels = if (relatednessOrderBar) { orderedLabels[match(colnames(MEs),orderedLabels[,2]),1] } else { xlabels.rankOrder }
		# handle single GWAS list input case, and reorder module bars if relatednessOrderBar==TRUE   # Feb 2, 2023 bugfix.
		if (as.numeric(length(MAGMAinputs))==1) { allBarData<-data.frame(oneInput=allBarData[if (relatednessOrderBar) { colnames(MEs) } else { 1:ncol(MEs) }]); colnames(allBarData)=as.character(MAGMAinputs)[1]; allBarData=t(allBarData); }  else { if (relatednessOrderBar) allBarData <- allBarData[,colnames(MEs)] }

	} #end if (!plotOnly)

	# Check that plot data exists, in case plotOnly==TRUE
	if (!exists("xlabels") | !exists("allBarData")) stop("Processed MAGMA Enrichment in Modules not found. Cannot plot.\nRerun with plotOnly=FALSE.\n\n")

##Simplest output
#	pdf(paste0("./",outFilePrefix,"MAGMA-Enr_BarPlot",outFileSuffix,".pdf"),height=8,width=16)
#	 par(mfrow=c(1,1))
#	 par(mar=c(4,6,4,2))
#	 barplot(allBarData,main =paste0("Enrichment of MAGMA-implicated Genetic Risk of Disease(s)\nbased on ",length(MAGMAinputs)," GWAS-derived MAGMA summary p (<=",maxP,") gene lists,\n(as published in Seyfried, et al, Cell Systems, 2017)"), ylab = "Mean-scaled Enr. Score",cex.names=0.70, cex.lab=1.75, width=0.8,las=2,cex.main=0.95,
#		 names.arg=xlabels, beside=TRUE,
#                 col=barcolors,
#                 legend.text=TRUE, args.legend=list(x="top", bty="n", inset=c(0, 0)), xpd=FALSE,ylim=c(min(na.omit(allBarData))-1,max(na.omit(allBarData))+1.3))
#	 abline(h=1.28, lty=2, col="red")
#	dev.off()

## Output with inset for color swatches for x axis labels 
	
	## Create ordered swatches for inset at x axis using ggplot bar
	 colorData=data.frame(Mnum=as.character(xlabels),yBlank=c(1), fill=gplots::col2hex(unique(WGCNA::labels2colors(as.numeric(gsub("M","",as.character(xlabels)))))), fillName=unique(WGCNA::labels2colors(as.numeric(gsub("M","",as.character(xlabels))))))
	# barplot(as.matrix(diff(1:93)), horiz=T, col=colorData$fillName, axes=F, xlab=NA)
	 p <- ggplot(colorData, aes(x=Mnum, y=yBlank)) + geom_bar(stat="identity", fill=colorData$fill, color="#000000", size=0.01, width = 0.8675, aes(fill=fillName)) +theme(axis.title = element_blank(), axis.text.x = element_text(face="bold", color="#000000", size=12, angle=90), axis.ticks.y=element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.background = element_rect(fill = "transparent",colour = NA)) +
	      scale_x_discrete( limits=colorData$Mnum)  + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.3), axis.text.y=element_text(color="#FFFFFF")) + labs(x="", y="") + scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0,1), breaks=c(0), label=c(""))
	
	
	pdf(file=paste0("./",outFilePrefix,"MAGMA-Enr_BarPlot",outFileSuffix,".pdf"),width=16,height=9, onefile=FALSE)  # onefile=FALSE for forcing no blank page 1.
	
	grid.newpage()
#	grid.rect(gp=gpar(col="black"))
	#grid.text("",
	#  y=unit(1, "npc") - unit(1, "lines"), gp=gpar(col="black"))
	#grid.rect(gp=gpar(col="green"))
	
	pushViewport(viewport(width=1, height=1, name="OuterFrame", just=c("center","center")))
	
	par(mar=c(5,6,4,2))
	plot.new()
	
	#	 par(omi = gridFIG(), new = TRUE)
		 par(fig = gridFIG(), new = TRUE)
		 opar <- par(lwd = 0.01)
	
		 barplot(allBarData,main =paste0("Enrichment of MAGMA-implicated Genetic Risk of Disease(s)\nbased on ",length(MAGMAinputs)," GWAS-derived MAGMA summary p (<=",maxP,") gene lists,\n(as published in Seyfried, et al, Cell Systems, 2017)"), ylab = "Mean-scaled Enr. Score",cex.names=0.70, cex.lab=1.75, width=0.8,las=2,cex.main=0.95,
			 names.arg=rep("",length(xlabels)),  #xlabels,
			 beside=TRUE, space=c(0,rep(c(rep(0,length(MAGMAinputs)-1),0.15*length(MAGMAinputs)),length(xlabels)-1),rep(0,length(MAGMAinputs)-1)), #border=NA,
	                 col=barcolors,
	                 legend.text=TRUE, args.legend=list(x="top", bty="n", inset=c(0, 0)), xpd=FALSE,ylim=c(min(na.omit(allBarData))-1,max(na.omit(allBarData))+1.3))
		 abline(h=qnorm(FDR,lower=F), lty=3, lwd=1, col="red")
		 calibrate::textxy(X=0,Y=qnorm(FDR,lower=F),paste0("FDR = ",FDR*100,"%"),cex=1)
	#upViewport()
	
	#For onscreen: pushViewport(viewport(width=0.92245, height=0.10, x=unit(9.255, "inch"), y=0.035)) #ok for out to screen, 
	#For 93 modules:
	#pushViewport(viewport(width=0.8498, height=0.08, x=unit(8.381, "inch"), y=0.0925)) #for out to PDF
	#For 23 modules:
	#pushViewport(viewport(width=0.8618, height=0.08, x=unit(8.381, "inch"), y=0.0925)) #for out to PDF
	#Automatic adjustment of viewport inset for swatches (works for 43-93 modules)
	#pushViewport(viewport(width=0.879-0.00031*length(xlabels), height=0.08, x=0.524, y=0.0925, just=c("center","center"))) #for out to PDF
	#(works for 13, 43, 93 modules)
	 pushViewport(viewport(width=0.879-(0.00031-0.00000571*(length(xlabels)-93))*length(xlabels), height=0.08, x=0.5235, y=0.0925, just=c("center","center"))) #for out to PDF
	 par(fig = gridFIG(), new = TRUE)
	#  just="bottom", name="Swatches"))
	# grid.rect(gp=gpar(col="blue"))
	 print(p, newpage=FALSE)
	
	dev.off()

return(list(allBarData=allBarData,xlabels=xlabels, all_output=all_output))
}  # end Function MAGMA.SPA
						 
