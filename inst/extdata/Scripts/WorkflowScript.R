
############################
#### set up environment ####
############################
library(MultiOmicsDataCompile); library(data.table)
homedir <- "<path to data base>"
setwd(homedir)

#####################################
#### Obtain platform information ####
#####################################
PlatAnnotInfo <- PlatformAnnotationLoad(PlatInfo=fread(file.path(homedir, "OverviewFiles", "GEOPlatformInfo.csv"), header = TRUE))
fwrite(PlatAnnotInfo, file.path(homedir, "OverviewFiles", "GEOPlatformAnnotation.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
PlatAnnotInfo <- fread(file.path(homedir, "OverviewFiles", "GEOPlatformAnnotation.txt"))

############################
#### Load overview file ####
############################
overview <- fread(file.path(homedir, "OverviewFiles", "GEODataOverview3.csv"), header = TRUE)
overview <- overview[GEO2R == TRUE,]
#### check for duplicates assay names ####
overview$AssayNames <-gsub(" ", "", paste(overview$ID, overview$FCColumnNames, overview$Tissue, overview$Disease, sep = "_"))
overview[duplicated(AssayNames),][!AssayNames == "__",] # check should return nothing

################################
#### Download data from GEO ####
################################
Compiled <- GEOCompile(DS=overview$ID,
                       gpl=overview$GPLNumber,
                       gsm=overview$ComparisonVector,
                       namestr=overview$FCColumnNames,
                       nameraw=overview$RawColumnNames,
                       PlatAnnotInfo = PlatAnnotInfo,
                       destdir = file.path(homedir, "GEOcache"),
                       filename = NULL,
                       writeDB=TRUE,
                       writeRaw=overview$DownloadRawData,
                       GenerateMetaData=overview$DownloadMetaData,
                       MetaDataPath = file.path(homedir, "ArrayMetaData"),
                       writeMetaData=TRUE,
                       DBPath=file.path(homedir, "AppData"),
                       Technology = overview$Technology,
                       renameRaw = FALSE,
                       subsetRaw = FALSE)

#################################################################################
#### double check that the fold change calculations were performed correctly ####
#################################################################################
checkFile <- GEO2RDirectionCheck(DBPath=file.path(homedir, "AppData"),
                                 DS=overview$ID,
                                 namestr=overview$FCColumnNames,
                                 gsm=overview$ComparisonVector,
                                 Technology = overview$Technology,
                                 GraphPath = file.path(homedir, "DirectionCheck"),
                                 subsetRaw = FALSE,
                                 writeRaw=overview$DownloadRawData,
                                 RawQCPath = file.path(homedir, "RawQC"))
checkFile[[1]][correlation < 0,]
fwrite(checkFile[[1]], file.path("./OverviewFiles/ArrayCheck.xls"), row.names = FALSE, quote = FALSE, sep = "\t")

###############################################
#### Externally analyzed RNAseq data setup ####
###############################################
library(org.Hs.eg.db); library(org.Mm.eg.db); library(org.Rn.eg.db)
ExternalDataHarmonize(Fpath = file.path(homedir, "ExternalAnalyzed"),
                      OutPath = file.path(homedir, "AppData"))

####################################################################################################
#### Create or update Summarized experiment objects that integrate all data sets: DEG data only #### Checks for comparisons that are already in the database
#################################################################################################### Duplicated comparisons are not deleted from "AppData"
(DESE <- DESEGenerate(DEGDatapath=file.path(homedir, "AppData"), SEPath = file.path(homedir, "ProcessFiles", "SumarizedExp_DB.rds") ))
(DESE <- readRDS(file.path(homedir, "ProcessFiles", "SumarizedExp_DB.rds")))

###########################
#### Download GTF file ####
###########################
options(timeout=24000); getOption('timeout')
download.file("https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.111.chr.gtf.gz",
              destfile = file.path(homedir, "OverviewFiles", "GTFHuman.gtf.gz"), quiet = FALSE)
download.file("https://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm39.111.chr.gtf.gz",
              destfile = file.path(homedir, "OverviewFiles", "GTFMouse.gtf.gz"), quiet = FALSE)

#########################################
#### Create or update Raw data files #### issue of duplicated genes here probably because of multiple probes measuring the same gene. dealt with this by averaging
#########################################
#### Add Array data ####
RawDataCompile(Fpath = file.path(homedir, "AppData"),
               outPath = file.path("./ProcessFiles/RawData.txt"),
               StartAt = 1,
               sleep = 30,
               GTFHumanFpath = file.path(homedir, "OverviewFiles", "GTFHuman.gtf.gz"),
               GTFMouseFpath = file.path(homedir, "OverviewFiles", "GTFMouse.gtf.gz") )
RawArrayComplete <- fread(file.path("./ProcessFiles/RawData.txt"))

############################################
#### Create meta data for Raw data file #### #### run check ####
############################################
metaDF <- MetDataCompile(RNAseqFilePath= file.path("RunInfo"), ArrayFilePath= file.path("ArrayMetaData"), overview=overview)
#### Save meta data data frame ####
metaDF$rownames <- rownames(metaDF)
fwrite(metaDF, file.path(homedir, "OverviewFiles", "metaData.xls"), row.names = FALSE, quote = FALSE, sep = "\t")
metaDF <- as.data.frame(fread(file.path(homedir, "OverviewFiles", "metaData.xls")))
rownames(metaDF) <- metaDF$rownames
metaDF$rownames <- NULL

#####################################################
#### Create raw data summarizedExperiment object ####
#####################################################
library(SummarizedExperiment)
RawSE <- GenerateRawSE(df = metaDF, ArrayDT=RawArrayComplete, overview=overview)
names(RawSE)
saveRDS(RawSE[["RawSE"]], file=file.path(homedir, "ProcessFiles", "SumarizedExp_RawDB.rds"))
RawSE <- readRDS(file.path(homedir, "ProcessFiles", "SumarizedExp_RawDB.rds"))
assay(RawSE)[1:5,1:5]; head(rowData(RawSE)); head(colData(RawSE))
saveRDS(RawSE[["RPKMSE"]], file=file.path(homedir, "ProcessFiles", "expression_norm.v2.RDS"))
expression_norm <- readRDS(file.path(homedir, "ProcessFiles", "expression_norm.v2.RDS"))
assay(expression_norm)[1:5,1:5]; head(rowData(expression_norm)); head(colData(expression_norm))

#############################################################
#### Update overview file with Externally processed data ####
#############################################################
overview <- fread(file.path(homedir, "OverviewFiles", "GEODataOverview3.csv"), header = TRUE)
overview <- overview[!ID == "",]
IntDT <- data.table(ID = gsub("_.+", "", names(assays(DESE))), FCColumnNames = gsub("^.+?_", "", names(assays(DESE))))
IntDT <- IntDT[!(IntDT$ID %in% overview$ID),] # [grepl("PID[0-9]+", ID),]
IntDT[, `:=`(DataType = "BulkExpressionProfile",	Technology = "RNAseq", DownloadMetaData = FALSE, DownloadRawData = FALSE, GEO2R = FALSE)]
overview <- rbind(overview, IntDT, fill = TRUE)
#### overwrite previous file ####
fwrite(overview, file.path(homedir, "OverviewFiles", "GEODataOverview3.xls"), row.names = FALSE, quote = FALSE, sep = "\t")
#### Update date and other columns by hand ####

################################################################################
################################################################################
#### Proteomic data workflow ###################################################
################################################################################
################################################################################

############################
#### Load overview file ####
############################
overview <- fread(file.path(homedir, "OverviewFiles", "GEODataOverview3.csv"), header = TRUE)
overviewpProteomics <- overview[DataType == "Proteomic",]

##################################
#### Download proteomics data ####
##################################
# devtools::install_version("dbplyr", version = "2.3.4")
# ProteomicsDataDownload(path = file.path(homedir, "Proteomic_1"), DS = overviewpProteomics$ID)

############################################
#### Loop through files and format data ####
############################################
# MZList <- FormatMaxQuant(path = file.path(homedir, "Proteomic_1"))
# names(MZList)

########################################
#### Save formatted proteomics data ####
########################################
# proteomicMZSave(MZList=MZList, path = file.path(homedir, "Proteomic_2"))

###########################################
#### Create experimental design object ####
###########################################
DesignDT <- DesignMatrixFromNames(Fpath = file.path(homedir, "Proteomic_2"))
fwrite(DesignDT, file.path(homedir, "OverviewFiles", "DesignMatrix.xls"), row.names = FALSE, quote = FALSE, sep = "\t")

################################################
#### load into summarized experiment object ####
################################################
library(SummarizedExperiment)
tissueSplitList <- ProtSELoad(DesignDT=DesignDT, Fpath=file.path(homedir, "Proteomic_2"))
assay(tissueSplitList[[2]])
rowData(tissueSplitList[[2]])

##################################################################
#### Perform overall protein abundance assessments ###############
##################################################################
#### compile data together to assess protein expression level ####
TotMel <- RowDataCompile(tissueSplitList=tissueSplitList)
TotMel

#############################################################################################
#### Filter for proteins that are identified in all replicates of at least one condition ####
#############################################################################################
dataFiltList <- DataFilter(dataSeList=tissueSplitList, thr = 0)
assay(dataFiltList[[3]])
rowData(dataFiltList[[3]])

###############################################
#### Perform N-peptides per protein cutoff ####
###############################################
PepCutOffList <- NPeptideThreshold(dataFiltList=dataFiltList, Npeptides = 2)
PepCutOffList[[2]]                  # return the percent of records remaining
dataPepCutOff <- PepCutOffList[[1]] # return the data
colnames(assay(dataPepCutOff[[1]]))

###############################
#### Perform normalization ####
###############################
MultiNormalizeList <- MultiNormalization(dataPepCutOff=dataPepCutOff)
names(MultiNormalizeList)

#### create normalization density plots ####
############################################
densityPlotList <- densityPlotFromList(MultiNormalizeList)
names(densityPlotList)
densityPlotList[["mean"]]
densityPlotList[["median"]]
densityPlotList[["vsn"]]
densityPlotList[["DEPvsn"]]
densityPlotList[["loess"]]
densityPlotList[["rlr"]]
densityPlotList[["smad"]]

###############################################################
#### Select normalization method ##############################
###############################################################
NormalizedSE <- MultiNormalizeList[["median"]]
assay(NormalizedSE[[3]])
rowData(NormalizedSE[[3]])

#############################
#### Impute missing data ####
#############################
#### Determine if any data sets have missing values that need to be imputed ####
(Missing <- DetermineMissing(data=NormalizedSE))
#### Impute missing values ####
set.seed(1)
NormImpAll <- DataImpute(dataFiltList = NormalizedSE, type = "MinProb")

##################################################################
#### remove samples if there are less than 50 samples in them ####
##################################################################
VarRMThreshList <- LowSampleCountRemove(VarRMList=NormImpAll, cut = 50)

#############################################################################
#### calculate Fold changes on imputed, thresholded, and normalized data ####
#############################################################################
FCcut <- log2(1); Pcut <- 0.05; sigCol = "BHCorrection" # "p.val" "p.adj" "BHCorrection"
Comparisons <- fread(file.path(homedir, "OverviewFiles", "ProteomicComparisons.xls"), header = TRUE)
compList2 <- Comparisons$Comparison
names(compList2) <- Comparisons$dataset
compList2 <- compList2[names(compList2) %in% names(VarRMThreshList)]
dataDiffNorm <- DEAnalysis(DataList = VarRMThreshList, type = "manual", ComparisonList=compList2)
dataDiffNorm

###############################
#### Save data to database ####
###############################
SaveToProteomicDB(SEList=dataDiffNorm, Path=file.path(homedir, "Proteomic_3"))

#########################################################################
#### Create a single proteomics DE data summarized experiment object ####
#########################################################################
(DEProtSE <- DESEProtGenerate(ProtDatapath=file.path(homedir, "Proteomic_3"), SEPath = file.path(homedir, "ProcessFiles", "SumarizedProtExp_DB.rds")))
(DEProtSE <- readRDS(file.path(homedir, "ProcessFiles", "SumarizedProtExp_DB.rds")))

###############################
#### Setup Protein name DB ####
###############################
Proteins <- ProteomicProteinName(fPath = "./Proteomic_3")
saveRDS(Proteins, file.path(homedir, "OverviewFiles", "ProteomicProteins.RDS"))

################################################################################
################################################################################
#### Load app ##################################################################
################################################################################
################################################################################
library(shiny)
AppSetup(homedir)
runApp("./Scripts")








