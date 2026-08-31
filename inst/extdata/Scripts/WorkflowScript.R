
############################
#### set up environment ####
############################
library(MultiOmicsDataCompile); library(data.table)
homedir <- "<path to data base>"

#####################################
#### Obtain platform information ####
#####################################
PlatAnnotInfo <- PlatformAnnotationLoad(PlatInfo=fread(file.path(homedir, "OverviewFiles", "GEOPlatformInfo.csv"), header = TRUE))
fwrite(PlatAnnotInfo, file.path(homedir, "OverviewFiles", "GEOPlatformAnnotation.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
PlatAnnotInfo <- fread(file.path(homedir, "OverviewFiles", "GEOPlatformAnnotation.txt"))

############################
#### Load overview file ####
############################
overview_all <- fread(file.path(homedir, "OverviewFiles", "GEODataOverview3.csv"), header = TRUE)
overview_all <- as.data.table(ValidateDatasetManifest(overview_all))
overview <- overview_all[GEO2R == TRUE,]
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
                       Species = overview$Species,
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
fwrite(checkFile[[1]], file.path(homedir, "OverviewFiles", "ArrayCheck.xls"), row.names = FALSE, quote = FALSE, sep = "\t")

###############################################
#### Externally analyzed RNAseq data setup ####
###############################################
library(org.Hs.eg.db); library(org.Mm.eg.db); library(org.Rn.eg.db)
ExternalDataHarmonize(Fpath = file.path(homedir, "ExternalAnalyzed"),
                      OutPath = file.path(homedir, "AppData"))

##############################################################################
#### Rebuild the differential-expression SummarizedExperiment from AppData ####
##############################################################################
(DESE <- DESEGenerate(DEGDatapath=file.path(homedir, "AppData"), SEPath = file.path(homedir, "ProcessFiles", "SumarizedExp_DB.rds") ))
(DESE <- readRDS(file.path(homedir, "ProcessFiles", "SumarizedExp_DB.rds")))

#########################################
#### Create or update raw data files ####
#########################################
RawDataCompile(Fpath = file.path(homedir, "AppData"),
               outPath = file.path(homedir, "ProcessFiles", "RawData.txt"),
               StartAt = 1,
               overview = overview_all)
RawArrayComplete <- fread(file.path(homedir, "ProcessFiles", "RawData.txt"))

##################################################################
#### Derive the sample conditions required by GenerateRawSE() ####
##################################################################
#### Conditions come from the manifest ComparisonVector, which positions one
#### character per sample of the stored raw matrix, so they are taken from an
#### explicit comparison code rather than inferred from sample order. Review
#### this table before continuing: samples that every comparison excludes are
#### labelled "Unassigned".
conditions <- SampleConditions(DBPath = file.path(homedir, "AppData"),
                               overview = overview_all)
fwrite(conditions, file.path(homedir, "OverviewFiles", "SampleConditions.xls"), row.names = FALSE, quote = FALSE, sep = "\t")

############################################
#### Create meta data for Raw data file #### #### run check ####
############################################
metaDF <- MetDataCompile(
  RNAseqFilePath = file.path(homedir, "RunInfo"),
  ArrayFilePath = file.path(homedir, "ArrayMetaData"),
  overview = overview_all,
  conditions = conditions
)
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
ExpressionCollection <- GenerateRawSE(df = metaDF, ArrayDT=RawArrayComplete, overview=overview_all)
names(ExpressionCollection)
saveRDS(ExpressionCollection, file=file.path(homedir, "ProcessFiles", "ExpressionCollection.rds"))
saveRDS(ExpressionCollection[["ExplorationSE"]], file=file.path(homedir, "ProcessFiles", "expression_norm.v2.RDS"))
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
Proteins <- ProteomicProteinName(fPath = file.path(homedir, "Proteomic_3"))
saveRDS(Proteins, file.path(homedir, "OverviewFiles", "ProteomicProteins.RDS"))

################################################################################
################################################################################
#### Load app ##################################################################
################################################################################
################################################################################
launchMultiOmicsExplorer(homedir)








