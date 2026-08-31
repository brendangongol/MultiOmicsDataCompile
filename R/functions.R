#' Create the MultiOmicsDataCompile directory structure
#'
#' Creates the directories used by the processing workflow and copies the
#' packaged application resources into them.
#'
#' @param homedir A path to the directory where the database will be stored.
#'
#' @return Invisibly, a logical vector indicating whether each resource copy
#'   succeeded.
#' @family setup functions
#' @importFrom R.utils copyDirectory
#' @export
makeDirectory <- function(homedir){
  dirs <- c(homedir,
            file.path(homedir, "ArrayMetaData"),
            file.path(homedir, "ProcessFiles"),
            file.path(homedir, "AppData"),
            file.path(homedir, "GEOcache"),
            file.path(homedir, "DirectionCheck"),
            file.path(homedir, "RawQC"),
            file.path(homedir, "Proteomic_1"),
            file.path(homedir, "Proteomic_3") )
  vapply(dirs, dir.create, logical(1), showWarnings = FALSE,
         recursive = TRUE)
  SysFpath <- system.file("extdata", package = "MultiOmicsDataCompile")
  from <- file.path(SysFpath, list.files(SysFpath))
  to <- file.path(homedir, list.files(SysFpath))
  copied <- vapply(seq_along(from), function(x){
    R.utils::copyDirectory(from=from[x], to=to[x], recursive=TRUE)
    TRUE
  }, logical(1))
  invisible(copied)
}

#' Compile microarray platform annotations
#'
#' Builds a common probe-to-gene annotation table for the GEO platforms listed
#' in `PlatInfo`. Annotation data are obtained from Bioconductor annotation
#' packages or GEO, depending on the platform.
#'
#' @param PlatInfo A data frame or data table with a `GPLNumber` column
#'   containing GEO platform accessions.
#'
#' @return A data table with probe `ID`, `ENTREZID`, `SYMBOL`, `GENENAME`, and
#'   `GPLID` columns.
#' @family transcriptomics functions
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @export
PlatformAnnotationLoad <- function(PlatInfo){
  AnnotationDT <- data.table()
  for(i in seq_len(nrow(PlatInfo))){
    if(PlatInfo$GPLNumber[i] == "GPL23126"){
      # BiocManager::install("clariomdhumantranscriptcluster.db")
      # library(clariomdhumantranscriptcluster.db)
      x <- .annotation_object("clariomdhumantranscriptcluster.db", "clariomdhumantranscriptclusterGENENAME")
      mapped_probes<-mappedkeys(x)
      GENENAME <- unlist(as.list(x[mapped_probes]))
      DTNames <- data.table(ID=names(GENENAME), GENENAME)
      x <- .annotation_object("clariomdhumantranscriptcluster.db", "clariomdhumantranscriptclusterENTREZID")
      mapped_probes<-mappedkeys(x)
      ENTREZ<-  unlist(as.list(x[mapped_probes]))
      DTENTREZ <- data.table(ID=names(ENTREZ), ENTREZID = ENTREZ)
      mer <- merge(DTNames, DTENTREZ, by = "ID")
      symbols <- as.character(mer$ENTREZ)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      GPL23126 <- merge(mer, goHuman, by = "ENTREZID")
      GPL23126 <- GPL23126[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL23126$GPLID <- "GPL23126"
      AnnotationDT <- rbind(AnnotationDT, GPL23126)
    } else if(PlatInfo$GPLNumber[i] == "GPL30511"){
      z <- getGEO("GPL30511")
      temp <- as.data.table(z@dataTable@table)
      setnames(temp, "Symbol", "SYMBOL")
      symbols <- temp$SYMBOL
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("ENTREZID"),"SYMBOL"))
      GPL30511 <- merge(temp, goHuman, by = "SYMBOL")
      GPL30511 <- GPL30511[,c("SYMBOL", "ID", "ENTREZID"), with = FALSE]
      GPL30511 <- GPL30511[!is.na(ENTREZID),]
      GPL30511$ENTREZID <- as.character(GPL30511$ENTREZID)
      Values <- GPL30511$ENTREZID
      mappings <- AnnotationDbi::select(org.Hs.eg.db, keys=Values, columns=c("GENENAME"),keytype="ENTREZID")
      GPL30511 <- merge(GPL30511, mappings, by = "ENTREZID")
      GPL30511 <- GPL30511[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL30511$GPLID <- "GPL30511"
      AnnotationDT <- rbind(AnnotationDT, GPL30511)
    } else if(PlatInfo$GPLNumber[i] == "GPL28577"){
      z <- getGEO("GPL28577")
      temp <- as.data.table(z@dataTable@table)
      symbols <- temp$ID
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("ENTREZID"),"SYMBOL"))
      GPL28577 <- merge(temp, goHuman, by.x = "ID", by.y = "SYMBOL", all.x = TRUE)
      GPL28577[, SYMBOL := ID]
      GPL28577 <- GPL28577[!is.na(ENTREZID),]
      GPL28577$ENTREZID <- as.character(GPL28577$ENTREZID)
      Values <- GPL28577$ENTREZID
      mappings <- AnnotationDbi::select(org.Hs.eg.db, keys=Values, columns=c("GENENAME"),keytype="ENTREZID")
      GPL28577 <- merge(GPL28577, mappings, by = "ENTREZID")
      GPL28577 <- GPL28577[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL28577$GPLID <- "GPL28577"
      AnnotationDT <- rbind(AnnotationDT, GPL28577)
    } else if(PlatInfo$GPLNumber[i] == "GPL29503"){
      z <- getGEO("GPL29503")
      temp <- as.data.table(z@dataTable@table)
      symbols <- temp$ID
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("ENTREZID", "GENENAME"),"SYMBOL"))
      GPL29503 <- merge(temp, goHuman, by.x = "ID", by.y = "SYMBOL", all.x = TRUE)
      GPL29503[, SYMBOL := ID]
      GPL29503 <- GPL29503[,c("ID", "ENTREZID", "SYMBOL", "GENENAME"), with = FALSE]
      GPL29503$GPLID <- "GPL29503"
      AnnotationDT <- rbind(AnnotationDT, GPL29503)
    } else if(PlatInfo$GPLNumber[i] == "GPL14951"){
      z <- getGEO("GPL14951")
      temp <- as.data.table(z@dataTable@table)
      setnames(temp, c("Symbol", "Entrez_Gene_ID"), c("SYMBOL", "ENTREZID"))
      GPL14951 <- temp
      GPL14951 <- GPL14951[,c("ID", "SYMBOL", "ENTREZID"), with = FALSE]
      GPL14951 <- unique(GPL14951[!is.na(ENTREZID),])
      GPL14951 <- GPL14951[!is.na(ENTREZID),]
      GPL14951$ENTREZID <- as.character(GPL14951$ENTREZID)
      Values <- GPL14951$ENTREZID
      mappings <- AnnotationDbi::select(org.Hs.eg.db, keys=Values, columns=c("GENENAME"),keytype="ENTREZID")
      GPL14951 <- merge(GPL14951, mappings, by = "ENTREZID")
      GPL14951 <- GPL14951[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL14951$GPLID <- "GPL14951"
      AnnotationDT <- rbind(AnnotationDT, GPL14951)
    } else if(PlatInfo$GPLNumber[i] == "GPL6244"){
      # BiocManager::install("hugene10sttranscriptcluster.db")
      # library(hugene10sttranscriptcluster.db)
      x <- .annotation_object("hugene10sttranscriptcluster.db", "hugene10sttranscriptclusterGENENAME")
      mapped_probes<-mappedkeys(x)
      GENENAME <- unlist(as.list(x[mapped_probes]))
      DTNames <- data.table(ID=names(GENENAME), GENENAME)
      x <- .annotation_object("hugene10sttranscriptcluster.db", "hugene10sttranscriptclusterENTREZID")
      mapped_probes<-mappedkeys(x)
      ENTREZ<-  unlist(as.list(x[mapped_probes]))
      DTENTREZ <- data.table(ID=names(ENTREZ), ENTREZID = ENTREZ)
      mer <- merge(DTNames, DTENTREZ, by = "ID")
      symbols <- as.character(mer$ENTREZ)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      GPL6244 <- merge(mer, goHuman, by = "ENTREZID")
      GPL6244 <- GPL6244[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL6244$GPLID <- "GPL6244"
      AnnotationDT <- rbind(AnnotationDT, GPL6244)
    } else if(PlatInfo$GPLNumber[i] == "GPL16686"){
      # BiocManager::install("hugene20sttranscriptcluster.db")
      # library(hugene20sttranscriptcluster.db)
      x <- .annotation_object("hugene20sttranscriptcluster.db", "hugene20sttranscriptclusterGENENAME")
      mapped_probes<-mappedkeys(x)
      GENENAME <- unlist(as.list(x[mapped_probes]))
      DTNames <- data.table(ID=names(GENENAME), GENENAME)
      x <- .annotation_object("hugene20sttranscriptcluster.db", "hugene20sttranscriptclusterENTREZID")
      mapped_probes<-mappedkeys(x)
      ENTREZ<-  unlist(as.list(x[mapped_probes]))
      DTENTREZ <- data.table(ID=names(ENTREZ), ENTREZID = ENTREZ)
      mer <- merge(DTNames, DTENTREZ, by = "ID")
      symbols <- as.character(mer$ENTREZ)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      GPL16686 <- merge(mer, goHuman, by = "ENTREZID")
      GPL16686 <- GPL16686[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL16686$GPLID <- "GPL16686"
      AnnotationDT <- rbind(AnnotationDT, GPL16686)
    } else if(PlatInfo$GPLNumber[i] == "GPL570"){
      z <- getGEO("GPL570")
      temp <- as.data.table(z@dataTable@table)
      temp <- temp[,c("ID", "Gene Symbol", "ENTREZ_GENE_ID"), with = FALSE]
      temp$`Gene Symbol` <- gsub(" ///.+", "", temp$`Gene Symbol`)
      temp$ENTREZ_GENE_ID <- gsub(" ///.+", "", temp$ENTREZ_GENE_ID)
      temp <- temp[!temp$`Gene Symbol` == "",]
      setnames(temp, c("Gene Symbol", "ENTREZ_GENE_ID"), c("SYMBOL", "ENTREZID"))
      GPL570 <- temp
      GPL570 <- GPL570[!is.na(ENTREZID),]
      GPL570 <- GPL570[!ENTREZID == "",]
      GPL570$ENTREZID <- as.character(GPL570$ENTREZID)
      Values <- GPL570$ENTREZID
      mappings <- AnnotationDbi::select(org.Hs.eg.db, keys=Values, columns=c("GENENAME"),keytype="ENTREZID")
      GPL570 <- unique(merge(GPL570, mappings, by = "ENTREZID", allow.cartesian = TRUE))
      GPL570 <- GPL570[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL570$GPLID <- "GPL570"
      AnnotationDT <- rbind(AnnotationDT, GPL570)
    } else if(PlatInfo$GPLNumber[i] == "GPL11532"){
      # BiocManager::install("hugene11sttranscriptcluster.db")
      # library(hugene11sttranscriptcluster.db)
      x <- .annotation_object("hugene11sttranscriptcluster.db", "hugene11sttranscriptclusterGENENAME")
      mapped_probes<-mappedkeys(x)
      GENENAME <- unlist(as.list(x[mapped_probes]))
      DTNames <- data.table(ID=names(GENENAME), GENENAME)
      x <- .annotation_object("hugene11sttranscriptcluster.db", "hugene11sttranscriptclusterENTREZID")
      mapped_probes<-mappedkeys(x)
      ENTREZ<-  unlist(as.list(x[mapped_probes]))
      DTENTREZ <- data.table(ID=names(ENTREZ), ENTREZID = ENTREZ)
      mer <- merge(DTNames, DTENTREZ, by = "ID")
      symbols <- as.character(mer$ENTREZ)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      GPL11532 <- merge(mer, goHuman, by = "ENTREZID")
      GPL11532 <- GPL11532[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL11532$GPLID <- "GPL11532"
      AnnotationDT <- rbind(AnnotationDT, GPL11532)
    } else if(PlatInfo$GPLNumber[i] == "GPL14877"){
      # BiocManager::install("hgu133plus2.db")
      # library(hgu133plus2.db)
      x <- .annotation_object("hgu133plus2.db", "hgu133plus2GENENAME")
      mapped_probes<-mappedkeys(x)
      GENENAME <- unlist(as.list(x[mapped_probes]))
      DTNames <- data.table(ID=names(GENENAME), GENENAME)
      x <- .annotation_object("hgu133plus2.db", "hgu133plus2ENTREZID")
      mapped_probes<-mappedkeys(x)
      ENTREZ<-  unlist(as.list(x[mapped_probes]))
      DTENTREZ <- data.table(ID=names(ENTREZ), ENTREZID = ENTREZ)
      mer <- merge(DTNames, DTENTREZ, by = "ID")
      symbols <- as.character(mer$ENTREZ)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      goHuman <- unique(goHuman)
      GPL14877 <- merge(mer, goHuman, by = "ENTREZID")
      GPL14877 <- GPL14877[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL14877$GPLID <- "GPL14877"
      AnnotationDT <- rbind(AnnotationDT, GPL14877)
    } else if(PlatInfo$GPLNumber[i] == "GPL2895"){
      # BiocManager::install("hwgcod.db")
      # library(hwgcod.db)
      x <- .annotation_object("hwgcod.db", "hwgcodGENENAME")
      mapped_probes<-mappedkeys(x)
      GENENAME <- unlist(as.list(x[mapped_probes]))
      DTNames <- data.table(ID=names(GENENAME), GENENAME)
      x <- .annotation_object("hwgcod.db", "hwgcodENTREZID")
      mapped_probes<-mappedkeys(x)
      ENTREZ<-  unlist(as.list(x[mapped_probes]))
      DTENTREZ <- data.table(ID=names(ENTREZ), ENTREZID = ENTREZ)
      mer <- merge(DTNames, DTENTREZ, by = "ID")
      symbols <- as.character(mer$ENTREZ)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      goHuman <- unique(goHuman)
      GPL2895 <- merge(mer, goHuman, by = "ENTREZID")
      GPL2895 <- GPL2895[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL2895$ID <- as.character(GPL2895$ID)
      GPL2895$ID <- gsub("GE", "", GPL2895$ID)
      GPL2895$GPLID <- "GPL2895"
      AnnotationDT <- rbind(AnnotationDT, GPL2895)
    } else if(PlatInfo$GPLNumber[i] == "GPL10558"){
      z <- getGEO("GPL10558")
      temp <- as.data.table(z@dataTable@table)
      setnames(temp, c("ID", "Entrez_Gene_ID"), c("ID", "ENTREZID"))
      symbols <- as.character(temp$ENTREZID)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL", "GENENAME"),"ENTREZID" ) )
      goHuman <- unique(goHuman)
      goHuman$ENTREZID <- as.character(goHuman$ENTREZID)
      temp$ENTREZID <- as.character(temp$ENTREZID)
      GPL10558 <- merge(temp, goHuman, by = "ENTREZID")
      GPL10558 <- GPL10558[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL10558 <- GPL10558[!is.na(ENTREZID),]
      GPL10558$GPLID <- "GPL10558"
      AnnotationDT <- rbind(AnnotationDT, GPL10558)
    } else if(PlatInfo$GPLNumber[i] == "GPL14550"){
      z <- getGEO("GPL14550")
      z <- z@dataTable@table[, c("ID", "GENE", "GENE_SYMBOL", "GENE_NAME")]
      z <- z[!z$GENE_NAME == "",]
      setnames(z, c("ID", "GENE", "GENE_SYMBOL", "GENE_NAME"), c("ID", "ENTREZID", "SYMBOL", "GENENAME"))
      GPL14550 <- as.data.table(z)
      GPL14550 <- GPL14550[,c("ID", "ENTREZID", "SYMBOL", "GENENAME"), with = FALSE]
      GPL14550$GPLID <- "GPL14550"
      AnnotationDT <- rbind(AnnotationDT, GPL14550)
    } else if(PlatInfo$GPLNumber[i] == "GPL17586"){
      # BiocManager::install("hta20transcriptcluster.db")
      # library(hta20transcriptcluster.db)
      x <- .annotation_object("hta20transcriptcluster.db", "hta20transcriptclusterGENENAME")
      mapped_probes<-mappedkeys(x)
      GENENAME <- unlist(as.list(x[mapped_probes]))
      DTNames <- data.table(ID=names(GENENAME), GENENAME)
      x <- .annotation_object("hta20transcriptcluster.db", "hta20transcriptclusterENTREZID")
      mapped_probes<-mappedkeys(x)
      ENTREZ<-  unlist(as.list(x[mapped_probes]))
      DTENTREZ <- data.table(ID=names(ENTREZ), ENTREZID = ENTREZ)
      mer <- merge(DTNames, DTENTREZ, by = "ID")
      symbols <- as.character(mer$ENTREZ)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      goHuman <- unique(goHuman)
      GPL17586 <- merge(mer, goHuman, by = "ENTREZID")
      GPL17586 <- GPL17586[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL17586$GPLID <- "GPL17586"
      AnnotationDT <- rbind(AnnotationDT, GPL17586)
    } else if(PlatInfo$GPLNumber[i] == "GPL8910"){
      z <- getGEO("GPL8910")
      temp <- as.data.table(z@dataTable@table)
      Values <- unique(temp$`Reference Accession`)
      mappings <- AnnotationDbi::select(org.Hs.eg.db, keys=Values, columns=c("ENTREZID","SYMBOL","GENENAME"),keytype="REFSEQ")
      t2 <- temp[,c("ID", "Reference Accession"), with = FALSE] %>% setnames("Reference Accession", "REFSEQ")
      GPL8910 <- merge(t2, mappings, by = "REFSEQ")
      GPL8910 <- GPL8910[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL8910$GPLID <- "GPL8910"
      AnnotationDT <- rbind(AnnotationDT, GPL8910)
    } else if(PlatInfo$GPLNumber[i] == "GPL127"){
      z <- getGEO("GPL127")
      temp <- as.data.table(z@dataTable@table)
      setnames(temp, c("ID", "GENE"), c("ID", "ENTREZID"))
      symbols <- as.character(temp$ENTREZID)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      goHuman <- unique(goHuman)
      goHuman$ENTREZID <- as.character(goHuman$ENTREZID)
      temp$ENTREZID <- as.character(temp$ENTREZID)
      GPL127 <- merge(temp, goHuman, by = "ENTREZID", all.x = TRUE)
      GPL127 <- GPL127[,c("ENTREZID", "ID", "GENE_NAME", "SYMBOL"), with = FALSE]
      setnames(GPL127, "GENE_NAME", "GENENAME")
      GPL127 <- GPL127[!ENTREZID == "",]
      GPL127 <- GPL127[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL127$GPLID <- "GPL127"
      AnnotationDT <- rbind(AnnotationDT, GPL127)
    } else if(PlatInfo$GPLNumber[i] == "GPL128"){
      z <- getGEO("GPL128")
      temp <- as.data.table(z@dataTable@table)
      setnames(temp, c("ID", "GENE"), c("ID", "ENTREZID"))
      symbols <- as.character(temp$ENTREZID)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      goHuman <- unique(goHuman)
      goHuman$ENTREZID <- as.character(goHuman$ENTREZID)
      temp$ENTREZID <- as.character(temp$ENTREZID)
      GPL128 <- merge(temp, goHuman, by = "ENTREZID", all.x = TRUE)
      GPL128 <- GPL128[,c("ENTREZID", "ID", "GENE_NAME", "SYMBOL"), with = FALSE]
      setnames(GPL128, "GENE_NAME", "GENENAME")
      GPL128 <- GPL128[!ENTREZID == "",]
      GPL128 <- GPL128[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL128$GPLID <- "GPL128"
      AnnotationDT <- rbind(AnnotationDT, GPL128)
    } else if(PlatInfo$GPLNumber[i] == "GPL131"){
      z <- getGEO("GPL131")
      temp <- as.data.table(z@dataTable@table)
      setnames(temp, c("ID", "GENE"), c("ID", "ENTREZID"))
      symbols <- as.character(temp$ENTREZID)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      goHuman <- unique(goHuman)
      goHuman$ENTREZID <- as.character(goHuman$ENTREZID)
      temp$ENTREZID <- as.character(temp$ENTREZID)
      GPL131 <- merge(temp, goHuman, by = "ENTREZID", all.x = TRUE)
      GPL131 <- GPL131[,c("ENTREZID", "ID", "GENE_NAME", "SYMBOL"), with = FALSE]
      setnames(GPL131, "GENE_NAME", "GENENAME")
      GPL131 <- GPL131[!ENTREZID == "",]
      GPL131 <- GPL131[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL131$GPLID <- "GPL131"
      AnnotationDT <- rbind(AnnotationDT, GPL131)
    } else if(PlatInfo$GPLNumber[i] == "GPL549"){
      z <- getGEO("GPL549")
      temp <- as.data.table(z@dataTable@table)
      setnames(temp, c("ID", "GENE"), c("ID", "ENTREZID"))
      symbols <- as.character(temp$ENTREZID)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"),"ENTREZID" ) )
      goHuman <- unique(goHuman)
      goHuman$ENTREZID <- as.character(goHuman$ENTREZID)
      temp$ENTREZID <- as.character(temp$ENTREZID)
      GPL549 <- merge(temp, goHuman, by = "ENTREZID", all.x = TRUE)
      GPL549 <- GPL549[,c("ENTREZID", "ID", "Gene Name", "SYMBOL"), with = FALSE]
      setnames(GPL549, "Gene Name", "GENENAME")
      GPL549 <- GPL549[!ENTREZID == "",]
      GPL549 <- GPL549[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL549$GPLID <- "GPL549"
      AnnotationDT <- rbind(AnnotationDT, GPL549)
    } else if(PlatInfo$GPLNumber[i] == "GPL96"){
      z <- getGEO("GPL96")
      temp <- as.data.table(z@dataTable@table)
      GPL96 <- temp[,c("ID", "Gene Symbol", "Gene Title", "ENTREZ_GENE_ID"), with = FALSE]
      setnames(GPL96, c("Gene Symbol", "Gene Title", "ENTREZ_GENE_ID"), c("SYMBOL", "GENENAME", "ENTREZID"))
      GPL96$ENTREZID <- gsub(" /// .+", "", GPL96$ENTREZID)
      GPL96$SYMBOL <- gsub(" /// .+", "", GPL96$SYMBOL)
      GPL96 <- GPL96[!ENTREZID == "",]
      GPL96 <- GPL96[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL96$GPLID <- "GPL96"
      AnnotationDT <- rbind(AnnotationDT, GPL96)
    } else if(PlatInfo$GPLNumber[i] == "GPL3050"){
      z <- getGEO("GPL3050")
      temp <- as.data.table(z@dataTable@table)
      temp <- temp[,c("ID", "ENSEMBL_ID"), with = FALSE]
      setnames(temp, c("ID", "ENSEMBL_ID"), c("ID", "ENSEMBL"))
      symbols <- as.character(temp$ENSEMBL)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL", "ENTREZID"), "ENSEMBL") )
      goHuman <- unique(goHuman)
      GPL3050 <- merge(temp, goHuman, by = "ENSEMBL")
      GPL3050 <- unique(GPL3050[,!c("ENSEMBL"), with = FALSE])
      GPL3050 <- GPL3050[!is.na(ENTREZID),]
      Values <- unique(GPL3050$ENTREZID)
      mappings <- AnnotationDbi::select(org.Hs.eg.db, keys=Values, columns=c("GENENAME"),keytype="ENTREZID")
      GPL3050 <- merge(GPL3050, mappings, by = "ENTREZID")
      GPL3050 <- GPL3050[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL3050$GPLID <- "GPL3050"
      AnnotationDT <- rbind(AnnotationDT, GPL3050)
    } else if(PlatInfo$GPLNumber[i] == "GPL20115"){
      z <- getGEO("GPL20115")
      temp <- as.data.table(z@dataTable@table)
      temp <- temp[, c("ID", "EntrezGeneID", "GeneName"), with = FALSE]
      setnames(temp, c("EntrezGeneID", "GeneName"), c("ENTREZID", "GENENAME"))
      temp$ENTREZID <- as.character(temp$ENTREZID)
      temp <- temp[!is.na(ENTREZID),]
      symbols <- as.character(temp$ENTREZID)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("SYMBOL"), "ENTREZID") )
      goHuman <- unique(goHuman)
      GPL20115 <- unique(merge(temp, goHuman, by = "ENTREZID"))
      GPL20115 <- GPL20115[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL20115$GPLID <- "GPL20115"
      AnnotationDT <- rbind(AnnotationDT, GPL20115)
    } else if(PlatInfo$GPLNumber[i] == "GPL18056"){
      z <- getGEO("GPL18056")
      temp <- as.data.table(z@dataTable@table)
      temp <- temp[,c("ID", "Gene Symbol", "Gene Name"), with = FALSE]
      setnames(temp, c("ID", "Gene Symbol", "Gene Name"), c("ID", "SYMBOL", "GENENAME"))
      symbols <- as.character(temp$SYMBOL)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, symbols, c("ENTREZID"), "SYMBOL") )
      goHuman <- unique(goHuman)
      GPL18056 <- unique(merge(temp, goHuman, by = "SYMBOL"))
      GPL18056 <- GPL18056[!is.na(ENTREZID),]
      GPL18056 <- GPL18056[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL18056$GPLID <- "GPL18056"
      AnnotationDT <- rbind(AnnotationDT, GPL18056)
    } else if(PlatInfo$GPLNumber[i] == "GPL19886"){
      z <- getGEO("GPL19886")
      temp <- as.data.table(z@dataTable@table)
      temp <- temp[, c("ID", "Gene symbol", "EntrezGeneID"), with = FALSE]
      setnames(temp, c("ID", "Gene symbol", "EntrezGeneID"), c("ID", "SYMBOL", "ENTREZID"))
      temp$ENTREZID <- as.character(temp$ENTREZID)
      Values <- temp$ENTREZID
      mappings <- AnnotationDbi::select(org.Hs.eg.db, keys=Values, columns=c("GENENAME"),keytype="ENTREZID")
      GPL19886 <- merge(temp, mappings, by = "ENTREZID")
      GPL19886 <- GPL19886[,c("ID", "ENTREZID", "SYMBOL", "GENENAME")]
      GPL19886$GPLID <- "GPL19886"
      AnnotationDT <- rbind(AnnotationDT, GPL19886)
    } else if (PlatInfo$GPLNumber[i] == "GPL23159"){
      z <- getGEO("GPL23159")
      z <- as.data.table(z@dataTable@table)
      REFSEQ <- gsub(" ", "", gsub("//.+", "" , z[[10]]))
      z$REFSEQ <- REFSEQ
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, REFSEQ, c("REFSEQ", "ENTREZID", "SYMBOL", "GENENAME"),"REFSEQ"))
      goHuman <- unique(goHuman[!is.na(ENTREZID),])
      z <- z[, c("ID", "REFSEQ"), with = FALSE]
      z <- z[!grepl("--CONTROL|--normgene", REFSEQ, ignore.case = TRUE),]
      GPL23159 <- merge(z, goHuman, by = "REFSEQ")
      GPL23159 <- GPL23159[,c("ID", "ENTREZID", "SYMBOL", "GENENAME"), with = FALSE]
      GPL23159$GPLID <- "GPL23159"
      AnnotationDT <- rbind(AnnotationDT, GPL23159)
    } else if (PlatInfo$GPLNumber[i] == "GPL19109"){
      z <- getGEO("GPL19109")
      z <- as.data.table(z@dataTable@table)
      setnames(z, "ENTREZ_GENE_ID", "ENTREZID")
      z <- z[!is.na(ENTREZID),]
      z$ENTREZID <- as.character(z$ENTREZID)
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, z$ENTREZID, c("SYMBOL", "GENENAME"),"ENTREZID"))
      goHuman <- unique(goHuman[!is.na(ENTREZID),])
      GPL19109 <- unique(merge(z, goHuman, by = "ENTREZID"))
      GPL19109 <- GPL19109[,c("ID", "ENTREZID", "SYMBOL", "GENENAME"), with = FALSE]
      GPL19109$GPLID <- "GPL19109"
      AnnotationDT <- rbind(AnnotationDT, GPL19109)
    } else if (PlatInfo$GPLNumber[i] == "GPL17692"){
      z <- getGEO("GPL17692")
      z <- as.data.table(z@dataTable@table)
      z <- z[,c("ID", "GB_ACC"), with = FALSE][grepl("NR_|NM_", GB_ACC),]
      setnames(z, "GB_ACC", "REFSEQ")
      z <- z[!grepl("--CONTROL|--normgene", REFSEQ, ignore.case = TRUE),]
      REFSEQ <- z$REFSEQ
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, REFSEQ, c("REFSEQ", "ENTREZID", "SYMBOL", "GENENAME"),"REFSEQ"))
      goHuman <- unique(goHuman[!is.na(ENTREZID),])
      GPL17692 <- merge(z, goHuman, by = "REFSEQ")
      GPL17692 <- GPL17692[,c("ID", "ENTREZID", "SYMBOL", "GENENAME"), with = FALSE]
      GPL17692$GPLID <- "GPL17692"
      AnnotationDT <- rbind(AnnotationDT, GPL17692)
    } else if (PlatInfo$GPLNumber[i] == "GPL13667"){
      z <- getGEO("GPL13667")
      z <- as.data.table(z@dataTable@table)
      z <- z[,c("ID", "Entrez Gene"), with = FALSE]
      z <- z[!grepl("---", `Entrez Gene`, ignore.case = TRUE),]
      z <- z[!duplicated(`Entrez Gene`),]
      z$`Entrez Gene` <- as.character(z$`Entrez Gene`)
      setnames(z, "Entrez Gene", "ENTREZID")
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, z$ENTREZID, c("ENTREZID", "SYMBOL", "GENENAME"),"ENTREZID"))
      goHuman <- goHuman[!is.na(SYMBOL),]
      GPL13667 <- merge(z, goHuman, by = "ENTREZID")
      GPL13667 <- GPL13667[,c("ID", "ENTREZID", "SYMBOL", "GENENAME"), with = FALSE]
      GPL13667$GPLID <- "GPL13667"
      AnnotationDT <- rbind(AnnotationDT, GPL13667)
    } else if (PlatInfo$GPLNumber[i] == "GPL10335"){
      z <- getGEO("GPL10335")
      temp <- as.data.table(z@dataTable@table)
      setnames(temp, c("ID", "GeneID"), c("ID", "ENTREZID"))
      temp <- temp[, c("ID", "ENTREZID"), with = FALSE]
      goHuman <- as.data.table( AnnotationDbi::select(org.Hs.eg.db, temp$ENTREZID, c("ENTREZID", "SYMBOL", "GENENAME"),"ENTREZID"))
      goHuman <- goHuman[!is.na(SYMBOL),]
      GPL10335 <- merge(temp, goHuman, by = "ENTREZID")
      GPL10335 <- GPL10335[,c("ID", "ENTREZID", "SYMBOL", "GENENAME"), with = FALSE]
      GPL10335$GPLID <- "GPL10335"
      AnnotationDT <- rbind(AnnotationDT, GPL10335)
      } else {
      print(paste("There is no annotation information available for:", PlatInfo$GPLNumber[i])) } }
  return(AnnotationDT)
}

#' Download and process GEO expression data
#'
#' Performs two-group differential-expression analysis for microarray or RNA-seq
#' GEO series and optionally saves harmonized differential-expression, raw, and
#' sample-metadata files.
#'
#' @param DS A character vector of GEO series accessions (GSE identifiers).
#' @param gpl A character vector of GEO platform accessions (GPL identifiers).
#' @param gsm A character vector of comparison codes. Each character identifies
#'   a control (`0`), treatment (`1`), or excluded (`X`) sample.
#' @param namestr A character vector containing comma-separated output column
#'   names for fold change, p-value, and adjusted p-value results.
#' @param nameraw An optional character vector containing comma- or
#'   semicolon-separated names for raw-data columns.
#' @param PlatAnnotInfo A platform annotation data table returned by
#'   [PlatformAnnotationLoad()].
#' @param destdir A path used to cache GEO downloads.
#' @param filename An optional character vector of local GEO filenames. Use
#'   `NULL` to let [GEOquery::getGEO()] retrieve each series.
#' @param writeDB A logical value indicating whether differential-expression
#'   results should be written to `DBPath`.
#' @param writeRaw A logical vector indicating whether raw data should be
#'   written for each dataset.
#' @param GenerateMetaData A logical vector indicating whether sample metadata
#'   should be generated for each dataset.
#' @param MetaDataPath A path where sample-metadata files are written.
#' @param writeMetaData A logical value indicating whether generated metadata
#'   should be written to disk.
#' @param DBPath A path where raw and differential-expression files are written.
#' @param Technology A character vector containing `"Array"` or `"RNAseq"` for
#'   each dataset.
#' @param Species An optional character vector containing the organism for each
#'   dataset. For RNA-seq data, the value is validated against the species and
#'   genome information returned by GEO when available.
#' @param shrinkLFC A logical value indicating whether RNA-seq log2 fold changes
#'   should be shrunk with DESeq2's normal-prior estimator for ranking and
#'   visualization. Statistical tests and adjusted p-values are unchanged.
#' @param renameRaw A logical value indicating whether raw-data columns should
#'   be renamed using `nameraw`.
#' @param subsetRaw A logical value indicating whether raw data should retain
#'   only samples used in the comparison.
#' @param sleep The number of seconds to pause between datasets.
#'
#' @return Invisibly, `NULL`. Results are written to the requested output
#'   directories.
#' @family transcriptomics functions
#' @importFrom stats model.matrix
#' @importFrom limma lmFit
#' @importFrom limma makeContrasts
#' @importFrom limma contrasts.fit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 results
#' @importFrom DESeq2 lfcShrink
#' @importFrom GEOquery getRNASeqData
#' @importFrom GEOquery getRNASeqQuantGenomeInfo
#' @importFrom GEOquery getGEO
#' @export
GEOCompile <- function(DS, gpl, gsm, namestr, nameraw, PlatAnnotInfo, destdir, filename=NULL,
                       writeDB=TRUE, writeRaw=TRUE, GenerateMetaData, MetaDataPath, writeMetaData=TRUE,
                       DBPath, Technology, Species = NULL, renameRaw=FALSE,
                       subsetRaw=FALSE, sleep = 0, shrinkLFC = TRUE){
  .validate_parallel_lengths(
    DS, gpl, gsm, namestr, nameraw, writeRaw, GenerateMetaData, Technology,
    names = c("gpl", "gsm", "namestr", "nameraw", "writeRaw",
              "GenerateMetaData", "Technology")
  )
  n_datasets <- length(DS)
  gpl <- .recycle_to_length(gpl, n_datasets)
  gsm <- .recycle_to_length(gsm, n_datasets)
  namestr <- .recycle_to_length(namestr, n_datasets)
  nameraw <- .recycle_to_length(nameraw, n_datasets)
  writeRaw <- .recycle_to_length(writeRaw, n_datasets)
  GenerateMetaData <- .recycle_to_length(GenerateMetaData, n_datasets)
  Technology <- .recycle_to_length(Technology, n_datasets)
  if (!is.null(Species)) {
    .validate_parallel_lengths(DS, Species, names = "Species")
    Species <- .recycle_to_length(.normalize_species(Species), n_datasets)
    if (anyNA(Species)) {
      stop("`Species` contains an unsupported organism.", call. = FALSE)
    }
  }
  for(a in seq_along(DS)){

    if(Technology[a] == "Array"){
      nam <- namestr[a]
      nam <- gsub(" ", "", nam)
      nam <- strsplit(nam, ",")[[1]]
      gset <- suppressMessages(getGEO(DS[a], filename = filename[a], destdir = destdir))
      if(length(gset) > 1) idx <- grep(gpl[a], attr(gset, "names")) else idx <- 1
      gset <- gset[[idx]]
      fvarLabels(gset) <- make.names(fvarLabels(gset))
      comp <- gsub(" ", "", gsm[a])
      comp <- gsub(",", "", comp)
      gsms <- paste0(comp)
      #### Set up raw names ####
      if(renameRaw){
        namraw1 <- nameraw[a]
        namraw1 <- gsub(" ", "", namraw1)
        namraw1 <- make.unique(strsplit(namraw1, ",|;")[[1]])
        namraw1 <- gsub("\n", "", namraw1)
        namraw1 <- namraw1[!strsplit(gsms, split="")[[1]] == "X"] }
      sml <- c()
      for(i in 1:nchar(gsms)){ sml[i] <- substr(gsms,i,i)}
      ex <- exprs(gset)
      comparison <- .parse_comparison_code(gsm[a], colnames(ex), DS[a])
      selected_labels <- comparison$labels[comparison$selected]
      qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
      LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] < 2)
      if(LogC){ ex[which(ex <= 0)] <- NaN
      exprs(gset) <- log2(ex) }
      gset_analysis <- gset[, comparison$selected]
      f1 <- factor(paste0("G", selected_labels), levels = c("G0", "G1"))
      gset_analysis$description2 <- f1
      design <- model.matrix(~description2 + 0, gset_analysis)
      colnames(design) <- levels(f1)
      fit <- lmFit(gset_analysis, design)
      cont.matrix <- makeContrasts(G1-G0, levels = design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2, 0.01)
      tT <- topTable(
        fit2, adjust.method = "fdr", sort.by = "B", number = Inf
      )
      #### subset ####
      ex2 <- data.table(subset(tT, select=c("ID", "logFC", "P.Value", "adj.P.Val")))
      ex2$ID <- as.character(ex2$ID)
      #### annotate with gene names ####
      plat <- PlatAnnotInfo[GPLID == gpl[a],][,!"GPLID", with = FALSE]
      if(nrow(plat) == 0){ print(paste("There is no annotation information available for", gpl[a])) }
      plat$ID <- as.character(plat$ID)
      plat <- plat[!duplicated(plat$ID),]
      ex2 <- merge(plat, ex2, by = "ID")
      ex2 <- ex2[,.SD[which.min(P.Value)], by = "SYMBOL"] # duplicated gene names are dealt with by taking the most significant record
      setnames(ex2, c("logFC", "P.Value", "adj.P.Val"), nam)
      ex2$ID <- as.character(ex2$ID)
      exraw <- data.table(ex)
      if(subsetRaw){ exraw <- exraw[,grepl("[0-9]", strsplit(as.character(gsms), split = "")[[1]]), with = FALSE] }
      if(renameRaw){ setnames(exraw, colnames(exraw), namraw1)
      } else { setnames(exraw, colnames(exraw), paste(DS[a], colnames(exraw), sep = "_")) }
      exraw$ID <- as.character(rownames(ex))
      #### annotate raw data with gene names ####
      exraw <- merge(plat, exraw, by = "ID")
      array_species <- if (is.null(Species)) NA_character_ else Species[a]
      provenance <- .expression_provenance(
        DS[a], species = array_species, platform = gpl[a],
        source = "NCBI GEO series matrix"
      )
      attr(ex2, "provenance") <- provenance
      attr(exraw, "provenance") <- provenance
      if(writeDB){
        fpath <- file.path(DBPath, .de_file_name(DS[a], namestr[a]))
        saveRDS(ex2, file = fpath) }
      if(writeRaw[a]){
        if(subsetRaw){
       fpath <- file.path(DBPath, .raw_file_name(DS[a], namestr[a], TRUE))
        } else {
          fpath <- file.path(DBPath, .raw_file_name(DS[a], namestr[a], FALSE))
        }
       saveRDS(exraw, file = fpath) }
      if(GenerateMetaData[a]){
        Pdat <- pData(gset)
        #### Add Meta data ####
        colMet <- as.data.table(Pdat)
        if(writeMetaData){
          fpath <- file.path(MetaDataPath, paste(paste(DS[a], gsub("-", "_", gsub("_.+", "", nam[1])), "MetaData", sep = "_"), ".txt", sep = ""))
          saveRDS(colMet, file = gsub(".txt", ".rds", fpath) ) }
        }
      }

    if(Technology[a] == "RNAseq"){
      #### Set up DEG names ####
      nam <- namestr[a]
      nam <- gsub(" ", "", nam)
      nam <- strsplit(nam, ",")[[1]]
      if (length(nam) != 3L) {
        stop(
          "`namestr` for ", DS[a],
          " must contain exactly three names: logFC, Pvalue, and AdjPValue.",
          call. = FALSE
        )
      }
      #### Set up raw names ####
      if(renameRaw){
        namraw1 <- nameraw[a]
        namraw1 <- gsub(" ", "", namraw1)
        namraw1 <- make.unique(strsplit(namraw1, ",|;")[[1]])
        namraw1 <- gsub("\n", "", namraw1) }
      # GEOquery discovers the current NCBI-generated raw-count and annotation
      # URLs, including their species and genome build, instead of relying on a
      # hard-coded human GRCh38 filename.
      rna_se <- GEOquery::getRNASeqData(DS[a])
      tbl_all <- as.matrix(SummarizedExperiment::assay(rna_se, "counts"))
      comparison <- .parse_comparison_code(gsm[a], colnames(tbl_all), DS[a])
      sml <- comparison$labels[comparison$selected]
      sel <- which(comparison$selected)
      exraw <- tbl_all
      tbl <- tbl_all[, sel, drop = FALSE]
      annot <- .standardize_geo_rnaseq_annotation(
        SummarizedExperiment::rowData(rna_se),
        rownames(tbl_all)
      )
      genome_info <- tryCatch(
        GEOquery::getRNASeqQuantGenomeInfo(DS[a]),
        error = function(e) c(species = NA_character_, genome = NA_character_)
      )
      source_species <- unname(genome_info[grepl("species", names(genome_info), ignore.case = TRUE)][1])
      source_genome <- unname(genome_info[grepl("genome|build", names(genome_info), ignore.case = TRUE)][1])
      if (!length(source_species)) source_species <- NA_character_
      if (!length(source_genome)) source_genome <- NA_character_
      if (!is.null(Species) && !is.na(source_species)) {
        normalized_source_species <- .normalize_species(source_species)
        if (!is.na(normalized_source_species) && Species[a] != normalized_source_species) {
          stop(
            "Species mismatch for ", DS[a], ": manifest says ", Species[a],
            " but GEO reports ", source_species, ".",
            call. = FALSE
          )
        }
      }
      gs <- factor(sml)
      groups <- make.names(c("Ctrl", "Tx"))
      levels(gs) <- groups
      sample_info <- data.frame(Group = gs, row.names = colnames(tbl))
      keep <- rowSums( tbl >= 10 ) >= min(table(gs))
      tbl <- tbl[keep, ]
      ds <- DESeqDataSetFromMatrix(countData=tbl, colData=sample_info, design= ~Group)
      ds <- DESeq(ds, test="Wald", sfType="poscount")
      r <- results(ds, contrast=c("Group", groups[2], groups[1]), alpha=0.05, pAdjustMethod ="fdr")
      mle_logfc <- stats::setNames(as.numeric(r$log2FoldChange), rownames(r))
      if (isTRUE(shrinkLFC)) {
        shrunk <- DESeq2::lfcShrink(ds, coef = 2L, type = "normal")
        r$log2FoldChange <- shrunk$log2FoldChange
      }
      tT <- as.data.frame(r[order(r$padj, na.last = TRUE), ])
      tT$log2FoldChangeMLE <- mle_logfc[match(rownames(tT), names(mle_logfc))]
      tT$ENTREZID <- rownames(tT)
      tT <- merge(tT, annot, by = "ENTREZID", all.x = TRUE, sort = FALSE)
      #### subset ####
      ex2 <- data.table::as.data.table(
        tT[, c("ENTREZID", "SYMBOL", "GENENAME", "log2FoldChange",
               "log2FoldChangeMLE", "pvalue", "padj")]
      )
      #### Adjust columnnames and order ####
      setnames(ex2, c("log2FoldChange", "pvalue", "padj"), nam)
      mle_name <- paste0(nam[1], "_MLE")
      setnames(ex2, "log2FoldChangeMLE", mle_name)
      ex2$ID <- ex2$ENTREZID
      ex2 <- ex2[,c("SYMBOL", "ID", "ENTREZID", "GENENAME", nam, mle_name), with = FALSE]
      #### Get Raw data ####
      if(subsetRaw){ exraw <- SummarizedExperiment::assay(ds)  }
      GeneID <- rownames(exraw)
      exraw <- as.data.table(exraw)
      #### Update column names ####
      if(renameRaw){
        namraw1 <- namraw1[gsub(".+_", "", namraw1) %in% colnames(exraw)]
        if(identical(colnames(exraw),gsub(".+_", "", namraw1))){
          setnames(exraw, colnames(exraw), namraw1)
        } else { print(paste("new colnames are not equivalent for", DS[a], nam[1])) }
      } else { setnames(exraw, colnames(exraw), paste(DS[a], colnames(exraw), sep = "_")) }
      exraw$ENTREZID <- as.character(GeneID)
      #### merge FC and raw data together ####
      raw_annotation <- data.table::as.data.table(annot, keep.rownames = FALSE)
      merRaw <- merge(raw_annotation, exraw, by = "ENTREZID", all.y = TRUE)
      provenance <- .expression_provenance(
        DS[a], species = source_species, genome = source_genome,
        source = "NCBI GEO generated RNA-seq counts"
      )
      attr(ex2, "provenance") <- provenance
      attr(merRaw, "provenance") <- provenance
      if(writeDB){
        fpath <- file.path(DBPath, .de_file_name(DS[a], namestr[a]))
        saveRDS(ex2, file = fpath) }
      if(writeRaw[a]){
        if(subsetRaw){
        fpath <- file.path(DBPath, .raw_file_name(DS[a], namestr[a], TRUE))
        } else {
          fpath <- file.path(DBPath, .raw_file_name(DS[a], namestr[a], FALSE))
        }
        saveRDS(merRaw, file = fpath) }
      if(GenerateMetaData[a] && writeMetaData){
        colMet <- data.table::as.data.table(
          as.data.frame(SummarizedExperiment::colData(rna_se)),
          keep.rownames = "geo_accession"
        )
        rna_metadata_path <- file.path(MetaDataPath, "RNAseq")
        dir.create(rna_metadata_path, recursive = TRUE, showWarnings = FALSE)
        fpath <- file.path(
          rna_metadata_path,
          paste(DS[a], gsub("-", "_", gsub("_.+", "", nam[1])),
                "MetaData", sep = "_")
        )
        saveRDS(colMet, file = paste0(fpath, ".rds"))
      }
      }
    print(paste("Completed", a, "of", length(DS)))
    if (sleep > 0) Sys.sleep(sleep) }
  }

#' Validate differential-expression direction
#'
#' Recomputes treatment-control summaries from saved expression data, records
#' their agreement with the compiled fold changes, and produces diagnostic
#' plots.
#'
#' @param DBPath A path containing raw and differential-expression data.
#' @param DS A character vector of GEO series accessions.
#' @param namestr A character vector containing the output column names used by
#'   [GEOCompile()].
#' @param gsm A character vector of control, treatment, and exclusion codes.
#' @param Technology A character vector containing `"Array"` or `"RNAseq"` for
#'   each dataset.
#' @param GraphPath A path where direction-check plots are written.
#' @param subsetRaw A logical value indicating whether the stored raw data were
#'   subset to samples in the comparison.
#' @param writeRaw A logical vector indicating which datasets have raw data.
#' @param RawQCPath A path where raw-data quality-control plots are written.
#'
#' @return A list with `checkFile`, a data table of validation statistics, and
#'   `plotList`, the generated plot objects.
#' @family transcriptomics functions
#' @export
GEO2RDirectionCheck <- function(DBPath, DS, namestr, gsm, Technology, GraphPath, subsetRaw = FALSE, writeRaw, RawQCPath){
  .validate_parallel_lengths(
    DS, namestr, gsm, Technology, writeRaw,
    names = c("namestr", "gsm", "Technology", "writeRaw")
  )
  namestr <- .recycle_to_length(namestr, length(DS))
  gsm <- .recycle_to_length(gsm, length(DS))
  Technology <- .recycle_to_length(Technology, length(DS))
  writeRaw <- .recycle_to_length(writeRaw, length(DS))
  dir.create(GraphPath, recursive = TRUE, showWarnings = FALSE)
  dir.create(RawQCPath, recursive = TRUE, showWarnings = FALSE)
  plotList <- list()
  plotnames <- NULL
  checkFile <- data.table()
  for(a in seq_along(DS)){
    if (!isTRUE(writeRaw[a])) {
      message("Skipping ", DS[a], " because no raw data were requested.")
      next
    }
    #### load data ####
    nam <- namestr[a]
    nam <- gsub(" ", "", nam)
    nam <- strsplit(nam, ",")[[1]]
    if(subsetRaw){
      rawfpath <- file.path(DBPath, .raw_file_name(DS[a], namestr[a], TRUE))
    } else {
      rawfpath <- file.path(DBPath, .raw_file_name(DS[a], namestr[a], FALSE))
    }
    fcfpath <- file.path(DBPath, .de_file_name(DS[a], namestr[a]))
    raw <- readRDS(rawfpath)
    FC <- readRDS(fcfpath)
    #### generate QC plots ####
    ###########################
    if(subsetRaw){
      if(writeRaw[a]){
        dt <- as.matrix(raw[,!colnames(raw) %in% c("ID", "ENTREZID", "SYMBOL", "GENENAME"), with = FALSE])
        # box-and-whisker plot
        plot_df <- suppressWarnings(reshape2::melt(dt))
        png(filename = paste(RawQCPath, "/BoxAndWhisker", "_", Technology[a], "_", paste(paste(DS[a], gsub("-", "_", gsub("_.+", "", nam[1])), sep = "_")), ".png", sep = ""), width = 2000, height = 1500, units = "px", res = 300)
        print(ggplot(data=plot_df, aes(x=Var2, y=value))+
          geom_boxplot()+
          theme(axis.text.x = element_text(colour = 'black',angle = 90))+
          ggtitle(paste(paste(paste(DS[a], gsub("-", "_", gsub("_.+", "", nam[1])), sep = "_")), "_", "box-and-whisker", "_", Technology[a], sep ="")))
        dev.off()
      }
      } else {
        if(writeRaw[a]){
          dt <- as.matrix(raw[,!colnames(raw) %in% c("ID", "ENTREZID", "SYMBOL", "GENENAME"), with = FALSE])
          # box-and-whisker plot
          plot_df <- suppressWarnings(reshape2::melt(dt))
          png(filename = paste(RawQCPath, "/BoxAndWhisker", "_", Technology[a], "_", DS[a], ".png", sep = ""), width = 2000, height = 1500, units = "px", res = 300)
          print(ggplot(data=plot_df, aes(x=Var2, y=value))+
            geom_boxplot()+
            theme(axis.text.x = element_text(colour = 'black',angle = 90))+
            ggtitle(paste(DS[a], "_", "box-and-whisker", "_", Technology[a], sep ="")))
          dev.off()
        }
    }
    #### Set up fold change information ####
    annotation_columns <- intersect(
      c("ID", "SYMBOL", "ENTREZID", "GENENAME"), colnames(raw)
    )
    raw2 <- raw[, setdiff(colnames(raw), annotation_columns), with = FALSE]
    if (subsetRaw) {
      labels <- strsplit(
        toupper(gsub("[[:space:],;]", "", as.character(gsm[a]))),
        "", fixed = TRUE
      )[[1]]
      selected_labels <- labels[labels != "X"]
      if (length(selected_labels) != ncol(raw2) ||
          !identical(sort(unique(selected_labels)), c("0", "1"))) {
        stop(
          "Stored subset and ComparisonVector do not agree for ", DS[a], ".",
          call. = FALSE
        )
      }
    } else {
      comparison <- .parse_comparison_code(gsm[a], colnames(raw2), DS[a])
      selected_labels <- comparison$labels[comparison$selected]
      raw2 <- raw2[, which(comparison$selected), with = FALSE]
    }
    manual_logfc <- NULL
    manual_symbols <- NULL
    if(Technology[a] == "RNAseq"){
      #### perform normalization for RNAseq samples ####
      tbl <- as.matrix(raw2)
      storage.mode(tbl) <- "numeric"
      # Row positions are carried through the count filter so that the
      # direction-only fold change stays aligned with the stored annotations.
      rownames(tbl) <- as.character(seq_len(nrow(tbl)))
      gs <- factor(selected_labels)
      groups <- make.names(c("Ctrl", "Tx"))
      levels(gs) <- groups
      sample_info <- data.frame(Group = gs, row.names = colnames(tbl))
      keep <- rowSums( tbl >= 10 ) >= min(table(gs))
      tbl <- tbl[keep, , drop = FALSE]
      if(!nrow(tbl)){
        stop("No features passed the count filter for ", DS[a], ".", call. = FALSE)
      }
      ds <- DESeqDataSetFromMatrix(countData=tbl, colData=sample_info, design= ~Group)
      #### Obtain normalized count values ####
      # Only size factors are needed here, so the full Wald fit is not run.
      ds <- DESeq2::estimateSizeFactors(ds, type = "poscounts")
      normalized <- DESeq2::counts(ds, normalized = TRUE)
      #### Run QC analysis on normalized data ####
      ############################################
          png(filename = paste(RawQCPath, "/BoxAndWhisker", "_", Technology[a], "_", paste(paste(DS[a], gsub("-", "_", gsub("_.+", "", nam[1])), sep = "_")), "_Normalized", ".png", sep = ""), width = 2000, height = 1500, units = "px", res = 300)
          par(mar=c(12,4,2,1))
          boxplot(normalized, boxwex=0.7, notch=TRUE, main=paste(DS[a], "_", "box-and-whisker", "_Normalized", "_", Technology[a], sep =""), outline=FALSE, las=2)
          dev.off()
      #### calculate an independent direction-only fold change ####
      onetx <- rowMeans(normalized[, selected_labels == "1", drop = FALSE])
      zeroCtrl <- rowMeans(normalized[, selected_labels == "0", drop = FALSE])
      manual_logfc <- log2(onetx + 0.5) - log2(zeroCtrl + 0.5)
      manual_symbols <- raw$SYMBOL[as.integer(rownames(normalized))]
      }
    if(Technology[a] == "Array"){
      #### put array values on the log2 scale when needed ####
      ex <- as.matrix(raw2)
      storage.mode(ex) <- "numeric"
      qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=TRUE))
      LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] < 2)
      if(LogC){
        ex[ex <= 0] <- NA_real_
        ex <- log2(ex)
      }
      onetx <- rowMeans(ex[, selected_labels == "1", drop = FALSE], na.rm = TRUE)
      zeroCtrl <- rowMeans(ex[, selected_labels == "0", drop = FALSE], na.rm = TRUE)
      manual_logfc <- onetx - zeroCtrl
      manual_symbols <- raw$SYMBOL
    }
    if(is.null(manual_logfc)){
      stop(
        "Unsupported Technology `", Technology[a], "` for ", DS[a],
        ". Use `Array` or `RNAseq`.",
        call. = FALSE
      )
    }
    manFC <- data.table(SYMBOL = manual_symbols, log = manual_logfc)
    mer <- merge(FC, manFC, by = "SYMBOL")
    FCcol <- nam[1]
    if (!FCcol %in% colnames(mer)) {
      stop("Fold-change column `", FCcol, "` was not found for ", DS[a], ".", call. = FALSE)
    }
    mer <- mer[!is.na(mer[[FCcol]]),]
    mer <- mer[!is.na(mer$log),]
    p <- ggplot(mer, aes(x=mer[[FCcol]], y=log)) + geom_point() + xlab(FCcol)
    #### write graph to file ####
    graphpath <- file.path(GraphPath, paste(a, "_", paste(DS[a], gsub("-", "_", gsub("_.+", "", nam[1])), sep = "_"), ".png", sep = "")) # ".xls"
    png(filename = graphpath, width = 2000, height = 1500, units = "px", res = 300); print(p); dev.off()
    plotList[[length(plotList) + 1L]] <- p
    p <- NULL
    plotnames[length(plotList)] <- tools::file_path_sans_ext(basename(fcfpath))
    #### fit a linear model ####
    linear <- lm(as.formula(paste( paste("`", FCcol, "`", sep = ""), "~", "log")), mer)
    LMIntercept <- unname(linear$coefficients[1])
    LMSlope <- unname(linear$coefficients[2])
    #### quantify differences ####
    corelation <- stats::cor(mer[[FCcol]], mer$log)
    temp <- data.table(dataset = paste(paste(DS[a], gsub("-", "_", gsub("_.+", "", nam[1])), sep = "_"), sep = ""),
                       correct = sum(nrow(mer[mer[[FCcol]] > 0 & log > 0,]), nrow(mer[mer[[FCcol]] < 0 & log < 0,])),
                       incorrect = sum(
                         nrow(mer[mer[[FCcol]] > 0 & log < 0,]),
                         nrow(mer[mer[[FCcol]] < 0 & log > 0,])
                       ),
                       correlation= corelation,
                       LMEstimate = LMIntercept,
                       LMSlope = LMSlope,
                       iteration = a)
    checkFile <- rbind(checkFile, temp )
    print(paste(a, "of", length(DS), "Completed")) }
  names(plotList) <- plotnames
  return(list(checkFile=checkFile, plotList = plotList)) }

#' Harmonize externally analyzed expression data
#'
#' Converts externally generated differential-expression and raw-count tables
#' to the gene identifiers and file naming conventions expected by the app.
#'
#' @param Fpath A path containing externally analyzed text files.
#' @param OutPath A path where harmonized RDS files are written.
#'
#' @return Invisibly, `NULL`. Harmonized files are written to `OutPath`.
#' @family transcriptomics functions
#' @export
ExternalDataHarmonize <- function(Fpath, OutPath){
  files <- list.files(Fpath)
    #### DE files ####
    ##################
    DEfiles <- files[grepl("Rbound", files)]
    for(a in seq_along(DEfiles)){
      temp <- fread(file.path(Fpath, DEfiles[a]))
      setnames(temp, c("pvalue", "FDR"), c("Pvalue", "AdjPValue"), skip_absent = TRUE)
      comps <- unique(temp$Subset_Comparison)
      for(b in seq_along(comps)){
        temp2 <- temp[Subset_Comparison == comps[b],]
        setnames(temp2, "GeneID", "ENSEMBL")
        if(sum(grepl("ENSMUSG", temp2$ENSEMBL)) > 0 ){ # Mouse
          symbols <- temp2$ENSEMBL
          go <- as.data.table( suppressMessages(AnnotationDbi::select(org.Mm.eg.db, symbols, c("GENENAME"), "ENSEMBL") ))
          temp2 <- merge(go, temp2, by = "ENSEMBL", all.y=TRUE)
          setnames(temp2, "ENSEMBL", "ID") }
        if(sum(grepl("ENSG", temp2$ENSEMBL)) > 0 ){ # Human
          symbols <- temp2$ENSEMBL
          go <- as.data.table( suppressMessages(AnnotationDbi::select(org.Hs.eg.db, symbols, c("GENENAME"), "ENSEMBL") ))
          temp2 <- merge(go, temp2, by = "ENSEMBL", all.y=TRUE)
          setnames(temp2, "ENSEMBL", "ID") }
        if(sum(grepl("ENSRNOG", temp2$ENSEMBL)) > 0 ){ # Human
          symbols <- temp2$ENSEMBL
          go <- as.data.table( suppressMessages(AnnotationDbi::select(org.Rn.eg.db, symbols, c("GENENAME"), "ENSEMBL") ))
          temp2 <- merge(go, temp2, by = "ENSEMBL", all.y=TRUE)
          setnames(temp2, "ENSEMBL", "ID") }
        #### set up file name ####
        Fname <- unique(temp2$Subset_Comparison)
        Fname <- paste(paste(gsub("_.+", "", DEfiles[a]), gsub("-", "_", gsub("_", "", Fname)), sep = "_"), ".txt", sep = "") # ".xls"
        #### save file ####
        temp2 <- temp2[,c("ID", "ENTREZID", "SYMBOL", "GENENAME", "logFC", "Pvalue", "AdjPValue"), with = FALSE]
        # fwrite(temp2, file.path(OutPath, Fname), row.names = FALSE, quote = FALSE, sep = "\t")
        saveRDS(temp2, file = gsub(".txt", ".rds", file.path(OutPath, Fname)) )
      }
      print(paste(a, "of", length(DEfiles), "DE files completed"))
    }
    #### FPKM Files ####
    ####################
    Rawfiles <- files[grepl("RPKMCodingGenes", files)]
    if(length(Rawfiles) > 0){
    for(a in seq_along(Rawfiles)){
      temp <- fread(file.path(Fpath, Rawfiles[a]))
      setnames(temp, "GeneID", "ENSEMBL")
      if(sum(grepl("ENSMUSG", temp$ENSEMBL)) > 0 ){ # Mouse
        symbols <- temp$ENSEMBL
        go <- as.data.table( suppressMessages(AnnotationDbi::select(org.Mm.eg.db, symbols, c("ENTREZID", "SYMBOL", "GENENAME"), "ENSEMBL") ))
        temp <- merge(go, temp, by = "ENSEMBL", all.y=TRUE)
        setnames(temp, "ENSEMBL", "ID") }
      if(sum(grepl("ENSG", temp$ENSEMBL)) > 0 ){ # Human
        symbols <- temp$ENSEMBL
        go <- as.data.table( suppressMessages(AnnotationDbi::select(org.Hs.eg.db, symbols, c("ENTREZID", "SYMBOL", "GENENAME"), "ENSEMBL") ))
        temp <- merge(go, temp, by = "ENSEMBL", all.y=TRUE)
        setnames(temp, "ENSEMBL", "ID") }
      if(sum(grepl("ENSRNOG", temp$ENSEMBL)) > 0 ){ # Rat
        symbols <- temp$ENSEMBL
        go <- as.data.table( suppressMessages(AnnotationDbi::select(org.Rn.eg.db, symbols, c("ENTREZID", "SYMBOL", "GENENAME"), "ENSEMBL") ))
        temp <- merge(go, temp, by = "ENSEMBL", all.y=TRUE)
        setnames(temp, "ENSEMBL", "ID") }
      Fname <- paste(gsub(".txt", "", Rawfiles[a]), "_RPKMRaw", ".rds", sep = "")
      saveRDS(temp, file = gsub(".txt", ".rds", file.path(OutPath, Fname)) )
      print(paste(a, "of", length(Rawfiles), "RPKM Raw files completed")) }
    } else {
      print("No files containing RPKM values were found")
    }
    #### Count files ####
    #####################
    Rawfiles <- files[grepl("CountCodingGenes", files)]
    if(length(Rawfiles) > 0){
    for(a in seq_along(Rawfiles)){
      temp <- fread(file.path(Fpath, Rawfiles[a]))
      setnames(temp, "GeneID", "ENSEMBL")
      if(sum(grepl("ENSMUSG", temp$ENSEMBL)) > 0 ){ # Mouse
        symbols <- temp$ENSEMBL
        go <- as.data.table( suppressMessages(AnnotationDbi::select(org.Mm.eg.db, symbols, c("ENTREZID", "SYMBOL", "GENENAME"), "ENSEMBL") ))
        temp <- merge(go, temp, by = "ENSEMBL", all.y=TRUE)
        setnames(temp, "ENSEMBL", "ID") }
      if(sum(grepl("ENSG", temp$ENSEMBL)) > 0 ){ # Human
        symbols <- temp$ENSEMBL
        go <- as.data.table( suppressMessages(AnnotationDbi::select(org.Hs.eg.db, symbols, c("ENTREZID", "SYMBOL", "GENENAME"), "ENSEMBL") ))
        temp <- merge(go, temp, by = "ENSEMBL", all.y=TRUE)
        setnames(temp, "ENSEMBL", "ID") }
      if(sum(grepl("ENSRNOG", temp$ENSEMBL)) > 0 ){ # Rat
        symbols <- temp$ENSEMBL
        go <- as.data.table( suppressMessages(AnnotationDbi::select(org.Rn.eg.db, symbols, c("ENTREZID", "SYMBOL", "GENENAME"), "ENSEMBL") ))
        temp <- merge(go, temp, by = "ENSEMBL", all.y=TRUE)
        setnames(temp, "ENSEMBL", "ID") }
      Fname <- paste(gsub(".txt", "", Rawfiles[a]), "_CountRaw", ".rds", sep = "")
      saveRDS(temp, file = gsub(".txt", ".rds", file.path(OutPath, Fname)) )
      print(paste(a, "of", length(Rawfiles), "Count Raw files completed")) }
    }  else {
      print("No files containing count values were found")
    }
  invisible(NULL)
}

#' Compile differential-expression data
#'
#' Combines harmonized transcriptomic differential-expression files into a
#' single [SummarizedExperiment::SummarizedExperiment()] with one assay per
#' comparison.
#'
#' @param DEGDatapath A path containing differential-expression RDS files.
#' @param SEPath A file path where the compiled object is saved.
#'
#' @return A `SummarizedExperiment` containing one three-column assay per
#'   comparison and feature annotations keyed by a species-qualified stable
#'   identifier. The database is rebuilt deterministically and saved to
#'   `SEPath`.
#' @family transcriptomics functions
#' @import data.table
#' @importFrom S4Vectors SimpleList
#' @export
DESEGenerate <- function(DEGDatapath, SEPath){
  .compile_de_experiments(DEGDatapath, SEPath)
}
#' Compile proteomic differential-expression data
#'
#' Combines processed proteomic differential-expression files into a single
#' [SummarizedExperiment::SummarizedExperiment()] with one assay per
#' comparison. When `SEPath` already holds a database, the datasets it does not
#' yet contain are appended to it and the existing protein annotations are
#' reused, so that every assay stays aligned to the same features.
#'
#' @param ProtDatapath A path containing proteomic differential-expression RDS
#'   files.
#' @param SEPath A file path where the compiled object is saved.
#'
#' @return A `SummarizedExperiment` containing the compiled assays and protein
#'   annotations. Each assay holds `logFC`, `AdjPValue`, `Pvalue`, and
#'   `BHCorrection` columns. The same object is saved to `SEPath`.
#' @family proteomics functions
#' @import data.table
#' @importFrom S4Vectors SimpleList
#' @export
DESEProtGenerate <- function(ProtDatapath, SEPath){
  files <- list.files(ProtDatapath, pattern = "[.]rds$", full.names = TRUE)
  if (!length(files)) {
    stop("No proteomics RDS files were found in `ProtDatapath`.", call. = FALSE)
  }
  existing <- if (file.exists(SEPath)) readRDS(SEPath) else NULL
  if (is.null(existing)) {
    message("No database found at ", SEPath, ". Compiling from scratch.")
    row_data <- .proteomic_row_data(.proteomic_symbol_annotation())
    existing_assays <- list()
  } else {
    message("Updating the existing proteomics database.")
    row_data <- .proteomic_row_data(
      SummarizedExperiment::rowData(existing, use.names = FALSE)
    )
    existing_assays <- as.list(SummarizedExperiment::assays(existing))
    already_present <- tools::file_path_sans_ext(basename(files)) %in%
      names(existing_assays)
    if (any(already_present)) {
      message(
        "Ignoring datasets already present in the database: ",
        paste(basename(files)[already_present], collapse = ", ")
      )
      files <- files[!already_present]
    }
  }
  feature_symbols <- rownames(row_data)
  new_assays <- list()
  for (i in seq_along(files)) {
    new_assays[[i]] <- .proteomic_assay(files[i], feature_symbols)
    names(new_assays)[i] <- tools::file_path_sans_ext(basename(files[i]))
    message("Completed ", i, " of ", length(files), ": ", basename(files[i]))
  }
  if (!length(files)) {
    message("The database is already up to date.")
  }
  SE <- SummarizedExperiment(
    assays = SimpleList(c(existing_assays, new_assays)),
    rowData = row_data
  )
  directory <- dirname(SEPath)
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  saveRDS(SE, file = SEPath)
  SE
}

#' Configure copied application files
#'
#' Replaces the database-path placeholder in the copied Shiny server and UI
#' scripts, including the workflow driver. Run this after [makeDirectory()].
#'
#' @param homedir A path to the database directory.
#'
#' @return Invisibly, `NULL`. The copied scripts are updated in place.
#' @family setup functions
#' @export
AppSetup <- function(homedir){
  script_paths <- file.path(
    homedir,
    "Scripts",
    c("server.R", "UI.R", "WorkflowScript.R")
  )
  missing_scripts <- script_paths[!file.exists(script_paths)]
  if (length(missing_scripts)) {
    stop("Missing application scripts: ", paste(missing_scripts, collapse = ", "),
         call. = FALSE)
  }
  database_path <- normalizePath(homedir, winslash = "/", mustWork = TRUE)
  for (script_path in script_paths) {
    contents <- readLines(script_path, warn = FALSE)
    contents <- gsub("<path to data base>", database_path, contents,
                     fixed = TRUE)
    writeLines(contents, script_path)
  }
  invisible(NULL)
  }

#' Download ProteomeXchange datasets
#'
#' Downloads files associated with one or more ProteomeXchange accessions and
#' moves them into accession-specific subdirectories.
#'
#' @param path A path used for downloads and the BiocFileCache cache.
#' @param DS A character vector of ProteomeXchange PXD accessions.
#' @param DownloadRawData A logical value indicating whether vendor raw-data
#'   files should also be downloaded.
#'
#' @return Invisibly, `NULL`. Downloaded files are written beneath `path`.
#' @family proteomics functions
#' @importFrom rpx PXDataset
#' @importFrom rpx pxfiles
#' @importFrom rpx pxget
#' @importFrom BiocFileCache BiocFileCache
#' @import data.table
#' @export
ProteomicsDataDownload <- function(path, DS, DownloadRawData = FALSE){
  bfc <- BiocFileCache(path, ask = FALSE)
  for(i in seq_along(DS)){
    px1 <- tryCatch(
      suppressWarnings(PXDataset(DS[i])),
      error = function(e) {
        message(DS[i], " could not be downloaded. Check the PXD accession.")
        NULL
      }
    )
    if(!is.null(px1)){
      if(sum(list.files(file.path(path)) %in% DS[i]) ==0){ dir.create(file.path(path, DS[i])) }
      files <- suppressMessages(pxfiles(px1))
      if(!DownloadRawData){ files <- files[!grepl(".raw$|raw|.baf|RAW|.xml|.sne|.mgf|.wiff|.group", files, ignore.case = TRUE)] }
      print(files)
      mztab <- pxget(px1, files, cache = bfc)
      #### move files into the appropriate directory ####
      print(paste("moving files to", DS[i], "directory"))
      file.copy(from = mztab, to   = file.path(path, DS[i], files))
      print("removing files from old directory")
      file.remove(mztab)
    }
    print(paste("Completed", i, "of", length(DS)))
  }
  invisible(NULL)
}

#' Format supported MaxQuant datasets
#'
#' Applies dataset-specific parsing rules to the supported MaxQuant exports in
#' `path`. Proteomics repositories do not use a single tabular convention, so
#' new datasets may require an additional formatting branch.
#'
#' @param path A path containing downloaded proteomics datasets.
#'
#' @return A named list of harmonized proteomics data tables.
#' @family proteomics functions
#' @import data.table
#' @export
FormatMaxQuant <- function(path){
  files <- gsub(".+/", "", list.dirs(path, recursive = FALSE) )
  MZList <- list()
  for(i in seq_along(files)){
    if(files[i] == "PXD005847"){
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      nam <- fread(fpath[grepl("experimentalDesignTemplate.txt", fpath)])
      data <- fread(fpath[grepl("proteinGroups.txt", fpath)])
      data <- data[,c("Majority protein IDs","Peptide counts (all)", colnames(data)[grep("Intensity", colnames(data), ignore.case = TRUE)],
                      "Protein names", "Gene names"), with = FALSE]
      data$Sample <- "HDL"
      data <- cbind(data[,!grepl("intensity", colnames(data), ignore.case = TRUE), with = FALSE],
                    data[,grepl("LFQ", colnames(data), ignore.case = TRUE), with = FALSE])
      data <- data[,c("Sample", "Majority protein IDs", "Gene names", "Protein names", "Peptide counts (all)",
                      colnames(data)[grepl("LFQ", colnames(data))]), with = FALSE]
      setnames(data, colnames(data)[1:5], c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))
      data$Numberofpeptides <- sapply(data$Numberofpeptides, function(x){ sum(as.numeric(strsplit(x, ";")[[1]])) })
      Cnames <- gsub("_.+", "", colnames(data) )
      Cnames <- makeunique::make_unique(Cnames, "_", wrap_in_brackets = FALSE)
      setnames(data, colnames(data), Cnames)
      MZList[[i]] <- data
      names(MZList)[i] <- files[i]
    }
    if(files[i] == "PXD008934"){
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      data <- fread(fpath[grepl("Prosser_HumanHearts_MassSpec.txt", fpath)])#fpath[2])
      data$Sample <- "Heart"
      data <- cbind(data[,!grepl("intensity", colnames(data), ignore.case = TRUE), with = FALSE],
                    data[,grepl("LFQ", colnames(data), ignore.case = TRUE), with = FALSE])
      data <- data[,c("Sample", "Majority protein IDs",  "Gene names", "Protein names","Peptides",
                      colnames(data)[grepl("LFQ", colnames(data))]), with = FALSE]
      setnames(data, colnames(data)[1:5], c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))
      meta <- fread(file.path(path, files[i], "sample_details.txt"), skip = 2)
      meta <- meta[!meta$`File #` == "",]
      meta <- meta[1:34,1:7, with = FALSE]
      Cnames <- gsub(" S", " ", colnames(data))
      Cnames <- mgsub(Cnames, meta$`File #`, meta$Etiology)
      Cnames <- makeunique::make_unique(Cnames, "_", wrap_in_brackets = FALSE)
      setnames(data, colnames(data), Cnames)
      MZList[[i]] <- data
      names(MZList)[i] <- files[i]
    }
    if(files[i] == "PXD011283"){
      # fpath <- file.path(path, files[i])
      # fpath <- file.path(fpath, list.files(fpath))
      # # unzip(fpath[1], exdir = file.path(path, files[i]))
      # fpath <- file.path(path, files[i])
      # fpath <- file.path(fpath, list.files(fpath))
      # fpath <- file.path(fpath[1], list.files(fpath[1]))
      # # meta <- fread(fpath[grepl("experimentalDesignTemplate.txt", fpath)])
      # data <- fread(fpath[grepl("peptides.txt", fpath)])
      # data$Sample <- "MV_urine"
      # data <- cbind(data[,!grepl("Experiment|Identification ", colnames(data), ignore.case = TRUE), with = FALSE],
      #               data[,grepl("Intensity MV", colnames(data), ignore.case = TRUE), with = FALSE])
      # data <- data[,c("Sample", "Proteins",  "Gene names", "Protein names","Length",
      #                 colnames(data)[grepl("Intensity MV", colnames(data))]), with = FALSE]
      # setnames(data, colnames(data)[1:5], c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))
      # Cnames <- colnames(data)
      # Cnames <- gsub("MV_urine_", "",  Cnames)
      # Cnames <- gsub("_.+", "",  Cnames)
      # colnames(data) <- makeunique::make_unique(Cnames, "_", wrap_in_brackets = FALSE)
      # MZList[[i]] <- data
      # names(MZList)[i] <- files[i]
    }
    if(files[i] == "PXD011839"){
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      # nam <- fread(fpath[grepl("sdrf.tsv", fpath)])
      data <- fread(fpath[grepl("proteinGroups.txt", fpath)])
      data <- data[,c("Majority protein IDs","Peptide counts (all)", colnames(data)[grep("Intensity", colnames(data), ignore.case = TRUE)],
                      "Protein names", "Gene names"), with = FALSE]
      data$Sample <- "Plasma"
      data <- cbind(data[,!grepl("intensity", colnames(data), ignore.case = TRUE), with = FALSE],
                    data[,grepl("LFQ", colnames(data), ignore.case = TRUE), with = FALSE])
      data <- data[,c("Sample", "Majority protein IDs", "Gene names", "Protein names", "Peptide counts (all)",
                      colnames(data)[grepl("LFQ", colnames(data))]), with = FALSE]
      setnames(data, colnames(data)[1:5], c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))
      data$Numberofpeptides <- sapply(data$Numberofpeptides, function(x){ sum(as.numeric(strsplit(x, ";")[[1]])) })
      data <- data[,!grepl("Liver_Rep", colnames(data)), with = FALSE]
      data <- data[,!grepl("_Top14Top6Depl_", colnames(data)), with = FALSE]
      data <- data[,!grepl("_Top6Top14Depl_", colnames(data)), with = FALSE]
      Cnames <- gsub("20170405.+_Human_", "", colnames(data))
      Cnames <- gsub("_F.+", "", Cnames)
      colnames(data) <- makeunique::make_unique(Cnames, "_", wrap_in_brackets = FALSE)
      data$ProteinID <- sapply(data$ProteinID, function(x){ strsplit(x, ";")[[1]][1] })
      data$GeneSymbol <- sapply(data$GeneSymbol, function(x){ strsplit(x, ";")[[1]][1] })
      data$Description <- sapply(data$Description, function(x){ strsplit(x, ";")[[1]][1] })
      MZList[[i]] <- data
      names(MZList)[i] <- files[i]
    }
    if(files[i] == "PXD022545"){
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      # system(paste("7za x", fpath[2]))
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      fpath <- file.path(fpath[4], list.files(fpath[4]))
      data <- fread(fpath[grepl("proteinGroups.txt", fpath)])
      data <- data[,c("Majority protein IDs","Peptide counts (all)", colnames(data)[grep("Intensity", colnames(data), ignore.case = TRUE)], "Fasta headers"), with = FALSE]
      data$GeneSymbol <- sapply(data$`Fasta headers`, function(x){test <- gsub("sp.+\\|| OS.+|HUMAN ", "", x); spl <- str_split(test, "_")[[1]]; spl[1] })
      data$Description <- sapply(data$`Fasta headers`, function(x){test <- gsub("sp.+\\|| OS.+|HUMAN ", "", x); spl <- str_split(test, "_")[[1]]; spl[2] })
      data$`Fasta headers` <- NULL
      data$Sample <- "retinum"
      data <- cbind(data[,!grepl("intensity", colnames(data), ignore.case = TRUE), with = FALSE],
                    data[,grepl("LFQ", colnames(data), ignore.case = TRUE), with = FALSE])
      data <- data[,c("Sample", "Majority protein IDs", "GeneSymbol", "Description", "Peptide counts (all)",
                      colnames(data)[grepl("LFQ", colnames(data))]), with = FALSE]
      setnames(data, colnames(data)[1:5], c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))
      data$Numberofpeptides <- sapply(data$Numberofpeptides, function(x){ sum(as.numeric(strsplit(x, ";")[[1]])) })
      data$ProteinID <- sapply(data$ProteinID, function(x){ strsplit(x, ";")[[1]][1] })
      data$GeneSymbol <- sapply(data$GeneSymbol, function(x){ strsplit(x, ";")[[1]][1] })
      data$Description <- sapply(data$Description, function(x){ strsplit(x, ";")[[1]][1] })
      colnames(data) <- gsub("_0", "_", colnames(data))
      MZList[[i]] <- data
      names(MZList)[i] <- files[i]
    }
    if(files[i] == "PXD029200"){
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      meta <- as.data.table(read_excel(fpath[grepl("Pride_explanatory_sample_file.xlsx", fpath)]))
      data <- suppressWarnings(as.data.table(read_xlsx(fpath[grepl("Data_from_maxquant_used_for_filtering.xlsx", fpath)])))
      data$Sample <- "Hepatocyte"
      data <- cbind(data[,!grepl("intensity", colnames(data), ignore.case = TRUE), with = FALSE],
                    data[,grepl("LFQ", colnames(data), ignore.case = TRUE), with = FALSE])
      data <- data[,c("Sample", "Protein IDs", "Gene names", "Protein names", "Peptide counts (all)",
                      colnames(data)[grepl("LFQ", colnames(data))]), with = FALSE]
      setnames(data, colnames(data)[1:5], c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))
      metasub <- meta[,c("Name in data analysis", "Day", "Treatment"), with = FALSE]
      metasub$Treatment <- gsub(" ", "", gsub("Low glucose, ", "", metasub$Treatment))
      metasub$Day <- paste("Day", metasub$Day, sep = "")
      metasub$SampleTx <- paste(metasub$Treatment, metasub$Day, sep = "_")
      metasub$colname <- paste("LFQ intensity", metasub$`Name in data analysis`)
      metasub$colname[1:3] <- c("LFQ intensity 011", "LFQ intensity 012", "LFQ intensity 013")
      setnames(data, metasub$colname, make.unique(metasub$SampleTx, "_"))
      data$Numberofpeptides <- sapply(data$Numberofpeptides, function(x){ sum(as.numeric(strsplit(x, ";")[[1]])) })
      Cnames <- colnames(data)
      colnames(data) <- gsub("_Day", "Day",  Cnames)
      colnames(data) <- makeunique::make_unique(gsub("_[0-9]+", "", colnames(data)), "_", wrap_in_brackets = FALSE)
      colnames(data)[6:length(colnames(data))] <- paste("LFQ intensity", colnames(data)[6:length(colnames(data))])
      MZList[[i]] <- data
      names(MZList)[i] <- files[i]
    }
    if(files[i] == "PXD030764"){
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      data <- fread(fpath[grepl("proteinGroups.txt", fpath)])
      data$Sample <- "Hepatocyte"
      data <- cbind(data[,!grepl("intensity", colnames(data), ignore.case = TRUE), with = FALSE],
                    data[,grepl("LFQ", colnames(data), ignore.case = TRUE), with = FALSE])
      data <- data[,c("Sample", "Protein IDs", "Gene names", "Protein names", "Peptide counts (all)",
                      colnames(data)[grepl("LFQ", colnames(data))]), with = FALSE]
      setnames(data, colnames(data)[1:5], c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))
      data$Numberofpeptides <- sapply(data$Numberofpeptides, function(x){ sum(as.numeric(strsplit(x, ";")[[1]])) })
      colnames(data)[6:length(colnames(data))] <- mgsub(colnames(data)[6:length(colnames(data))], c("M", "P", "O"), c("myristic", "palmitic", "oleic"))
      Cnames <- makeunique::make_unique(gsub("_[0-9]+", "", colnames(data)), "_", wrap_in_brackets = FALSE)
      Cnames <- gsub("LFQ intensity palmiticoleic12", "LFQ intensity palmiticoleic12_1", Cnames)
      Cnames <- gsub("LFQ intensity palmiticmyristic11", "LFQ intensity palmiticmyristic11_1", Cnames)
      setnames(data, colnames(data), Cnames)
      data <- data[!GeneSymbol == "",]
      data$ProteinID <- sapply(data$ProteinID, function(x){ strsplit(x, ";")[[1]][1] })
      data$GeneSymbol <- sapply(data$GeneSymbol, function(x){ strsplit(x, ";")[[1]][1] })
      data$Description <- sapply(data$Description, function(x){ strsplit(x, ";")[[1]][1] })
      MZList[[i]] <- data
      names(MZList)[i] <- files[i]
    }
    if(files[i] == "PXD011563"){
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      data <- fread(fpath[grepl("proteinGroups.txt", fpath)])
      data <- data[,c("Majority protein IDs","Peptide counts (all)", colnames(data)[grep("Intensity", colnames(data), ignore.case = TRUE)],
                      "Protein names", "Gene names"), with = FALSE]
      data$Sample <- "Microparticle"
      data <- cbind(data[,!grepl("intensity", colnames(data), ignore.case = TRUE), with = FALSE],
                    data[,grepl("LFQ", colnames(data), ignore.case = TRUE), with = FALSE])
      data <- data[,c("Sample", "Majority protein IDs", "Gene names", "Protein names", "Peptide counts (all)",
                      colnames(data)[grepl("LFQ", colnames(data))]), with = FALSE]
      setnames(data, colnames(data)[1:5], c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))
      data$Numberofpeptides <- sapply(data$Numberofpeptides, function(x){ sum(as.numeric(strsplit(x, ";")[[1]])) })
      Cnames <- gsub("_.+", "", colnames(data) )
      Cnames <- makeunique::make_unique(Cnames, "_", wrap_in_brackets = FALSE)
      setnames(data, colnames(data), Cnames)
      MZList[[i]] <- data
      names(MZList)[i] <- files[i]
    }
    if(files[i] == "PXD024734"){
      fpath <- file.path(path, files[i])
      fpath <- file.path(fpath, list.files(fpath))
      # nam <- fread(fpath[grepl("mzRange.txt", fpath)])
      data <- fread(fpath[grepl("proteinGroups.txt", fpath)])
      # setnames(data, c("Fasta headers"), c("Protein names"))
      # data$`Gene names` <- data$`Protein names`
      data <- data[,c("Majority protein IDs","Peptide counts (all)", colnames(data)[grep("Intensity", colnames(data), ignore.case = TRUE)],
                      "Protein names", "Gene names"), with = FALSE]
      data$Sample <- "Adipose"
      data <- cbind(data[,!grepl("intensity", colnames(data), ignore.case = TRUE), with = FALSE],
                    data[,grepl("LFQ", colnames(data), ignore.case = TRUE), with = FALSE])
      data <- data[,c("Sample", "Majority protein IDs", "Gene names", "Protein names", "Peptide counts (all)",
                      colnames(data)[grepl("LFQ", colnames(data))]), with = FALSE]
      setnames(data, colnames(data)[1:5], c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))
      data$Numberofpeptides <- sapply(data$Numberofpeptides, function(x){ sum(as.numeric(strsplit(x, ";")[[1]])) })
      Cnames <- gsub("[0-9]+", "", colnames(data) )
      Cnames <- makeunique::make_unique(Cnames, "_", wrap_in_brackets = FALSE)
      setnames(data, colnames(data), Cnames)
      MZList[[i]] <- data
      names(MZList)[i] <- files[i]
    }
  }
  MZList <- MZList[!names(MZList) == ""]
  MZList <- MZList[!is.na(names(MZList))]
  return(MZList)
}

#' Save formatted proteomics data
#'
#' @param MZList A named list of formatted proteomics data tables.
#' @param path A path where tab-separated `.xls` files are written.
#'
#' @return Invisibly, `NULL`. One file is written for each element of `MZList`.
#' @family proteomics functions
#' @import data.table
#' @export
proteomicMZSave <- function(MZList, path){
  for(i in seq_along(MZList)){ fwrite(MZList[[i]], file.path(path, paste(names(MZList)[i], ".xls", sep = "")), row.names=FALSE, quote=FALSE, sep="\t") }
  invisible(NULL)
}

#' Build a proteomics design matrix from file names
#'
#' Derives dataset, sample label, condition, and replicate fields from the file
#' names in a directory of formatted proteomics tables.
#'
#' @param Fpath A path containing formatted proteomics files.
#'
#' @return A data table describing the experimental design.
#' @family proteomics functions
#' @import data.table
#' @export
DesignMatrixFromNames <- function(Fpath){
  files <- list.files(Fpath)
  DesignDT <- data.table()
  for(i in seq_along(files)){
    Cnames <- colnames(fread(file.path(Fpath, files[i])))
    Cnames <- Cnames[!(Cnames %in% c("Sample", "ProteinID", "GeneSymbol", "Description", "Numberofpeptides"))]
    DT <- data.table(dataset = gsub(".xls", "", files[i]),
                     label = Cnames,
                     condition = gsub("_[0-9]+$", "", gsub(".+ ", "", Cnames)),
                     replicate = gsub(".+_", "", Cnames) )
    DesignDT <- rbind(DesignDT, DT)
  }
  replicate <- gsub(".+_", "", makeunique::make_unique(DesignDT$condition, "_", wrap_in_brackets = FALSE))
  replicate[grepl("[a-z]", replicate, ignore.case = TRUE)] <- "1"
  DesignDT$replicate <- replicate
  return(DesignDT)}

#' Load proteomics data as SummarizedExperiment objects
#'
#' Reads each formatted proteomics table and combines it with the matching rows
#' in a design table.
#'
#' @param DesignDT A data table containing dataset, sample label, condition, and
#'   replicate columns.
#' @param Fpath A path containing formatted proteomics files.
#'
#' @return A named list of `SummarizedExperiment` objects.
#' @family proteomics functions
#' @import data.table
#' @export
ProtSELoad <- function(DesignDT, Fpath){
  .require_optional_package("DEP", "ProtSELoad()")
  SEList <- list()
  files <- list.files(Fpath)
  for(i in seq_along(files)){
    data <- fread(file.path(Fpath, files[i]))
    data_unique <- DEP::make_unique(data, "GeneSymbol", "ProteinID", delim = ";")
    LFQ_columns <- grep("LFQ.", colnames(data_unique))
    if(length(LFQ_columns) == 0){
      LFQ_columns <- grep("Intensity.", colnames(data_unique), ignore.case = TRUE)
    }
    experimental_design <- DesignDT[DesignDT$dataset == gsub(".xls", "", files[i]),]
    experimental_design$dataset <- NULL
    for(a in seq_along(LFQ_columns)){
    data_unique[[LFQ_columns[a]]] <- suppressWarnings(as.numeric(data_unique[[LFQ_columns[a]]]))
    }
    data_se <- DEP::make_se(as.data.frame(data_unique), LFQ_columns, as.data.frame(experimental_design))
    SEList[[i]] <- data_se
    names(SEList)[i] <- gsub(".xls", "", files[i])
  }
  return(SEList) }

#' Reshape assay data to long format
#'
#' @param tissueSplitList A list of `SummarizedExperiment` objects.
#'
#' @return A data table containing gene symbols, sample values, and parsed cell,
#'   tissue, and diet labels.
#' @family proteomics functions
#' @import data.table
#' @export
RowDataCompile <- function(tissueSplitList){
  TotMel <- data.table()
  for(i in seq_along(tissueSplitList)){
    temp <- as.data.frame(assay(tissueSplitList[[i]]))
    temp$GeneSymbol <- rownames(temp)
    temp <- as.data.table(temp)
    mel <- suppressWarnings(suppressMessages(reshape2::melt(temp)))
    mel$Comb <- paste(mel$variable, gsub("_.+", "", names(tissueSplitList)[i]), sep = ":")
    mel$Cell <- gsub("_.+", "", names(tissueSplitList)[i])
    TotMel <- rbind(TotMel, mel)
  }
  TotMel$CellTissue <- paste(TotMel$Cell, gsub("_.+", "", TotMel$variable), sep = ":")
  TotMel$Tissue <- gsub("_.+", "", TotMel$variable)
  TotMel$Diet <- gsub("Plasma_|Liver_", "", TotMel$variable)
  return(TotMel) }

#' Filter features by missing-value count
#'
#' Retains features whose missing-value count does not exceed `thr` in every
#' experimental condition.
#'
#' @param dataSeList A list of `SummarizedExperiment` objects created with DEP.
#' @param thr A non-negative numeric missing-value threshold per condition.
#'
#' @return A named list of filtered `SummarizedExperiment` objects.
#' @family proteomics functions
#' @importFrom assertthat assert_that
#' @importFrom tidyr gather spread
#' @export
DataFilter <- function(dataSeList, thr = 0){
  dataFiltList <- list()
  for(i in seq_along(dataSeList)){
    se = dataSeList[[i]]
    if (is.integer(thr)){
      thr <- as.numeric(thr)}
    assertthat::assert_that(inherits(se, "SummarizedExperiment"), is.numeric(thr), length(thr) == 1)
    if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
      stop("'name' and/or 'ID' columns are not present in '", deparse(substitute(se)), "'\nRun make_unique() and make_se() to obtain the required columns", call. = FALSE)
    }
    if (any(!c("label", "condition", "replicate") %in% colnames(colData(se)))) {
      stop("'label', 'condition' and/or 'replicate' columns are not present in '", deparse(substitute(se)), "'\nRun make_se() or make_se_parse() to obtain the required columns", call. = FALSE)
    }
    max_repl <- max(colData(se)$replicate)
    if (thr < 0 | thr > max_repl) { stop("invalid filter threshold applied", "\nRun filter_missval() with a threshold ranging from 0 to ", max_repl)
    }
    bin_data <- assay(se)
    idx <- is.na(assay(se))
    bin_data[!idx] <- 1
    bin_data[idx] <- 0
    rownames(bin_data) <- 1:nrow(bin_data)
    keep <- data.frame(
      rowname = rownames(bin_data), bin_data, check.names = FALSE
    ) %>%
      gather(ID, value, -rowname) %>% left_join(., data.frame(colData(se)), by = "ID") %>%
      group_by(rowname, condition) %>% summarize(miss_val = n() - sum(value)) %>%
      filter(miss_val <= thr) %>% spread(condition, miss_val)
    dataFiltList[[i]] <- se[as.numeric(keep$rowname), ]
    }
  names(dataFiltList) <- names(dataSeList)
  return(dataFiltList) }

#' Filter proteins by peptide count
#'
#' @param dataFiltList A list of `SummarizedExperiment` objects with a
#'   `Numberofpeptides` row-data column.
#' @param Npeptides The minimum number of peptides required for a protein.
#'
#' @return A list with `data`, the filtered experiments, and `remaining`, a data
#'   table containing the proportion of proteins retained per dataset.
#' @family proteomics functions
#' @export
NPeptideThreshold <- function(dataFiltList, Npeptides = 2){
  PepCutOff <- dataFiltList
  PercentRemain <- NULL
  for(i in seq_along(PepCutOff)){
    PepCutOff[[i]] <- PepCutOff[[i]][rowData(PepCutOff[[i]])$Numberofpeptides >= Npeptides,]
    PercentRemain[i] <- nrow(assay(PepCutOff[[i]]))/nrow(assay(dataFiltList[[i]]))
  }
  remaining <- data.table(dataset = names(PepCutOff), percentRemaining = PercentRemain)
  return(list(data = PepCutOff, remaining = remaining)) }

#' Apply multiple proteomics normalizations
#'
#' Applies mean, median, variance-stabilizing, cyclic loess, robust linear
#' regression, and robust mean/scale normalizations to facilitate method
#' comparison.
#'
#' @param dataPepCutOff A list of peptide-filtered `SummarizedExperiment`
#'   objects.
#'
#' @return A named list whose elements (`mean`, `median`, `vsn`, `DEPvsn`,
#'   `loess`, `rlr`, and `smad`) each contain a normalized experiment list.
#' @family proteomics functions
#' @export
MultiNormalization <- function(dataPepCutOff){
  .require_optional_package("DEP", "MultiNormalization()")
  .require_optional_package("NormalyzerDE", "MultiNormalization()")
  #### MEAN ####
  MeanNormList <- dataPepCutOff
  for(i in seq_along(MeanNormList)){
    assayTemp <- 2^assay(MeanNormList[[i]])
    norm <- NormalyzerDE::meanNormalization(assayTemp)
    rownames(norm) <- rownames(assayTemp)
    assay(MeanNormList[[i]]) <- norm }
  #### MEDIAN ####
  MedianNormList <- dataPepCutOff
  for(i in seq_along(MedianNormList)){
    assayTemp <- 2^assay(MedianNormList[[i]])
    norm <- NormalyzerDE::medianNormalization(assayTemp)
    rownames(norm) <- rownames(assayTemp)
    assay(MedianNormList[[i]]) <- norm }
  #### VSN ####
  VSNNormList <- dataPepCutOff
  for(i in seq_along(VSNNormList)){
    assayTemp <- 2^assay(VSNNormList[[i]])
    norm <- NormalyzerDE::performVSNNormalization(assayTemp)
    rownames(norm) <- rownames(assayTemp)
    assay(VSNNormList[[i]]) <- norm }
  ### DEP package VSN ####
  DEPVSNNormList <- dataPepCutOff
  for(i in seq_along(DEPVSNNormList)){
    DEPVSNNormList[[i]] <- DEP::normalize_vsn(DEPVSNNormList[[i]]) }
  #### Loess ####
  LoessNormList <- dataPepCutOff
  for(i in seq_along(LoessNormList)){
    assayTemp <- assay(LoessNormList[[i]])
    norm <- NormalyzerDE::performCyclicLoessNormalization(
      assayTemp, noLogTransform = TRUE
    )
    rownames(norm) <- rownames(assayTemp)
    assay(LoessNormList[[i]]) <- norm }
  #### RLR ####
  RLRNormList <- dataPepCutOff
  for(i in seq_along(RLRNormList)){
    assayTemp <- assay(RLRNormList[[i]])
    norm <- NormalyzerDE::performGlobalRLRNormalization(
      assayTemp, noLogTransform = TRUE
    )
    rownames(norm) <- rownames(assayTemp)
    assay(RLRNormList[[i]]) <- norm }
  #### SMAD ####
  SMADNormList <- dataPepCutOff
  for(i in seq_along(SMADNormList)){
    assayTemp <- assay(SMADNormList[[i]])
    norm <- NormalyzerDE::performSMADNormalization(
      assayTemp, noLogTransform = TRUE
    )
    rownames(norm) <- rownames(assayTemp)
    assay(SMADNormList[[i]]) <- norm }
  return(list(mean=MeanNormList, median=MedianNormList, vsn=VSNNormList, DEPvsn=DEPVSNNormList, loess=LoessNormList, rlr=RLRNormList, smad=SMADNormList))
}

#' Plot assay-value densities
#'
#' @param MultiNormalizeList A named list of normalization results, where each
#'   element is a list of `SummarizedExperiment` objects.
#'
#' @return A named list of `ggplot` density plots.
#' @family exploration functions
#' @export
densityPlotFromList <- function(MultiNormalizeList){
  densityPlotList <- list()
  for(i in seq_along(MultiNormalizeList)){
    temp <- RowDataCompile(tissueSplitList=MultiNormalizeList[[i]])
    densityPlotList[[i]] <- ggplot(temp, aes(x=value, color=variable))+
      geom_density()+
      theme(legend.position="none",
            panel.spacing = unit(0.1, "lines"),
            strip.text.x = element_text(size = 8) )+
      facet_wrap(~CellTissue)
  }
  names(densityPlotList) <- names(MultiNormalizeList)
  return(densityPlotList) }

#' Identify datasets containing missing values
#'
#' @param data A list of `SummarizedExperiment` objects.
#'
#' @return An integer vector containing the positions of list elements with at
#'   least one missing assay value. The vector is empty when none are missing.
#' @family proteomics functions
#' @export
DetermineMissing <- function(data){
  missing <- vapply(
    data,
    function(x) anyNA(SummarizedExperiment::assay(x)),
    logical(1)
  )
  which(missing)
}

#' @rdname DetermineMissing
#' @export
DetermineMising <- function(data){
  .Deprecated("DetermineMissing", package = "MultiOmicsDataCompile")
  DetermineMissing(data)
}

#' Impute missing proteomics data
#'
#' @param dataFiltList A list of filtered `SummarizedExperiment` objects.
#' @param type One of `"MinProb"`, `"man"`, or `"knn"`, identifying the
#'   imputation strategy.
#'
#' @return A named list of imputed `SummarizedExperiment` objects.
#' @family proteomics functions
#' @export
DataImpute <- function(dataFiltList, type = "MinProb"){
  type <- match.arg(type, c("MinProb", "man", "knn"))
  ImputeList <- list()
  for(i in seq_along(dataFiltList)){
    if(type == "MinProb"){
      # Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
      ImputeList[[i]] <- impute2(dataFiltList[[i]], fun = "MinProb")}
    if(type == "man"){
      # Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
      ImputeList[[i]] <- impute2(dataFiltList[[i]], fun = "man", shift = 1.8, scale = 0.3) }
    if(type == "knn"){
      # Impute missing data using the k-nearest neighbour approach (for MAR)
      ImputeList[[i]] <- impute2(dataFiltList[[i]], fun = "knn", rowmax = 0.9) }
  }
  names(ImputeList) <- names(dataFiltList)
  return(ImputeList) }

#' Impute a SummarizedExperiment assay
#'
#' Records feature-level imputation metadata and imputes missing assay values
#' with an MSnbase method or a manually parameterized Gaussian distribution.
#'
#' @param se A `SummarizedExperiment` object.
#' @param fun An imputation method accepted by `MSnbase::impute()`, or `"man"`
#'   for manual Gaussian imputation.
#' @param ... Additional arguments passed to `MSnbase::impute()` or, for
#'   `fun = "man"`, the numeric `shift` and `scale` parameters.
#'
#' @return The `SummarizedExperiment` with an imputed assay and `imputed` and
#'   `num_NAs` row-data columns.
#' @family proteomics functions
#' @importFrom assertthat assert_that
#' @export
impute2 <- function (se, fun = c("bpca", "knn", "QRILC", "MLE", "MinDet","MinProb", "man", "min", "zero", "mixed", "nbavg"), ...)  {
  assertthat::assert_that(inherits(se, "SummarizedExperiment"), is.character(fun))
  fun <- match.arg(fun)
  if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '", deparse(substitute(se)), "'\nRun make_unique() and make_se() to obtain the required columns", call. = FALSE)
  }
  if (!any(is.na(assay(se)))) {
    warning("No missing values in '", deparse(substitute(se)), "'. ", "Returning the unchanged object.", call. = FALSE)
    return(se)
  }
  rowData(se)$imputed <- apply(is.na(assay(se)), 1, any)
  rowData(se)$num_NAs <- rowSums(is.na(assay(se)))
  if (fun == "man") { se <- manual_impute2(se, ...)
  } else {
    .require_optional_package("MSnbase", "impute2()")
    MSnSet_data <- as(se, "MSnSet")
    MSnSet_imputed <- MSnbase::impute(MSnSet_data, method = fun, ...)
    assay(se, withDimnames=FALSE) <- MSnbase::exprs(MSnSet_imputed)
  }
  return(se)}

#' Manually impute missing assay values
#'
#' @param se A `SummarizedExperiment` object.
#' @param shift Number of observed standard deviations by which to left-shift
#'   the imputation distribution.
#' @param scale Scale factor applied to the observed standard deviation.
#'
#' @return The input object with missing assay values imputed.
#' @noRd
manual_impute2 <- function(se, shift = 1.8, scale = 0.3) {
  stopifnot(
    is.numeric(shift), length(shift) == 1L, is.finite(shift),
    is.numeric(scale), length(scale) == 1L, is.finite(scale), scale > 0
  )
  values <- SummarizedExperiment::assay(se)
  observed <- values[!is.na(values)]
  if (length(observed) < 2L) {
    stop("At least two observed assay values are required for manual imputation.",
         call. = FALSE)
  }
  observed_sd <- stats::sd(observed)
  if (!is.finite(observed_sd) || observed_sd == 0) {
    stop("Observed assay values must have non-zero variance.", call. = FALSE)
  }
  missing <- is.na(values)
  values[missing] <- stats::rnorm(
    sum(missing),
    mean = mean(observed) - shift * observed_sd,
    sd = scale * observed_sd
  )
  SummarizedExperiment::assay(se, withDimnames = FALSE) <- values
  se
}

#' Remove datasets with too few quantified features
#'
#' @param VarRMList A named list of `SummarizedExperiment` objects.
#' @param cut The minimum number of assay rows required to retain a dataset.
#'
#' @return The subset of `VarRMList` whose experiments contain more than `cut`
#'   features.
#' @family proteomics functions
#' @export
LowSampleCountRemove <- function(VarRMList, cut){
  Keep <- vapply(
    VarRMList,
    function(x) nrow(SummarizedExperiment::assay(x)) > cut,
    logical(1)
  )
  print(paste(names(VarRMList)[!Keep], "samples removed"))
  return(VarRMList[Keep])}

#' @rdname LowSampleCountRemove
#' @export
LowSampleCountRmove <- function(VarRMList, cut){
  .Deprecated("LowSampleCountRemove", package = "MultiOmicsDataCompile")
  LowSampleCountRemove(VarRMList, cut)
}

#' Calculate differential protein abundance
#'
#' @param DataList A named list of `SummarizedExperiment` objects.
#' @param type One of `"control"`, `"all"`, or `"manual"`, identifying the DEP
#'   comparison mode.
#' @param ComparisonList A named list or vector of manual comparisons. Names
#'   must match elements of `DataList` when `type = "manual"`.
#'
#' @return A named list of differential-analysis `SummarizedExperiment`
#'   objects, including Benjamini-Hochberg adjusted p-values in `rowData`.
#' @family proteomics functions
#' @export
DEAnalysis <- function(DataList, type = "manual", ComparisonList = NULL){
  .require_optional_package("DEP", "DEAnalysis()")
  type <- match.arg(type, c("manual", "control", "all"))
  if (type == "manual" && is.null(ComparisonList)) {
    stop("`ComparisonList` is required when `type = \"manual\"`.", call. = FALSE)
  }
  data_diff <- list()
  for(i in seq_along(DataList)){
    if(type == "control"){
      data_diff[[i]] <- DEP::test_diff(DataList[[i]], type = "control", control = "Plasma_ND") }
    if(type == "all"){
      data_diff[[i]] <- DEP::test_diff(DataList[[i]], type = "all") }
    if(type == "manual"){
      comps <- ComparisonList[names(ComparisonList) ==  names(DataList)[i]]#[[1]] #gsub("_normalized", "", gsub("_Raw", "", gsub("scaled_", "", names(DataList)[i] ) )) ][[1]]
      for(b in seq_along(comps)){
        data_diff[[length(data_diff)+1]] <- suppressWarnings(DEP::test_diff(DataList[[i]], type = "manual", test = comps[b]))
        names(data_diff)[length(data_diff)] <- paste(names(DataList)[i], comps[b], sep = "-")
      } } }
  #### perform manual P-value correction ####
  for(i in seq_along(data_diff)){
    dat <- rowData(data_diff[[i]])
    if(sum(grepl("p.val", colnames(dat))) == 1){
      datsub <-   as.data.frame(dat[,grepl("p.val", colnames(dat))])
      colnames(datsub) <- colnames(dat)[grepl("p.val", colnames(dat))]
    } else { datsub <- dat[,grepl("p.val", colnames(dat))] }
    for(a in seq_len(ncol(datsub))){
      adjusted <- p.adjust(datsub[[a]], method = "BH", n = length(datsub[[a]]))
      if(a == 1){
        adjustPvals <- as.data.frame(adjusted)
      } else { adjustPvals <- cbind(adjustPvals, as.data.frame(adjusted)) }  }
    colnames(adjustPvals) <- gsub("p.val", "BHCorrection", colnames(datsub))
    rowData(data_diff[[i]]) <- cbind(dat, adjustPvals)   }
  return(data_diff) }

#' Save proteomics experiments
#'
#' @param SEList A named list of `SummarizedExperiment` objects.
#' @param Path A path where RDS files are written.
#'
#' @return Invisibly, `NULL`. One RDS file is written for each list element.
#' @family proteomics functions
#' @export
SaveToProteomicDB <- function(SEList, Path){
  for(i in seq_along(SEList)){ saveRDS(SEList[[i]],  file.path(Path, paste(names(SEList[i]), ".rds", sep = ""))) }
  invisible(NULL)
}

#' Load a proteomics experiment
#'
#' @param Path A path containing saved proteomics RDS files.
#' @param Fname The filename of the dataset to load. When omitted, the first RDS
#'   file in `Path` is used.
#'
#' @return A one-element named list containing the selected
#'   `SummarizedExperiment` object.
#' @family proteomics functions
#' @export
ProteomicSELoad <- function(Path, Fname = NULL){
  if (is.null(Fname)) {
    files <- list.files(Path, pattern = "[.]rds$", ignore.case = TRUE)
    if (!length(files)) {
      stop("No RDS files were found in `Path`.", call. = FALSE)
    }
    Fname <- files[[1L]]
  }
  SelectedList <- list()
  SelectedList[[1]] <- readRDS(file.path(Path, Fname))
  names(SelectedList) <- sub("[.]rds$", "", basename(Fname),
                             ignore.case = TRUE)
  return(SelectedList) }

#' Annotate significant differential abundance
#'
#' @param dataDiff A named list of differential-analysis
#'   `SummarizedExperiment` objects.
#' @param alpha A numeric p-value significance threshold.
#' @param lfc A numeric absolute log2 fold-change threshold.
#' @param sigCol One of `"p.val"`, `"p.adj"`, or `"BHCorrection"`, identifying
#'   the significance column family.
#'
#' @return A named list of experiments with significance columns added to
#'   `rowData`.
#' @family proteomics functions
#' @export
SigDEAnnotate <- function(dataDiff, alpha = 0.05, lfc = log2(2), sigCol="p.adj"){
  dep <- list()
  for(i in seq_along(dataDiff)){ dep[[i]] <- add_rejections2(dataDiff[[i]], alpha = alpha, lfc = lfc, sigCol=sigCol) }
  names(dep) <- names(dataDiff)
  return(dep) }

#' Annotate significance in one differential-abundance result
#'
#' @param diff A differential-analysis `SummarizedExperiment` object.
#' @param alpha A numeric p-value significance threshold.
#' @param lfc A numeric absolute log2 fold-change threshold.
#' @param sigCol One of `"p.val"`, `"p.adj"`, or `"BHCorrection"`, identifying
#'   the significance column family.
#'
#' @return `diff` with contrast-specific significance columns added to
#'   `rowData`.
#' @family proteomics functions
#' @importFrom assertthat assert_that
#' @export
add_rejections2 <- function (diff, alpha = 0.05, lfc = 1, sigCol="p.adj")  { # "p.val"
  if (is.integer(alpha)){alpha <- as.numeric(alpha)}
  if (is.integer(lfc)){lfc <- as.numeric(lfc)}
  assertthat::assert_that(inherits(diff, "SummarizedExperiment"), is.numeric(alpha), length(alpha) == 1, is.numeric(lfc), length(lfc) == 1)
  row_data <- rowData(diff, use.names = FALSE) %>% as.data.frame()
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '", deparse(substitute(diff)), "'\nRun make_unique() and make_se() to obtain the required columns", call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '", deparse(substitute(diff)), "'\nRun test_diff() to obtain the required columns", call. = FALSE)
  }
  cols_p <- grep(sigCol, colnames(row_data))
  cols_diff <- grep("_diff", colnames(row_data))
  if (length(cols_p) == 1) {
    rowData(diff)$significant <- row_data[, cols_p] <= alpha & abs(row_data[, cols_diff]) >= lfc
    rowData(diff)$contrast_significant <- rowData(diff, use.names = FALSE)$significant
    colnames(rowData(diff))[ncol(rowData(diff, use.names = FALSE))] <- gsub(sigCol, "significant", colnames(row_data)[cols_p])
  }
  if (length(cols_p) > 1) {
    p_reject <- row_data[, cols_p] <= alpha
    p_reject[is.na(p_reject)] <- FALSE
    diff_reject <- abs(row_data[, cols_diff]) >= lfc
    diff_reject[is.na(diff_reject)] <- FALSE
    sign_df <- p_reject & diff_reject
    sign_df <- cbind(sign_df, significant = apply(sign_df, 1, function(x) any(x)))
    colnames(sign_df) <- gsub(paste("_", sigCol, sep = ""), "_significant", colnames(sign_df))
    sign_df <- cbind(name = row_data$name, as.data.frame(sign_df))
    rowData(diff) <- merge(rowData(diff, use.names = FALSE),
                           sign_df, by = "name")
  }
  return(diff) }

#' Create volcano plots for multiple comparisons
#'
#' @param DEPList A named list of differential-analysis
#'   `SummarizedExperiment` objects.
#' @param ComparisonList Deprecated and ignored; comparison names are derived
#'   from `DEPList`.
#' @param ymin,ymax Numeric y-axis limits.
#' @param xmin,xmax Numeric x-axis limits.
#' @param addNames A logical value indicating whether significant feature names
#'   should be added.
#' @param Adjusted A logical value indicating whether DEP-adjusted p-values
#'   should be plotted.
#' @param BHadjust A logical value indicating whether separately calculated
#'   Benjamini-Hochberg-adjusted p-values should be plotted.
#'
#' @return A named list of `ggplot` objects.
#' @family exploration functions
#' @export
VolcanoPlot <- function(DEPList, ComparisonList = NULL, ymin = 0, ymax = 15,
                        xmin = -3, xmax = 3, addNames = TRUE,
                        Adjusted = TRUE, BHadjust = FALSE){
  PlotList <- list()
  comps <- names(DEPList)# for(i in seq_along(DEPList)){ # comps <- ComparisonList[names(ComparisonList) == names(DEPList)[i]][[1]]#  gsub("_normalized", "", gsub("_Raw", "", gsub("scaled_", "", names(DEPList)[i] ) )) ][[1]]
  for(b in seq_along(comps)){
    tempLit <- list()
    nam <- names(DEPList)[b]#paste(names(DEPList)[i], comps[b], sep = "-")
    tempLit[[1]] <- plot_volcano2(dep = DEPList[[b]], contrast = gsub(".+-", "", comps[b]), label_size = 3, add_names = addNames, adjusted = Adjusted, BHadjusted = BHadjust) + ggtitle(nam) +
      ylim(ymin,ymax) + xlim(xmin,xmax)
    names(tempLit) <- nam
    PlotList <- c(PlotList, tempLit)
  } # }
  return(PlotList) }

#' Create a differential-abundance volcano plot
#'
#' @param dep A differential-analysis `SummarizedExperiment` object.
#' @param contrast A single comparison name present in the object's row data.
#' @param label_size The numeric size of feature labels.
#' @param add_names A logical value indicating whether significant feature names
#'   should be added.
#' @param adjusted A logical value indicating whether DEP-adjusted p-values
#'   should be used.
#' @param plot A logical value indicating whether to return a plot (`TRUE`) or
#'   the underlying plotting data (`FALSE`).
#' @param BHadjusted A logical value indicating whether separately calculated
#'   Benjamini-Hochberg-adjusted p-values should be used.
#'
#' @return A `ggplot` object when `plot = TRUE`; otherwise, a data frame with
#'   protein, fold-change, p-value, and significance columns.
#' @family exploration functions
#' @importFrom assertthat assert_that
#' @export
plot_volcano2 <- function (dep, contrast, label_size = 3, add_names = TRUE,
                           adjusted = FALSE, plot = TRUE,
                           BHadjusted = FALSE) {
  .require_optional_package("DEP", "plot_volcano2()")
  if (is.integer(label_size))
    label_size <- as.numeric(label_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(contrast), length(contrast) == 1, is.numeric(label_size),
                          length(label_size) == 1, is.logical(add_names), length(add_names) ==
                            1, is.logical(adjusted), length(adjusted) == 1, is.logical(plot),
                          length(plot) == 1)
  row_data <- rowData(dep, use.names = FALSE)
  if(any(!c("name", "ID") %in% colnames(row_data))){
    stop(paste0("'name' and/or 'ID' columns are not present in '", deparse(substitute(dep)), "'.\nRun make_unique() to obtain required columns."), call. = FALSE)
  }
  if(length(grep("_p.adj|_diff", colnames(row_data))) < 1){
    stop(paste0("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '", deparse(substitute(dep)), "'.\nRun test_diff() to obtain the required columns."), call. = FALSE)
  }
  if(length(grep("_significant", colnames(row_data))) < 1){
    stop(paste0("'[contrast]_significant' columns are not present in '", deparse(substitute(dep)), "'.\nRun add_rejections() to obtain the required columns."), call. = FALSE)
  }
  if(length(grep(paste(contrast, "_diff", sep = ""), colnames(row_data))) ==  0) {
    valid_cntrsts <- row_data %>% data.frame() %>% select(ends_with("_diff")) %>% colnames(.) %>% gsub("_diff", "", .)
    valid_cntrsts_msg <- paste0("Valid contrasts are: '", paste0(valid_cntrsts, collapse = "', '"), "'")
    stop("Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n",
         valid_cntrsts_msg, call. = FALSE) }
  diff <- grep(paste(contrast, "_diff", sep = ""), colnames(row_data))
  if (adjusted) { p_values <- grep(paste(contrast, "_p.adj", sep = ""), colnames(row_data))
  } else if(BHadjusted){ p_values <- grep(paste(contrast, "_BHCorrection", sep = ""), colnames(row_data))
  } else {
    p_values <- grep(paste(contrast, "_p.val", sep = ""), colnames(row_data)) }
  signif <- grep(paste(contrast, "_significant", sep = ""), colnames(row_data))
  df <- data.frame(x = row_data[, diff], y = -log10(row_data[, p_values]), significant = row_data[, signif], name = row_data$name) %>%
    filter(!is.na(significant)) %>% arrange(significant)
  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)
  p <- ggplot(df, aes(x, y)) + geom_vline(xintercept = 0) +
    geom_point(aes(col = significant)) + geom_text(data = data.frame(),
                                                   aes(x = c(Inf, -Inf), y = c(-Inf, -Inf), hjust = c(1, 0), vjust = c(-1, -1), label = c(name1, name2), size = 5,
                                                       fontface = "bold")) + labs(title = contrast, x = expression(log[2] ~ "Fold change")) + DEP::theme_DEP1() + theme(legend.position = "none") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "grey"))
  if (add_names) {
    p <- p + ggrepel::geom_text_repel(data = filter(df, significant),
                                      aes(label = name), size = label_size, box.padding = unit(0.1, "lines"), point.padding = unit(0.1, "lines"),
                                      segment.size = 0.5) }
  if (adjusted) { p <- p + labs(y = expression(-log[10] ~ "Adjusted p-value"))
  } else if(BHadjusted){ p <- p + labs(y = expression(-log[10] ~ "BH Adjusted p-value"))
  } else { p <- p + labs(y = expression(-log[10] ~ "P-value")) }
  if (plot) { return(p) } else {
    df <- df %>% select(name, x, y, significant) %>% arrange(desc(x))
    colnames(df)[c(1, 2, 3)] <- c("protein", "log2_fold_change", "p_value_-log10")
    if (adjusted) { colnames(df)[3] <- "adjusted_p_value_-log10" }
    if (BHadjusted) { colnames(df)[3] <- "BH_adjusted_p_value_-log10" }
    return(df) } }

#' List protein names in saved proteomics datasets
#'
#' @param fPath A path containing proteomic `SummarizedExperiment` RDS files.
#'
#' @return A sorted character vector of unique protein names.
#' @family proteomics functions
#' @export
ProteomicProteinName <- function(fPath){
  ProtDatasetSelection <- list.files(fPath)
  ProteomicNames <- NULL
  for(i in seq_along(ProtDatasetSelection)){
    ProteomicNames <- unique( c(ProteomicNames, row.names(rowData(readRDS(file.path(fPath, ProtDatasetSelection[i]))))))
  }
  ProteomicNames <- ProteomicNames[order(ProteomicNames)]
  return(ProteomicNames)
}



#################################################################
#### Function that combines all raw data into one data table ####
#################################################################
#' Compile raw transcript-expression data
#'
#' Combines harmonized array and RNA-seq expression files into a gene-level
#' table while preserving species-specific stable feature identifiers.
#'
#' @param Fpath A path containing raw expression RDS files.
#' @param outPath A file path where the compiled tab-separated table is written.
#' @param StartAt The one-based file index at which processing starts.
#' @param sleep Deprecated and ignored. Output is now written once after all
#'   selected datasets have been incorporated.
#' @param GTFHumanFpath,GTFMouseFpath Deprecated and ignored. The active
#'   workflow no longer calculates FPKM from genomic spans.
#' @param overview An optional dataset manifest used to identify each dataset's
#'   technology and species. When omitted, species is inferred from identifiers
#'   and technology is inferred conservatively from the filename.
#'
#' @return Invisibly, the compiled data table.
#' @family transcriptomics functions
#' @import data.table
#' @export
RawDataCompile <- function(Fpath, outPath,
                           StartAt = 1, sleep = 0, GTFHumanFpath = NULL,
                           GTFMouseFpath = NULL, overview = NULL){
  .raw_data_compile_impl(
    Fpath = Fpath, outPath = outPath, StartAt = StartAt, overview = overview
  )
}
#############################################################
#### Function that compiles all meta data into one table ####
#############################################################
#' Compile sample metadata
#'
#' Harmonizes SRA run tables and GEO microarray sample metadata, standardizes
#' field names, and adds disease and tissue annotations from the overview table.
#' Supplying `conditions` also attaches the sample-level `condition` column that
#' [GenerateRawSE()] requires for RNA-seq datasets. Neither repository records a
#' usable condition, so it is taken from the manifest's comparison codes rather
#' than inferred from sample order.
#'
#' @param RNAseqFilePath A path containing RNA-seq metadata tables.
#' @param ArrayFilePath A path containing saved microarray metadata objects.
#' @param overview A data table with dataset `ID`, `Tissue`, and `Disease`
#'   columns.
#' @param conditions An optional sample-condition table produced by
#'   [SampleConditions()], or a path to the directory of raw-expression files
#'   from which one should be derived. When `NULL`, no `condition` column is
#'   added and [GenerateRawSE()] will stop for any RNA-seq dataset.
#'
#' @return A data frame with one row per raw-data sample.
#' @family transcriptomics functions
#' @import data.table
#' @export
MetDataCompile <- function(RNAseqFilePath, ArrayFilePath, overview,
                           conditions = NULL){
  #### RNAseq data ####
  files <- list.files(RNAseqFilePath, full.names = TRUE)
  files <- files[file.exists(files) & !dir.exists(files)]
  dt1 <- data.table()
  for(i in seq_along(files)){
    t <- fread(files[i]);
    t$Dataset <- sub("_.*$", "", basename(files[i]))
    t$RawColumnNames <- paste(t$Dataset, t$`Sample Name`, sep = "_")
    dt1 <- rbind(dt1, t, fill = TRUE)}
  if ("Assay Type" %in% names(dt1)) {
    dt1 <- dt1[`Assay Type` == "RNA-Seq",]
  }
  #### Array data ####
  files <- list.files(ArrayFilePath, pattern = "[.]rds$", full.names = TRUE)
  files <- files[file.exists(files) & !dir.exists(files)]
  dt2 <- data.table()
  for(i in seq_along(files)){
    # t <- fread(files[i]);
    t <- readRDS(files[i]);
    t$Dataset <- sub("_.*$", "", basename(files[i]))
    t$RawColumnNames <- paste(t$Dataset, t$geo_accession, sep = "_")
    dt2 <- rbind(dt2, t, fill = TRUE)}
  dt2 <- as.data.table(dt2)
  # Assigning a column to an empty data table would create a phantom
  # all-missing sample row, so only label the array table when it has rows.
  if (nrow(dt2)) {
    dt2$`Assay Type` <- "Array"
  }
  #### Combine together ####
  dt <- rbind(dt1, dt2, fill = TRUE)
  if (!nrow(dt)) {
    stop("No RNA-seq or array metadata files were found.", call. = FALSE)
  }
  #### clean up column names ####
  cname <- colnames(dt)
  cname <- gsub(":.+", "", cname)
  cname <- gsub("erythrocyte  eicosapentaenoic acid (% of total lipids)", "erythrocyte  eicosapentaenoic acid % of total lipids", cname, fixed = TRUE)
  cname <- gsub("erythrocyte arachidonic acid (% of total lipids)", "erythrocyte arachidonic acid % of total lipids", cname, fixed = TRUE)
  cname <- gsub("erythrocyte docosahexaenoic acid (% of total lipids)", "erythrocyte docosahexaenoic acid % of total lipids", cname, fixed = TRUE)
  cname <- gsub("liver arachidonic acid (% of diacylglycerols)", "liver arachidonic acid % of diacylglycerols", cname, fixed = TRUE)
  cname <- gsub("liver arachidonic acid (% of phospholipids)", "liver arachidonic acid % of phospholipids", cname, fixed = TRUE)
  cname <- gsub("liver arachidonic acid (% of triacylglycerols)", "liver arachidonic acid % of triacylglycerols", cname, fixed = TRUE)
  cname <- gsub("liver docosahexaenoic acid (% of diacylglycerols)", "liver docosahexaenoic acid % of diacylglycerols", cname, fixed = TRUE)
  cname <- gsub("liver docosahexaenoic acid (% of phospholipids)", "liver docosahexaenoic acid % of phospholipids", cname, fixed = TRUE)
  cname <- gsub("liver docosahexaenoic acid (% of triacylglycerols)", "liver docosahexaenoic acid % of triacylglycerols", cname, fixed = TRUE)
  cname <- gsub("liver eicosapentaenoic acid (% of diacylglycerols)", "liver eicosapentaenoic acid % of diacylglycerols", cname, fixed = TRUE)
  cname <- gsub("liver eicosapentaenoic acid (% of phospholipids)", "liver eicosapentaenoic acid % of phospholipids", cname, fixed = TRUE)
  cname <- gsub("liver eicosapentaenoic acid (% of triacylglycerols)", "liver eicosapentaenoic acid % of triacylglycerols", cname, fixed = TRUE)
  cname <- gsub(" \\(.+\\)", "", cname)
  cname <- gsub("\\..+", "", cname)
  cname <- make.names(tolower(cname), unique = FALSE)
  setnames(dt, colnames(dt), cname)
  dt <- dt[,!colnames(dt) %in% c("characteristics_ch1"), with = FALSE,]
  #### merge duplicated column names into one column ####
  cdups <- unique(colnames(dt)[duplicated(colnames(dt))])
  for(a in seq_along(cdups)){
    #### remove duplicated column names ####
    MerAtrribute <- NULL
    temp <- dt[,colnames(dt) %in% cdups[a], with = FALSE]
    for(i in seq_len(nrow(temp))){
      rec <- NULL
      for(b in seq_len(ncol(temp[i,]))){
        rec <- c(rec, temp[i,][[b]])#, temp[i,][[2]])
      }
      if(sum(is.na(rec)) == length(rec)){
        MerAtrribute[i] <- NA
      } else {
        MerAtrribute[i] <- paste(rec[!is.na(rec)], collapse = " / ")
      }
    }
    dt <- dt[,!colnames(dt) %in% cdups[a], with = FALSE]
    dt[[cdups[a]]] <- MerAtrribute
        setnames(dt, colnames(dt), sub("[.][0-9]+$", "", colnames(dt)))
  }
  duplicated_samples <- unique(dt$rawcolumnnames[duplicated(dt$rawcolumnnames)])
  if (length(duplicated_samples)) {
    stop(
      "Duplicated sample identifiers were found in the metadata: ",
      paste(utils::head(duplicated_samples, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  df <- as.data.frame(dt)
  rownames(df) <- df$rawcolumnnames
  #### remove columns containing only missing values ####
  df <- df[,!apply(df, 2, function(x){sum(is.na(x))}) == nrow(df)]
  #### Add Tissue and disease annotations ####
  df <- df[df$dataset %in% unique(overview$ID),]
  df$disease <- NULL; df$tissue <- NULL
  df$disease <- NA; df$tissue <- NA
  overview_df <- as.data.frame(overview, stringsAsFactors = FALSE)
  O2 <- unique(overview_df[, intersect(
    c("ID", "Tissue", "Disease", "Species", "Technology"),
    names(overview_df)
  ), drop = FALSE])
  #### add tissue and disease columns ####
  tempID <- unique(O2$ID)
  for(i in seq_along(tempID)){
    tempDF <- df[df$dataset == tempID[i],]
    if(nrow(tempDF) > 0){
      tempAn <- O2[O2$ID == tempID[i], , drop = FALSE]
      tempDF$disease <- tempAn$Disease[1]
      tempDF$tissue <- tempAn$Tissue[1]
      if ("Species" %in% names(tempAn)) tempDF$species <- tempAn$Species[1]
      if ("Technology" %in% names(tempAn)) tempDF$technology <- tempAn$Technology[1]
      df <- rbind(df[!df$dataset == tempID[i],], tempDF) }
  }
  #### Attach sample conditions ####
  df <- .attach_sample_conditions(df, conditions, overview)
  return(df)
}

# Join manifest-derived sample conditions onto compiled sample metadata.
.attach_sample_conditions <- function(df, conditions, overview) {
  if (is.null(conditions)) {
    message(
      "No `conditions` were supplied. GenerateRawSE() requires a `condition` ",
      "column for every RNA-seq dataset; derive one with SampleConditions()."
    )
    return(df)
  }
  if (is.character(conditions)) {
    if (length(conditions) != 1L) {
      stop("`conditions` must be a single path or a data frame.", call. = FALSE)
    }
    conditions <- SampleConditions(DBPath = conditions, overview = overview)
  }
  conditions <- as.data.frame(conditions, stringsAsFactors = FALSE)
  required <- c("rawcolumnnames", "condition")
  missing <- setdiff(required, names(conditions))
  if (length(missing)) {
    stop(
      "`conditions` is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(conditions$rawcolumnnames)) {
    stop("`conditions` must contain one row per sample.", call. = FALSE)
  }
  supplied <- as.character(
    conditions$condition[match(df$rawcolumnnames, conditions$rawcolumnnames)]
  )
  if ("condition" %in% names(df)) {
    existing <- as.character(df$condition)
    df$condition <- ifelse(is.na(supplied), existing, supplied)
  } else {
    df$condition <- supplied
  }
  unmatched <- df$rawcolumnnames[is.na(df$condition)]
  if (length(unmatched)) {
    message(
      length(unmatched),
      " sample(s) have no condition assignment: ",
      paste(utils::head(unmatched, 5L), collapse = ", "),
      if (length(unmatched) > 5L) " ..." else ""
    )
  }
  df
}

########################################################
#### Generate raw data summarized experiment object ####
########################################################
#' Create raw and normalized expression experiments
#'
#' Aligns sample metadata with a compiled expression table and creates separate
#' per-dataset RNA-seq and microarray experiments. RNA-seq experiments retain
#' raw counts, DESeq2 size-factor normalized counts, and transformed values.
#' A separate exploration experiment contains within-dataset z-scores so that
#' incompatible measurement units are not presented as one abundance assay.
#'
#' @param df A data frame containing sample metadata, with row names matching
#'   raw-data column names.
#' @param ArrayDT A data table returned by [RawDataCompile()].
#' @param overview A data table containing dataset technology, species, tissue,
#'   and disease metadata.
#'
#' @return A list containing `RNAseqSE` and `ArraySE`, named lists of
#'   per-dataset experiments; `ExplorationSE`, a within-dataset standardized
#'   experiment for qualitative visualization; and compatibility aliases
#'   `RawSE`, `ExpressionSE`, and `RPKMSE`. `RPKMSE` no longer contains RPKM or
#'   FPKM values and is retained only to ease migration.
#' @family transcriptomics functions
#' @export
GenerateRawSE <- function(df, ArrayDT, overview){
  .generate_expression_collection(df, ArrayDT, overview)
}
