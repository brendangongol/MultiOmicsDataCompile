
###########################
#### MultiOmics Server ####
###########################

#############################################
#### load Objects and set up environment ####
#############################################
library(shiny); library(shinyWidgets); library(SummarizedExperiment)
library(data.table); library(textclean); library(dplyr)
library(ggplot2); library(plotly); library(RColorBrewer); library(stringr)
suppressWarnings(library(DEP))
library(grid)
# Analysis functions live in the package so that the app and the batch
# workflow cannot drift apart.
library(MultiOmicsDataCompile)

#### loading tables ####
########################
#### Load genes ####
homedir <- getOption(
  "MultiOmicsDataCompile.data_dir",
  Sys.getenv("MULTIOMICS_DATA_DIR", unset = "<path to data base>")
)
if (!dir.exists(homedir)) {
  stop("Set `MultiOmicsDataCompile.data_dir` or `MULTIOMICS_DATA_DIR` to a valid data directory.")
}
DESE <- readRDS(file.path(homedir, "ProcessFiles", "SumarizedExp_DB.rds"))
AvailComps <- data.table(Dataset = gsub("_.+", "", names(assays(DESE))), Comparison = gsub("^.+?_", "", names(assays(DESE))) )
genesAll <- rowData(DESE)$SYMBOL
#### Overview table ####
overview <- fread(file.path(homedir, "OverviewFiles", "GEODataOverview3.csv"), header = TRUE)
overview <- overview[!ID == "",]
overview$Website <- paste("<a href=", overview$Website, "target=_blank>", "GEO_website", "</a>")
overview[grepl("href=  target=_blank", overview$Website),]$Website <- NA
#### Raw data files ####
raw_df = readRDS(file.path(homedir, "ProcessFiles", "expression_norm.v2.RDS"))
genesAll = rowData(raw_df)$SYMBOL
dataAll = unique(colData(raw_df)$dataset)%>%na.omit()
diseaseAll = unique(overview$Disease)
techAll = unique(overview$Technology)
tissueAll = unique(overview$Tissue)
#### Metabolomics data ####
metabolomicsData <- readRDS(file.path(homedir, "Metabolomics", "mtbls298.de.RDS"))
names(metabolomicsData) <- c("Basal_Artery_Basal_Vein", "Insulin_Artery_Insulin_Vein")
#### Proteomics tab ####
DEProtSE <- readRDS(file.path(homedir, "ProcessFiles", "SumarizedProtExp_DB.rds"))
Proteins <- readRDS(file.path(homedir, "OverviewFiles", "ProteomicProteins.RDS"))

AssaySymbols <- function(se, assay_matrix) {
  symbols <- as.character(rowData(se)$SYMBOL[match(rownames(assay_matrix), rownames(se))])
  symbols[is.na(symbols) | !nzchar(symbols)] <- rownames(assay_matrix)[is.na(symbols) | !nzchar(symbols)]
  toupper(symbols)
}

###############################
#### Summary plot function ####
###############################
SummaryPlot <- function(over2, features, attribute){
  if(features == "Studies"){
    unique <- over2[!duplicated(ID),]
    if(attribute == "Tissue"){
      unique2 <- unique
      remove <- NULL
      for(i in seq_len(nrow(unique2))){
        spl <- str_split(unique2$Tissue[i], ",")[[1]]
        if(length(spl) > 1){
          spl <- gsub(" ", "", spl)
          t <- unique2[i,]
          dttemp <- data.table(ID=t$ID,GPLNumber=t$GPLNumber,ComparisonVector=t$ComparisonVector,RawColumnNames=t$RawColumnNames,
                               FCColumnNames=t$FCColumnNames,GenerateSpotfire=t$GenerateSpotfire,DownloadData=t$DownloadData,Title=t$Title,Year=t$Year,`Profiling Resource`=t$`Profiling Resource`,
                               Tissue = spl,Species = t$Species,Disease=t$Disease,`Donor count`=t$`Donor count`,Website=t$Website,DataType=t$DataType,
                               Technology=t$Technology,DataIncorporated=t$DataIncorporated, Description=t$Description, Resources=t$Resources )
          unique2 <- rbind(unique2, dttemp)
          remove <- c(remove, i) }}
      if(!is.null(remove)){ unique2 <- unique2[!c(remove),] }
      unique2$Tissue <- gsub("Wholeblood", "Whole blood", unique2$Tissue)
      DTDiseaseBar <- unique2[,.(Count=length(ID)), by = Tissue]
      p <- ggplotly(ggplot(DTDiseaseBar, aes(x=Tissue, y=Count, fill=Tissue))+
                      geom_bar(stat = "identity", color="black") +
                      scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(8, "Accent"), RColorBrewer::brewer.pal(12, "Paired"))) +
                      theme(legend.position = "none",
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            panel.background = element_blank(), axis.line = element_line(colour="black"),
                            axis.text.x = element_text(angle = 90, hjust = 1)) +
                      ylab("Number of Studies") + xlab(colnames(DTDiseaseBar)[1]), tooltip = c("x", "y") ) }
    if(attribute == "Disease"){
      DTDiseaseBar <- unique[,.(Count=length(ID)), by = Disease]
      p <- ggplotly(
        ggplot(DTDiseaseBar, aes(x=Disease, y=Count, fill=Disease))+
          geom_bar(stat = "identity", color="black") +
          scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(8, "Accent"), RColorBrewer::brewer.pal(12, "Paired"))) +
          theme(legend.position = "none",
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour="black"),
                axis.text.x = element_text(angle = 90, hjust = 1)) +
          ylab("Number of Studies") + xlab(colnames(DTDiseaseBar)[1]), tooltip = c("x", "y") ) }
    if(attribute == "Species"){
      unique2 <- unique
      remove <- NULL
      for(i in seq_len(nrow(unique2))){
        spl <- str_split(unique2$Species[i], ",")[[1]]
        if(length(spl) > 1){
          spl <- gsub(" ", "", spl)
          t <- unique[i,]
          dttemp <- data.table(ID=t$ID,GPLNumber=t$GPLNumber,ComparisonVector=t$ComparisonVector,RawColumnNames=t$RawColumnNames,
                               FCColumnNames=t$FCColumnNames,GenerateSpotfire=t$GenerateSpotfire,DownloadData=t$DownloadData,Title=t$Title,Year=t$Year,`Profiling Resource`=t$`Profiling Resource`,
                               Tissue = t$Tissue,Species = spl,Disease=t$Disease,`Donor count`=t$`Donor count`,Website=t$Website,DataType=t$DataType,
                               Technology=t$Technology,DataIncorporated=t$DataIncorporated, Description=t$Description, Resources=t$Resources)
          unique2 <- rbind(unique2, dttemp)
          remove <- c(remove, i) } }
      if(!is.null(remove)){ unique2 <- unique2[!c(remove),] }
      DTDiseaseBar <- unique2[,.(Count=length(ID)), by = Species]
      p <- ggplotly( ggplot(DTDiseaseBar, aes(x=Species, y=Count, fill=Species))+
                       geom_bar(stat = "identity", color="black") +
                       scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(8, "Accent"), RColorBrewer::brewer.pal(12, "Paired"))) +
                       theme(legend.position = "none",
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour="black"),
                             axis.text.x = element_text(angle = 90, hjust = 1)) +
                       ylab("Number of Studies") + xlab(colnames(DTDiseaseBar)[1]) , tooltip = c("x", "y") ) }
    if(attribute == "Technology"){
      DTDiseaseBar <- unique[,.(Count=length(ID)), by = Technology]
      p <- ggplotly(
        ggplot(DTDiseaseBar, aes(x=Technology, y=Count, fill=Technology))+
          geom_bar(stat = "identity", color="black") +
          scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(8, "Accent"), RColorBrewer::brewer.pal(12, "Paired"))) +
          theme(legend.position = "none",
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour="black"),
                axis.text.x = element_text(angle = 90, hjust = 1)) +
          ylab("Number of Studies") + xlab(colnames(DTDiseaseBar)[1]), tooltip = c("x", "y") ) }
    if(attribute == "Year"){
      DTDiseaseBar <- unique[,.(Count=length(ID)), by = Year][,Year := as.character(Year)]
      p <- ggplotly(
        ggplot(DTDiseaseBar, aes(x=Year, y=Count, fill=Year))+
          geom_bar(stat = "identity", color="black") +
          scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(8, "Accent"), RColorBrewer::brewer.pal(12, "Paired"))) +
          theme(legend.position = "none",
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour="black"),
                axis.text.x = element_text(angle = 90, hjust = 1)) +
          ylab("Number of Studies") + xlab(colnames(DTDiseaseBar)[1]), tooltip = c("x", "y") ) } }
  return(p) }

###################################
#### DE table display function ####
###################################
DETableDisplay <- function(DESE, Dataset, Comparison, return = "NA"){
  if(!Dataset == "Select Dataset"){
    Fname <- paste(Dataset, Comparison, sep = "_")
    assaySelect <- assays(DESE)[grepl(Fname, gsub("-", "_", names(assays(DESE))))]
    tempDT <- assaySelect[[1]]
    tempDT$SYMBOL <- AssaySymbols(DESE, tempDT)
    tempDT <- as.data.table(tempDT)
    FDT <- tempDT[complete.cases(tempDT),]
    if(return == "Genes"){
      return(FDT$SYMBOL)
    } else { return(as.data.frame(FDT)) } } }

###############################
#### Volcano plot function ####
###############################
VolcanoPlotR <- function(DESE,Dataset,Comparison,sigColSelect,FCcut,Pcut,filterBy,Gene){
  if(is.numeric(FCcut) && length(FCcut) == 1 && is.finite(FCcut) && FCcut >= 0){
    if(!Dataset == "Select Dataset"){
      Fname <- paste(Dataset, Comparison, sep = "_")
      assaySelect <- assays(DESE)[grepl(Fname, gsub("-", "_", names(assays(DESE))))]
      tempDT <- assaySelect[[1]]
      tempDT$SYMBOL <- AssaySymbols(DESE, tempDT)
      tempDT <- as.data.table(tempDT)
      tempDT <- tempDT[complete.cases(tempDT),]
      tempDT[,Selected:= "Not selected"]
      tempDT[,Selected:=ifelse( ( abs(tempDT[["logFC"]]) > FCcut & tempDT[[sigColSelect]] < Pcut ), "Selected", "Not selected")]
      if(sigColSelect == "Pvalue"){
        if(filterBy == "Significance"){
          p <- ggplot(tempDT, aes(x=logFC,y=-log10(Pvalue), color = Selected, text = SYMBOL)) +
            xlab("log2(FC)")+ ylab(paste("-log10(", sigColSelect, ")", sep = "")) +
            ggtitle(names(assaySelect)) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour="black")) +
            geom_point(data=subset(tempDT, Selected == "Not selected"),
                       aes(x=logFC,y=-log10(Pvalue) ),
                       shape=21, size=2, colour="gray", fill = "gray") +
            geom_point(data=subset(tempDT, Selected == "Selected"),
                       aes(x=logFC,y=-log10(Pvalue) ),
                       shape=21, size=3, colour="black", fill = "#E31A1C")
        } else {
          p <- ggplot(tempDT, aes(x=logFC,y=-log10(Pvalue), color = Selected, text = SYMBOL)) +
            xlab("log2(FC)")+ ylab(paste("-log10(", sigColSelect, ")", sep = "")) +
            ggtitle(names(assaySelect)) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour="black")) +
            geom_point(data=subset(tempDT, !SYMBOL == toupper(Gene)),
                       aes(x=logFC,y=-log10(Pvalue) ),
                       shape=21, size=2, colour="gray", fill = "gray") +
            geom_point(data=subset(tempDT, SYMBOL == toupper(Gene)),
                       aes(x=logFC,y=-log10(Pvalue) ),
                       shape=21, size=3, colour="black", fill = "#E31A1C") }
      } else if(sigColSelect == "AdjPValue"){
        if(filterBy == "Significance"){
          p <- ggplot(tempDT, aes(x=logFC,y=-log10(AdjPValue), color = Selected, text = SYMBOL)) +
            xlab("log2(FC)")+ ylab(paste("-log10(", sigColSelect, ")", sep = "")) +
            ggtitle(names(assaySelect)) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour="black")) +
            geom_point(data=subset(tempDT, Selected == "Not selected"),
                       aes(x=logFC,y=-log10(AdjPValue) ),
                       shape=21, size=2, colour="gray", fill = "gray") +
            geom_point(data=subset(tempDT, Selected == "Selected"),
                       aes(x=logFC,y=-log10(AdjPValue) ),
                       shape=21, size=3, colour="black", fill = "#E31A1C")
        } else {
          p <- ggplot(tempDT, aes(x=logFC,y=-log10(AdjPValue), color = Selected, text = SYMBOL)) +
            xlab("log2(FC)")+ ylab(paste("-log10(", sigColSelect, ")", sep = "")) +
            ggtitle(names(assaySelect)) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour="black")) +
            geom_point(data=subset(tempDT, !SYMBOL == toupper(Gene)),
                       aes(x=logFC,y=-log10(AdjPValue) ),
                       shape=21, size=2, colour="gray", fill = "gray") +
            geom_point(data=subset(tempDT, SYMBOL == toupper(Gene)),
                       aes(x=logFC,y=-log10(AdjPValue) ),
                       shape=21, size=3, colour="black", fill = "#E31A1C") }}
      return(p) } } }

############################
#### DEG count function ####
############################
DEGcount <- function(DESE = DESE, Dataset,Comparison,sigColSelect = "Pvalue",filterBy,Gene,FCcut,Pcut){
  if(is.numeric(FCcut) && length(FCcut) == 1 && is.finite(FCcut) && FCcut >= 0){
    if(!Dataset == "Select Dataset"){
      Fname <- paste(Dataset, Comparison, sep = "_")
      assaySelect <- assays(DESE)[grepl(Fname, gsub("-", "_", names(assays(DESE)))   )  ]
      tempDT <- assaySelect[[1]]
      tempDT$SYMBOL <- AssaySymbols(DESE, tempDT)
      tempDT <- as.data.table(tempDT)
      tempDT <- tempDT[complete.cases(tempDT),]
      #### update colors ####
      tempDT[,Selected:= "Not selected"]
      tempDT[,Selected:=ifelse( ( abs(tempDT[["logFC"]]) > FCcut & tempDT[[sigColSelect]] < Pcut ), "Selected", "Not selected")]
      if(filterBy == "Significance"){
        statement <- as.data.frame(list(c(
          paste("Total genes:", nrow(tempDT)),
          paste("Total Selected:", nrow(tempDT[Selected == "Selected",])),
          paste("Total Not Selected", nrow(tempDT[Selected == "Not selected",])),
          paste("Total Up-regulated Selected:", nrow(tempDT[Selected == "Selected" & logFC > 0,])),
          paste("Total Down-regulated Selected:", nrow(tempDT[Selected == "Selected" & logFC < 0,])))  ))
        colnames(statement) <- " "
      } else {
        statement <- as.data.frame(paste("Selected gene is Highlighted"))
        colnames(statement) <- " " }
      return(statement) } } }

######################################
#### Multi-Study Heatmap function ####
######################################
CrossDataHeat <- function(DESE, GeneSelection, Dataselection, ScaleData, plottype = "Heat", FCCutoff, PCutoff, SigCol){
  if(is.numeric(FCCutoff) && length(FCCutoff) == 1 && is.finite(FCCutoff) && FCCutoff >= 0){
    if(!is.null(Dataselection)){
      #### Select Genes ####
      SEsub <- DESE[rowData(DESE)$SYMBOL %in% GeneSelection]
      #### select assays ####
      assaySelect <- assays(SEsub)[names(assays(SEsub)) %in% Dataselection]
      if (!length(assaySelect)) return(NULL)
      #### loop through remaining assays and combine results ####
      feature_labels <- make.unique(AssaySymbols(SEsub, assaySelect[[1]]))
      for(i in seq_along(assaySelect)){ dat <- assaySelect[[i]]
      rownames(dat) <- feature_labels
      colnames(dat) <- paste(colnames(dat), names(assaySelect)[i], sep = "_")
      if(i == 1){ compiledDF <-dat
      } else { compiledDF <- merge(compiledDF, dat, by = "row.names", all = TRUE)
      row.names(compiledDF) <- compiledDF$Row.names
      compiledDF$Row.names <- NULL  } }
      #### create heatmap ####
      ########################
      if(plottype == "Heat"){
        selcomps <- unique(mgsub(colnames(compiledDF), c("logFC_", "AdjPValue_", "Pvalue_"), c("","","")))
        DFsub3 <- data.table()
        for(b in seq_along(selcomps)){
          temp <- compiledDF[,grepl(selcomps[b], colnames(compiledDF))]
          temp <- temp[,grepl(paste("logFC", SigCol, sep = "|"), colnames(temp))]
          colnames(temp) <- c("logFC", "Significance")
          temp$SYMBOL <- rownames(temp)
          temp$Comparison <- selcomps[b]
          DFsub3 <- rbind(DFsub3, temp)
        }
        d1 <- DFsub3[abs(logFC) > FCCutoff & Significance < PCutoff,]
        d2 <- DFsub3[!(paste(DFsub3$SYMBOL, DFsub3$Comparison, sep = "") %in% paste(d1$SYMBOL, d1$Comparison, sep = "")),][, `:=` (logFC = NA, Significance = NA)][]
        d3 <- rbind(d1, d2)
        d3 <- d3[,c("logFC", "SYMBOL", "Comparison"), with = FALSE]
        DFsub <- reshape2::dcast(d3, SYMBOL~Comparison, value.var = "logFC")
        rownames(DFsub) <- DFsub$SYMBOL
        DFsub$SYMBOL <- NULL
        DFsub2 <- DFsub
        #### scale data ####
        if(ScaleData){
          transMat <- t(DFsub2)
          Rnames <- row.names(transMat)
          transMat <- suppressWarnings(apply(transMat, 2, as.numeric))
          rowScaled <- apply(as.matrix(transMat), 2, scale)
          rownames(rowScaled) <- Rnames
          DFsub2 <- as.data.frame(t(rowScaled))
        }
        DFsub2$names <- row.names(DFsub)
        mel <- as.data.table(reshape2::melt(DFsub2))
        if(ScaleData){
          setnames(mel, c("names", "variable", "value"), c("Gene", "Dataset", "Scaled log2(Fold Change)"))
          p <- ggplot(mel, aes(x = Dataset, y=Gene, fill = `Scaled log2(Fold Change)`)) +
            geom_tile(color = "white", lwd = 0.75, linetype = 1) +
            scale_fill_gradient2(low = "#075AFF",
                                 high= "#FF0000") +
            coord_fixed() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
        } else {
          setnames(mel, c("names", "variable", "value"), c("Gene", "Dataset", "log2(Fold Change)"))
          p <- ggplot(mel, aes(x = Dataset, y=Gene, fill = `log2(Fold Change)`)) +
            geom_tile(color = "white", lwd = 0.75, linetype = 1) +
            scale_fill_gradient2(low = "#075AFF",
                                 high= "#FF0000") +
            coord_fixed() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
        }
        return(p) }
      #### Obtain counts of selected genes and return a barplot ####
      if(plottype == "Bar"){
        selcomps <- unique(mgsub(colnames(compiledDF), c("logFC_", "AdjPValue_", "Pvalue_"), c("","","")))
        DFsub3 <- data.table()
        for(b in seq_along(selcomps)){
          temp <- compiledDF[,grepl(selcomps[b], colnames(compiledDF))]
          temp <- temp[,grepl(paste("logFC", SigCol, sep = "|"), colnames(temp))]
          colnames(temp) <- c("logFC", "Significance")
          temp$SYMBOL <- rownames(temp)
          temp$Comparison <- selcomps[b]
          DFsub3 <- rbind(DFsub3, temp)
        }
        great <- DFsub3[logFC > FCCutoff & Significance < PCutoff,]
        FCGreat <- great[, .(AveFC = mean(logFC, na.rm = TRUE)), by = "SYMBOL"]
        setnames(FCGreat, c("SYMBOL"), c("GeneName"))
        if(nrow(great) > 0){
          great <- data.table(table(great$SYMBOL)) %>% setnames(c("V1", "N"), c("GeneName", "Nup"))
          great <- merge(great, FCGreat, by = "GeneName")
          GreatCountAll <- great
          great <- great[order(great$Nup, decreasing = TRUE),]
        }
        less <- DFsub3[logFC < -FCCutoff & Significance < PCutoff,]
        FCless <- less[, .(AveFC = mean(logFC, na.rm = TRUE)), by = "SYMBOL"]
        setnames(FCless, c("SYMBOL"), c("GeneName"))
        if(nrow(less) > 0){
          less <- data.table(table(less$SYMBOL)) %>% setnames(c("V1", "N"), c("GeneName", "Ndown"))
          less <- merge(less, FCless, by = "GeneName")
          lessCountAll <- less
          less <- less[order(less$Ndown, decreasing = TRUE),]
        }
        if((nrow(great) > 0 & nrow(less) > 0)){
          final <- merge(great, less, by = "GeneName", all = TRUE)
          #### Merge Fold Change ####
          final$AveFC <- 0
          for(b in seq_len(nrow(final))){
            t <- c(final[b,]$AveFC.x, final[b,]$AveFC.y)
            t <- t[!is.na(t)]
            if(length(t) > 1){
              F2 <- data.table()
              for(c in seq_along(t)){
                F2 <- rbind(F2, final[b,])
              }
              F2$AveFC <- t
              final <- rbind(final, F2)
            } else{ final$AveFC[b] <- t[!is.na(t)]
            } }
          final$AveFC.x <- NULL; final$AveFC.y <- NULL
          final <- final[!AveFC == 0,]
          final <- unique(final)
          #### Merge counts ####
          final$direction <- "NA"
          final$counts <- 0
          for(b in seq_len(nrow(final))){
            up <- final$Nup[b]
            down <- final$Ndown[b]
            fc <- final$AveFC[b]
            if(fc > 0){
              final$direction[b] <- "up"
              final$counts[b] <- up
            } else {
              final$direction[b] <- "down"
              final$counts[b] <- down
            } }
          final$Nup <- NULL; final$Ndown <- NULL
          final <- rbind(final[direction == "up",][order(-counts, -AveFC),],
                         final[direction == "down",][order(counts, -AveFC),])
        }
        if(nrow(great) > 0 & nrow(less) == 0){
          final <- great
          final$direction <- "up"
          setnames(final, c("Nup"), c("counts"))
          final <- final[order(-counts, AveFC),]
        }
        if(nrow(great) == 0 & nrow(less) > 0){
          final <- less
          final$direction <- "down"
          setnames(final, c("Ndown"), c("counts"))
          final <- final[order(-counts, AveFC),]
        }
        if(!is.null(final)){
          final$GeneName <- factor(final$GeneName, levels = unique(final$GeneName ))
          p <- ggplot(final, aes(GeneName, AveFC, fill = direction)) +
            geom_bar(position = "dodge", stat = "identity") +
            geom_text(aes(label=counts), vjust=-0.5) +
            theme(axis.text.y = element_text(angle = 0, hjust = 1),
                  axis.text.x = element_text(angle = 90, hjust = 1),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour="black") ) +
            ylab("Average log2(Fold change)")+
            xlab("Gene Name")+
            ggtitle(paste("average fold change across comparisons")) +
            scale_fill_manual(values = RColorBrewer::brewer.pal(11, "RdYlBu")[c(10,2)])
        }
        return(p) }
      if(plottype == "Table"){
        selcomps <- unique(mgsub(colnames(compiledDF), c("logFC_", "AdjPValue_", "Pvalue_"), c("","","")))
        DFsub3 <- data.table()
        for(b in seq_along(selcomps)){
          temp <- compiledDF[,grepl(selcomps[b], colnames(compiledDF))]
          temp <- temp[,grepl(paste("logFC", SigCol, sep = "|"), colnames(temp))]
          colnames(temp) <- c("logFC", "Significance")
          temp$SYMBOL <- rownames(temp)
          temp$Comparison <- selcomps[b]
          DFsub3 <- rbind(DFsub3, temp)
        }
        d1 <- DFsub3[abs(logFC) > FCCutoff & Significance < PCutoff,]
        CompDF <- d1[!is.na(d1$logFC),]
        CompDF <- CompDF[order(SYMBOL),]
        return(as.data.frame(CompDF))
      } } } }

####################################
#### Check NA boxplot NA values ####
####################################
check_na = function(raw_df, select_gene, select_data){
  if(is.null(select_gene) || !length(select_gene)) return(FALSE)
  if(is.null(select_data) || !length(select_data)) return(FALSE)
  samples <- colnames(raw_df)[colData(raw_df)$dataset %in% select_data]
  features <- which(rowData(raw_df)$SYMBOL %in% select_gene)
  if(!length(samples) || !length(features)) return(TRUE)
  any(is.na(assay(raw_df)[features, samples, drop = FALSE]))
}

##########################################################
#### Shared reshaping for the expression violin plots ####
##########################################################
#### Sample annotations are attached only when the compiled metadata
#### actually contain them, so a missing optional field cannot break a plot.
ViolinSampleColumns <- c(
  "dataset", "rawcolumnnames", "condition", "disease", "tissue",
  "species", "technology"
)

ExpressionLongFormat <- function(input_df, select_gene, samples){
  features <- which(rowData(input_df)$SYMBOL %in% select_gene)
  if(!length(features) || !length(samples)) return(NULL)
  long <- MeltExpressionAssay(
    input_df[features, samples],
    assay_name = "Expression",
    row_columns = "SYMBOL",
    col_columns = ViolinSampleColumns
  )
  long <- long[!is.na(long$Expression), , drop = FALSE]
  if(!nrow(long)) return(NULL)
  if(!"condition" %in% names(long)) long$condition <- "All samples"
  long$condition[is.na(long$condition)] <- "Unassigned"
  if(!"dataset" %in% names(long)) long$dataset <- "All datasets"
  long
}

ViolinPalette <- function(n){
  colorRampPalette(brewer.pal(9, "Set3"))(max(as.integer(n), 1L))
}

################################################################
#### Violin plot of standardized expression: single dataset ####
################################################################
box_violin_plot_v2.single = function(input_df, select_gene, select_proj, log10){
  samples <- colnames(input_df)[colData(input_df)$dataset %in% select_proj]
  long <- ExpressionLongFormat(input_df, select_gene, samples)
  if(is.null(long)) return(NULL)
  ggplot(data = long, aes(x = condition, y = Expression, fill = condition))+
    geom_violin(scale = "width", trim = TRUE, alpha = 0.7)+
    geom_boxplot(outlier.shape = NA, coef = 0, fill = "gray", width = 0.3)+
    coord_flip()+
    ylab("Within-dataset expression z-score")+
    xlab("")+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(colour = 'black', size = 10),
          axis.text.x = element_text(colour = 'black', size = 10, vjust = 1, hjust = 1),
          legend.position = "none",
          strip.text.y = element_text(angle = 0, size = 10))+
    stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.2, position = "dodge2")+
    scale_fill_manual(values = ViolinPalette(length(unique(long$condition)))) +
    facet_grid(SYMBOL~., scales = "free", space = "free", drop = FALSE)+
    scale_x_discrete(drop = FALSE)
}

##################################################################
#### Violin plot of standardized expression: multiple datasets ####
##################################################################
box_violin_plot_v2.multi = function(input_df, overview_file, select_gene, select_disease, select_tech, select_tissue, log10){
  select_group = overview_file%>%dplyr::filter(Disease%in%select_disease & Technology%in%select_tech & Tissue%in%select_tissue)%>%dplyr::select(ID)%>%unlist%>%as.vector()%>%unique()
  samples <- colnames(input_df)[colData(input_df)$dataset %in% select_group]
  long <- ExpressionLongFormat(input_df, select_gene, samples)
  if(is.null(long)) return(NULL)
  #### One panel per dataset ####
  long$facet_name <- long$dataset
  ggplot(data = long, aes(x = condition, y = Expression, fill = condition))+
    geom_violin(scale = "width", trim = TRUE, alpha = 0.7)+
    geom_boxplot(outlier.shape = NA, coef = 0, fill = "gray", width = 0.3)+
    ylab("Within-dataset expression z-score")+
    xlab("")+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(colour = 'black', size = 10),
          axis.text.x = element_text(colour = 'black', size = 10),
          legend.position = "none",
          strip.text.y = element_text(angle = 0, size = 10))+
    stat_boxplot(geom = 'errorbar', linetype = 1, width = 0.2, position = "dodge2")+
    scale_fill_manual(values = ViolinPalette(length(unique(long$condition)))) +
    facet_wrap(~facet_name, scales = "free", dir = "v", ncol = 1, strip.position = "right")
}

#################################
#### Load proteomics dataset ####
#################################
#### Loading and significance annotation both come from the package so that
#### the app and the batch workflow apply identical rules.
ProteomicsDataLoad <- function(Path, Fname, alpha, lfc, sigCol){
  SelectedList <- ProteomicSELoad(Path = Path, Fname = paste0(Fname, ".rds"))
  SigDEAnnotate(SelectedList, alpha = alpha, lfc = lfc, sigCol = sigCol)
}

##########################################
#### Proteomics Volcano plot function ####
##########################################
#### Contrast names are read from the stored row data because one saved
#### experiment can hold several comparisons. Plotting itself is delegated to
#### MultiOmicsDataCompile::plot_volcano2().
ProteomicVolcanoPlots <- function(DEPList, addNames = TRUE, sigCol){
  adjusted <- identical(sigCol, "p.adj")
  BHadjusted <- identical(sigCol, "BHCorrection")
  PlotList <- list()
  for(i in seq_along(DEPList)){
    annotations <- colnames(rowData(DEPList[[i]]))
    contrasts <- unique(mgsub(
      annotations[grepl("_vs_", annotations)],
      c("_CI.L", "_CI.R", "_diff", "_p.adj", "_p.val", "_BHCorrection", "_significant"),
      rep("", 7)
    ))
    for(b in seq_along(contrasts)){
      nam <- paste(names(DEPList)[i], contrasts[b], sep = "-")
      PlotList[[nam]] <- plot_volcano2(
        DEPList[[i]], contrast = contrasts[b], label_size = 3,
        add_names = addNames, adjusted = adjusted, BHadjusted = BHadjusted
      ) + ggtitle(nam)
    }
  }
  return(PlotList)
}

########################################################################
#### Count the number of selected differentially expressed proteins ####
########################################################################
ProteomicDEGcount <- function(ProtDE = t){
  tempDT <- rowData(ProtDE[[1]])
  statement <- as.data.frame(list(c(
    paste("Total genes:", nrow(tempDT)),
    paste("Total Selected:", nrow(tempDT[tempDT$significant == TRUE,])),
    paste("Total Not Selected", nrow(tempDT[tempDT$significant == FALSE,])),
    paste("Total Up-regulated Selected:", nrow(tempDT[tempDT$significant == TRUE & tempDT[[grep("diff", colnames(tempDT))]] > 0,])),
    paste("Total Down-regulated Selected:", nrow(tempDT[tempDT$significant == TRUE & tempDT[[grep("diff", colnames(tempDT))]] < 0,])))  ))
  colnames(statement) <- " "
  return(statement)
}

########################################################
#### Proteomic table of the number of selected DEPs ####
########################################################
ProteomicTableDisplay <- function(Path=file.path(homedir, "Proteomic_3"),
                                  Fname=isolate(input$ProteomicDataset),
                                  alpha = isolate(input$ProteomicPval),
                                  lfc = isolate(input$ProteomicFC),
                                  sigCol=isolate(input$ProteomicHypothesisTestDE) ){
  if(!Fname == "Select Dataset"){
    ProtLis = ProteomicsDataLoad(Path=Path, Fname=Fname, alpha=alpha, lfc=lfc, sigCol=sigCol)
    temp <- as.data.frame(rowData(ProtLis[[1]]))
    return(temp)
  } }

################################################
#### Proteomic Multi-Study Heatmap function ####
################################################
CrossDataHeatProteomic <- function(DESE, GeneSelection, Dataselection, ScaleData, plottype = "Heat", FCCutoff, PCutoff, SigCol){
  if(!is.null(Dataselection)){
    #### Select Genes ####
    SEsub <- DESE[rowData(DESE)$SYMBOL %in% GeneSelection]
    #### select assays ####
    assaySelect <- assays(SEsub)[names(assays(SEsub)) %in% Dataselection]
    #### loop through remaining assays and combine results ####
    for(i in seq_along(assaySelect)){ dat <- assaySelect[[i]]
    colnames(dat) <- paste(colnames(dat), names(assaySelect)[i], sep = "_")
    if(i == 1){ compiledDF <-dat
    } else { compiledDF <- merge(compiledDF, dat, by = "row.names", all = TRUE)
    row.names(compiledDF) <- compiledDF$Row.names
    compiledDF$Row.names <- NULL  } }
    #### create heatmap ####
    ########################
    if(plottype == "Heat"){
      selcomps <- unique(mgsub(colnames(compiledDF), c("logFC_", "AdjPValue_", "Pvalue_", "BHCorrection_"), c("","","","")))
      DFsub3 <- data.table()
      for(b in seq_along(selcomps)){
        temp <- compiledDF[,grepl(selcomps[b], colnames(compiledDF))]
        temp <- temp[,grepl(paste("logFC", SigCol, sep = "|"), colnames(temp))]
        colnames(temp) <- c("logFC", "Significance")
        temp$SYMBOL <- rownames(temp)
        temp$Comparison <- selcomps[b]
        DFsub3 <- rbind(DFsub3, temp)
      }
      d1 <- DFsub3[abs(logFC) > FCCutoff & Significance < PCutoff,]
      d2 <- DFsub3[!(paste(DFsub3$SYMBOL, DFsub3$Comparison, sep = "") %in% paste(d1$SYMBOL, d1$Comparison, sep = "")),][, `:=` (logFC = NA, Significance = NA)][]
      d3 <- rbind(d1, d2)
      d3 <- d3[,c("logFC", "SYMBOL", "Comparison"), with = FALSE]
      DFsub <- reshape2::dcast(d3, SYMBOL~Comparison, value.var = "logFC")
      rownames(DFsub) <- DFsub$SYMBOL
      DFsub$SYMBOL <- NULL
      DFsub2 <- DFsub
      #### scale data ####
      if(ScaleData){
        transMat <- t(DFsub2)
        Rnames <- row.names(transMat)
        transMat <- suppressWarnings(apply(transMat, 2, as.numeric))
        rowScaled <- apply(as.matrix(transMat), 2, scale)
        rownames(rowScaled) <- Rnames
        DFsub2 <- as.data.frame(t(rowScaled))
      }
      DFsub2$names <- row.names(DFsub)
      mel <- as.data.table(reshape2::melt(DFsub2))
      if(ScaleData){
        setnames(mel, c("names", "variable", "value"), c("Gene", "Dataset", "Scaled log2(Fold Change)"))
        p <- ggplot(mel, aes(x = Dataset, y=Gene, fill = `Scaled log2(Fold Change)`)) +
          geom_tile(color = "white", lwd = 0.75, linetype = 1) +
          scale_fill_gradient2(low = "#075AFF",
                               high= "#FF0000") +
          coord_fixed() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
      } else {
        setnames(mel, c("names", "variable", "value"), c("Gene", "Dataset", "log2(Fold Change)"))
        p <- ggplot(mel, aes(x = Dataset, y=Gene, fill = `log2(Fold Change)`)) +
          geom_tile(color = "white", lwd = 0.75, linetype = 1) +
          scale_fill_gradient2(low = "#075AFF",
                               high= "#FF0000") +
          coord_fixed() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
      }
      return(p) }
    #### Obtain counts of selected genes and return a barplot ####
    if(plottype == "Bar"){
      selcomps <- unique(mgsub(colnames(compiledDF), c("logFC_", "AdjPValue_", "Pvalue_", "BHCorrection_"), c("","","","")))
      DFsub3 <- data.table()
      for(b in seq_along(selcomps)){
        temp <- compiledDF[,grepl(selcomps[b], colnames(compiledDF))]
        temp <- temp[,grepl(paste("logFC", SigCol, sep = "|"), colnames(temp))]
        colnames(temp) <- c("logFC", "Significance")
        temp$SYMBOL <- rownames(temp)
        temp$Comparison <- selcomps[b]
        DFsub3 <- rbind(DFsub3, temp)
      }
      great <- DFsub3[logFC > FCCutoff & Significance < PCutoff,]
      FCGreat <- great[, .(AveFC = mean(logFC, na.rm = TRUE)), by = "SYMBOL"]
      setnames(FCGreat, c("SYMBOL"), c("GeneName"))
      if(nrow(great) > 0){
        great <- data.table(table(great$SYMBOL)) %>% setnames(c("V1", "N"), c("GeneName", "Nup"))
        great <- merge(great, FCGreat, by = "GeneName")
        GreatCountAll <- great
        great <- great[order(great$Nup, decreasing = TRUE),]#[1:25,]
      }
      less <- DFsub3[logFC < -FCCutoff & Significance < PCutoff,]
      FCless <- less[, .(AveFC = mean(logFC, na.rm = TRUE)), by = "SYMBOL"]
      setnames(FCless, c("SYMBOL"), c("GeneName"))
      if(nrow(less) > 0){
        less <- data.table(table(less$SYMBOL)) %>% setnames(c("V1", "N"), c("GeneName", "Ndown"))
        less <- merge(less, FCless, by = "GeneName")
        lessCountAll <- less
        less <- less[order(less$Ndown, decreasing = TRUE),]#[1:25,]
      }
      if((nrow(great) > 0 & nrow(less) > 0)){
        final <- merge(great, less, by = "GeneName", all = TRUE)
        #### Merge Fold Change ####
        final$AveFC <- 0
        for(b in seq_len(nrow(final))){
          t <- c(final[b,]$AveFC.x, final[b,]$AveFC.y)
          t <- t[!is.na(t)]
          if(length(t) > 1){
            F2 <- data.table()
            for(c in seq_along(t)){
              F2 <- rbind(F2, final[b,])
            }
            F2$AveFC <- t
            final <- rbind(final, F2)
          } else{ final$AveFC[b] <- t[!is.na(t)]
          } }
        final$AveFC.x <- NULL; final$AveFC.y <- NULL
        final <- final[!AveFC == 0,]
        final <- unique(final)
        #### Merge counts ####
        final$direction <- "NA"
        final$counts <- 0
        for(b in seq_len(nrow(final))){
          up <- final$Nup[b]
          down <- final$Ndown[b]
          fc <- final$AveFC[b]
          if(fc > 0){
            final$direction[b] <- "up"
            final$counts[b] <- up
          } else {
            final$direction[b] <- "down"
            final$counts[b] <- down
          } }
        final$Nup <- NULL; final$Ndown <- NULL
        final <- rbind(final[direction == "up",][order(-counts, -AveFC),],
                       final[direction == "down",][order(counts, -AveFC),])
      }
      if(nrow(great) > 0 & nrow(less) == 0){
        final <- great
        final$direction <- "up"
        setnames(final, c("Nup"), c("counts"))
        final <- final[order(-counts, AveFC),]
      }
      if(nrow(great) == 0 & nrow(less) > 0){
        final <- less
        final$direction <- "down"
        setnames(final, c("Ndown"), c("counts"))
        final <- final[order(-counts, AveFC),]
      }
      if(!is.null(final)){
        final$GeneName <- factor(final$GeneName, levels = unique(final$GeneName ))
        p <- ggplot(final, aes(GeneName, AveFC, fill = direction)) +
          geom_bar(position = "dodge", stat = "identity") +
          geom_text(aes(label=counts), vjust=-0.5) +
          theme(axis.text.y = element_text(angle = 0, hjust = 1),
                axis.text.x = element_text(angle = 90, hjust = 1),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour="black") ) +
          ylab("Average log2(Fold change)")+
          xlab("Gene Name")+
          ggtitle(paste("average fold change across comparisons")) +
          scale_fill_manual(values = RColorBrewer::brewer.pal(11, "RdYlBu")[c(10,2)])
      }
      return(p) }
    if(plottype == "Table"){
      selcomps <- unique(mgsub(colnames(compiledDF), c("logFC_", "AdjPValue_", "Pvalue_", "BHCorrection_"), c("","","", "")))
      DFsub3 <- data.table()
      for(b in seq_along(selcomps)){
        temp <- compiledDF[,grepl(selcomps[b], colnames(compiledDF))]
        temp <- temp[,grepl(paste("logFC", SigCol, sep = "|"), colnames(temp))]
        colnames(temp) <- c("logFC", "Significance")
        temp$SYMBOL <- rownames(temp)
        temp$Comparison <- selcomps[b]
        DFsub3 <- rbind(DFsub3, temp)
      }
      d1 <- DFsub3[abs(logFC) > FCCutoff & Significance < PCutoff,]
      CompDF <- d1[!is.na(d1$logFC),]
      CompDF <- CompDF[order(SYMBOL),]
      return(as.data.frame(CompDF))
    } }
}

##########################################
#### Return metabolite DE information ####
##########################################
MetaboliteVolcno <- function(met, DataSelect, FCcut, Pcut, sigColSelect, Return){
  if(!DataSelect == "Select Dataset"){
    select <- as.data.table(met[[DataSelect]])
    tempDT <- select[,c("Metabolite", "Fold_Change", "log2FoldChange", "t_value", "pval", "padj")]
    tempDT[,Selected:= "Not selected"]
    tempDT[,Selected:=ifelse( ( abs(tempDT[["log2FoldChange"]]) > FCcut & tempDT[[sigColSelect]] < Pcut ), "Selected", "Not selected")]
    if(sigColSelect == "pval"){
      p <- ggplot(tempDT, aes(x=log2FoldChange,y=-log10(pval), color = Selected)) +
        xlab("log2(FC)")+ ylab(paste("-log10(", sigColSelect, ")", sep = "")) +
        ggtitle(names(DataSelect)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour="black")) +
        geom_point(data=subset(tempDT, Selected == "Not selected"),
                   aes(x=log2FoldChange,y=-log10(pval) ),
                   shape=21, size=2, colour="gray", fill = "gray") +
        geom_point(data=subset(tempDT, Selected == "Selected"),
                   aes(x=log2FoldChange,y=-log10(pval) ),
                   shape=21, size=3, colour="black", fill = "#E31A1C")
    }
    if(sigColSelect == "padj"){
      p <- ggplot(tempDT, aes(x=log2FoldChange,y=-log10(padj), color = Selected)) +
        xlab("log2(FC)")+ ylab(paste("-log10(", sigColSelect, ")", sep = "")) +
        ggtitle(names(DataSelect)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour="black")) +
        geom_point(data=subset(tempDT, Selected == "Not selected"),
                   aes(x=log2FoldChange,y=-log10(padj) ),
                   shape=21, size=2, colour="gray", fill = "gray") +
        geom_point(data=subset(tempDT, Selected == "Selected"),
                   aes(x=log2FoldChange,y=-log10(padj) ),
                   shape=21, size=3, colour="black", fill = "#E31A1C")
    }
    statement <- as.data.frame(list(c(
      paste("Total genes:", nrow(tempDT)),
      paste("Total Selected:", nrow(tempDT[Selected == "Selected",])),
      paste("Total Not Selected", nrow(tempDT[Selected == "Not selected",])),
      paste("Total Up-regulated Selected:", nrow(tempDT[Selected == "Selected" & log2FoldChange > 0,])),
      paste("Total Down-regulated Selected:", nrow(tempDT[Selected == "Selected" & log2FoldChange < 0,])))  ))
    colnames(statement) <- " "
    if(Return == "DT"){return(tempDT)}
    if(Return == "plot"){return(p)}
    if(Return == "Counts"){return(statement)}
  }
}

#################################
#### Load methylation tables ####
#################################
methylationDTload <- function(dataset){
  if(!dataset == "Select Dataset"){
    return(readRDS(file.path(homedir, "Methylation", paste(dataset, ".rds", sep = "" ) )))
  }
}

################################################################################
#### Perform Computations ######################################################
################################################################################
server <- function(input, output, session) {

  ##########################
  #### Overview Barplot ####
  ##########################
  output$DiseaseBar <- renderPlotly({
    SummaryPlot(over2=overview, features="Studies", attribute=input$attributes) })

  ########################
  #### Overview Table ####
  ########################
  output$OverviewTable <- DT::renderDT(
    overview[!duplicated(overview$ID),c("ID","Title","Year","GPLNumber","Profiling Resource","Tissue","Disease","Species","Website","DataType", "Description", "Resources")],
    filter = 'top', options = list(pageLength = 15, scrollX = TRUE, scrollY = '400px', autoWidth = TRUE, dom = 'ltipr'), escape = FALSE)

  ################
  #### DE Tab ####
  ################
  observeEvent(input$Dataset,{
    updateSelectInput(session,'Comparison',
                      choices=unique(AvailComps$Comparison[AvailComps$Dataset %in% input$Dataset] ))  })

  output$DEMessage <- reactive({
    req(input$submit)
    if(input$FC < 1){print("The fold change value you entered is less than 1. Please enter a fold change value greater than or equal to 1.")}  })

  output$volcano <- renderPlot({
    req(input$submit)
    VolcanoPlotR(DESE = DESE,
                 Dataset = isolate(input$Dataset),
                 Comparison = isolate(input$Comparison),
                 sigColSelect = isolate(input$HypothesisTestDE),
                 filterBy = isolate(input$Filter),
                 Gene = isolate(input$Genes),
                 FCcut=isolate(input$FC),
                 Pcut=isolate(input$Pval)) })

  output$text <- renderTable({
    req(input$submit)
    DEGcount(DESE = DESE,
             Dataset = isolate(input$Dataset),
             Comparison = isolate(input$Comparison),
             sigColSelect = isolate(input$HypothesisTestDE),
             filterBy = isolate(input$Filter),
             Gene = isolate(input$Genes),
             FCcut=isolate(input$FC),
             Pcut=isolate(input$Pval))  })

  output$DETable <- DT::renderDT(
    DETableDisplay(DESE = DESE, Dataset = input$Dataset, Comparison = input$Comparison
    ), filter = 'top', options = list(pageLength = 15, scrollX = TRUE, scrollY = '400px', autoWidth = TRUE, dom = 'ltipr'), escape = FALSE )

  ########################
  #### DE Heatmap tab ####
  ########################
  updateSelectizeInput(session, 'GenesHeat', choices = unique(genesAll), server = TRUE, selected = c("ACTR2","ALAS1","NFKB2","PSMB9","RFC2","TCF4","TRAF6","ABLIM1","ACSM3","ADD3","ADH5","ALDH9A1","ANXA1","ANXA3","AOX1","AP1B1","AR","ARF4","ATP2C1","ATP6AP2","BAIAP2","CCDC6","CCT8","CD48",
                                                                                                     "CHERP","COL4A1","COL6A3","CPOX","CREB1","CTBS","CYB5A","DAD1","DLGAP4","DOCK4","DST","EIF2B1","EMD","EPM2AIP1","FILIP1L","GGCX","GHR","HMG20B","HMGCR","HNRNPA2B1","IGF1","IMPA1","IQGAP1","ITGA9",
                                                                                                     "LPP","MSH6","NARS2","NENF","NOL7","PCK1","PCYT2","PFKFB1","PIGK","PLCB3","PMM1","POLRMT","PON3","PROS1","PSMB3","PTPN2","RAD52","RBM34","RLF","RRM1","SDHC","SELENBP1","SERPINB8","SLC12A2",
                                                                                                     "SLC4A4","SNX10","SP100","SREBF2","SUCLG1","TM7SF2","TMEM97","TNK2","TRIO","TSC22D2","TXN2","USP14","VPS41","XPA","ZBTB17","ZER1","AKAP13","DAB2","FMO1","AKR1A1","COX8A","EMP2","GSTM1","HINT1",
                                                                                                     "NDUFA2","PTPN18","ANKRD46","HMGCL","SLC37A4","CHSY1","ACOT2","ERCC5","HSPE1"))
  output$heatMessage <- reactive({if(input$FCMultiBar < 1){print("The fold change value you entered is less than 1. Please enter a fold change value greater than or equal to 1.")} })

  output$heat <- renderPlot({
    req(input$Heatsubmit)
    CrossDataHeat(DESE = DESE,
                  GeneSelection=isolate(input$GenesHeat),
                  Dataselection=isolate(input$DEdatasetHeat),
                  ScaleData = input$ScaleData,
                  plottype = "Heat",
                  FCCutoff = isolate(input$FCMultiBar),
                  PCutoff = isolate(input$PvalMultiBar),
                  SigCol = isolate(input$HypothesisTest) ) },  height = 1000, width=1300)

  output$heatBar <- renderPlot({
    req(input$Heatsubmit)
    CrossDataHeat(DESE = DESE,
                  GeneSelection=isolate(input$GenesHeat),
                  Dataselection=isolate(input$DEdatasetHeat),
                  ScaleData = input$ScaleData,
                  plottype = "Bar",
                  FCCutoff = isolate(input$FCMultiBar),
                  PCutoff = isolate(input$PvalMultiBar),
                  SigCol = isolate(input$HypothesisTest) ) },  height = 500, width=1200)

  output$heatText <- DT::renderDT(
    CrossDataHeat(DESE = DESE,
                  GeneSelection=input$GenesHeat,
                  Dataselection=input$DEdatasetHeat,
                  ScaleData = input$ScaleData,
                  plottype = "Table",
                  FCCutoff = isolate(input$FCMultiBar),
                  PCutoff = isolate(input$PvalMultiBar),
                  SigCol = isolate(input$HypothesisTest)
    ), filter = 'top', options = list(pageLength = 30, scrollX = TRUE, scrollY = '800px', autoWidth = TRUE, dom = 'ltipr'), escape = FALSE )

  ##########################################################
  ### Violin plot of multiple genes in a single project ####
  ##########################################################
  updateSelectizeInput(session, 'Gene_single', choices = unique(genesAll), server = TRUE, selected = c("GAPDH"))
  select_gene = function(){return(input$Gene_single)}
  anyna_single = reactive({check_na(raw_df, select_gene(), input$Dataset_single)})

  output$text_single <- renderText({
    if(anyna_single()){"Some of the selected gene(s) not available in the selected dataset"}
    else{"Here is the plot!"} })

  violin_single=reactive({box_violin_plot_v2.single(raw_df, select_gene(), input$Dataset_single, FALSE)})

  output$contents_single = renderPlot(violin_single())

  output$single_data_plot <- renderUI({
    plotOutput("contents_single", height = length(select_gene())*250 ) })

  output$download_single <- downloadHandler(
    filename = function() { paste("Violin_plot", input$Dataset_single, input$extension_single, sep = ".") },
    content = function(file){ ggsave(file, violin_single(), device = input$extension_single, width = 9, height = length(select_gene())*3) }
  )

  ##########################################
  ### Violin plot for multiple projects ####
  ##########################################
  updateSelectizeInput(session, 'Gene_multi', choices = unique(genesAll), server = TRUE, selected = c("GAPDH"))
  select_dis.multi = function(){return(input$Disease_multi)}
  select_tech.multi = function(){return(input$Tech_multi)}
  select_tissue.multi = function(){return(input$Tissue_multi)}
  select_df.multi = function(){return(overview%>%dplyr::filter(Disease%in%select_dis.multi() & Technology%in%select_tech.multi() & Tissue%in%select_tissue.multi())%>%dplyr::select(ID)%>%unlist%>%as.vector()%>%unique())}
  anyna_multi = reactive({check_na(raw_df, input$Gene_multi, select_df.multi())})
  output$text_multi <- renderText({
    if(anyna_multi()){"The selected gene not available in some of the selected dataset(s)"}
    else{"Here is the plot!"} })

  violin_multi=reactive({box_violin_plot_v2.multi(raw_df, overview, input$Gene_multi,  select_dis.multi(), select_tech.multi(),  select_tissue.multi(), FALSE)})
  output$contents_multi = renderPlot(violin_multi())

  output$multi_data_plot <- renderUI({ plotOutput("contents_multi", height = length(select_df.multi())*150) })

  output$download_multi <- downloadHandler(
    filename = function() { paste("Violin_plot", input$Gene_multi, input$extension_multi, sep = ".") },
    content = function(file){ ggsave(file, violin_multi(), device = input$extension_multi, width = 9, height = length(select_df.multi())*3) })

  ########################
  #### Proteomics tab ####
  ########################
  ProtData = reactive({
    req(input$ProteomicSubmit)
    ProteomicsDataLoad(Path=file.path(homedir, "Proteomic_3"), Fname=isolate(input$ProteomicDataset),
                       alpha = isolate(input$ProteomicPval), lfc = isolate(input$ProteomicFC),
                       sigCol=isolate(input$ProteomicHypothesisTestDE))   })

  output$Proteomicvolcano <- renderPlot({
    req(input$ProteomicSubmit)
    ProteomicVolcanoPlots(DEPList = ProtData(), addNames = TRUE, sigCol = isolate(input$ProteomicHypothesisTestDE))  })

  output$Proteomictext <- renderTable({
    req(input$ProteomicSubmit)
    ProteomicDEGcount(ProtDE = ProtData())  })

  output$ProteomicDETable <- DT::renderDT(
    ProteomicTableDisplay(Path=file.path(homedir, "Proteomic_3"), Fname=input$ProteomicDataset,
                          alpha = input$ProteomicPval, lfc = input$ProteomicFC,
                          sigCol=input$ProteomicHypothesisTestDE), filter = 'top', options = list(pageLength = 15, scrollX = TRUE, scrollY = '400px', autoWidth = TRUE, dom = 'ltipr'), escape = FALSE )

  ##################################
  #### Proteomic DE Heatmap tab ####
  ##################################
  updateSelectizeInput(session, 'ProteomicGenesHeat', choices = Proteins, server = TRUE, selected = c("EXOC6B","AKAP12","CTNNBIP1","VCPIP1","AP1AR","CNBP","KRT5","PLP2","RPS28", "ACTR2","ALAS1","NFKB2","PSMB9","RFC2","TCF4","TRAF6","ABLIM1","ACSM3","ADD3","ADH5","ALDH9A1","ANXA1","ANXA3","AOX1","AP1B1","AR","ARF4","ATP2C1","ATP6AP2","BAIAP2","CCDC6","CCT8","CD48",
                                                                                                      "CHERP","COL4A1","COL6A3","CPOX","CREB1","CTBS","CYB5A","DAD1","DLGAP4","DOCK4","DST","EIF2B1","EMD","EPM2AIP1","FILIP1L","GGCX","GHR","HMG20B","HMGCR","HNRNPA2B1","IGF1","IMPA1","IQGAP1","ITGA9",
                                                                                                      "LPP","MSH6","NARS2","NENF","NOL7","PCK1","PCYT2","PFKFB1","PIGK","PLCB3","PMM1","POLRMT","PON3","PROS1","PSMB3","PTPN2","RAD52","RBM34","RLF","RRM1","SDHC","SELENBP1","SERPINB8","SLC12A2",
                                                                                                      "SLC4A4","SNX10","SP100","SREBF2","SUCLG1","TM7SF2","TMEM97","TNK2","TRIO","TSC22D2","TXN2","USP14","VPS41","XPA","ZBTB17","ZER1","AKAP13","DAB2","FMO1","AKR1A1","COX8A","EMP2","GSTM1","HINT1",
                                                                                                      "NDUFA2","PTPN18","ANKRD46","HMGCL","SLC37A4","CHSY1","ACOT2","ERCC5","HSPE1"))

  output$Proteomicheat <- renderPlot({
    req(input$ProteomicHeatsubmit)
    CrossDataHeatProteomic(DESE = DEProtSE,
                           GeneSelection = input$ProteomicGenesHeat,
                           Dataselection = input$ProteomicDEdatasetHeat,
                           ScaleData = input$ProteomicScaleData,
                           plottype = "Heat",
                           FCCutoff = isolate(input$ProteomicFCMultiBar),
                           PCutoff = isolate(input$ProteomicPvalMultiBar),
                           SigCol = isolate(input$ProteomicHypothesisTest)
    )  },  height = 1000, width=1300)

  output$ProteomicheatBar <- renderPlot({
    req(input$ProteomicHeatsubmit)
    CrossDataHeatProteomic(DESE = DEProtSE,
                           GeneSelection = input$ProteomicGenesHeat,
                           Dataselection = input$ProteomicDEdatasetHeat,
                           ScaleData = input$ProteomicScaleData,
                           plottype = "Bar",
                           FCCutoff = isolate(input$ProteomicFCMultiBar),
                           PCutoff = isolate(input$ProteomicPvalMultiBar),
                           SigCol = isolate(input$ProteomicHypothesisTest)
    )  },  height = 1000, width=1300)

  output$ProteomicheatText <- DT::renderDT(
    CrossDataHeatProteomic(DESE = DEProtSE,
                           GeneSelection = input$ProteomicGenesHeat,
                           Dataselection = input$ProteomicDEdatasetHeat,
                           ScaleData = input$ProteomicScaleData,
                           plottype = "Table",
                           FCCutoff = isolate(input$ProteomicFCMultiBar),
                           PCutoff = isolate(input$ProteomicPvalMultiBar),
                           SigCol = isolate(input$ProteomicHypothesisTest)
    ), filter = 'top', options = list(pageLength = 30, scrollX = TRUE, scrollY = '800px', autoWidth = TRUE, dom = 'ltipr'), escape = FALSE )

  ########################
  #### Metabolics tab ####
  ########################
  output$Metabolomicsvolcano <- renderPlot({
    req(input$MetabolomicsSubmit)
    MetaboliteVolcno(met = metabolomicsData,
                     DataSelect = input$MetabolomicsDataset,
                     FCcut = input$MetabolomicsFC,
                     Pcut = input$MetabolomicsPval,
                     sigColSelect = input$MetabolomicsHypothesisTestDE,
                     Return = "plot") })

  output$Metabolomicstext <- renderTable({
    req(input$MetabolomicsSubmit)
    MetaboliteVolcno(met = metabolomicsData,
                     DataSelect = input$MetabolomicsDataset,
                     FCcut = input$MetabolomicsFC,
                     Pcut = input$MetabolomicsPval,
                     sigColSelect = input$MetabolomicsHypothesisTestDE,
                     Return = "Counts") })

  output$MetabolomicsDETable <- DT::renderDT({
    MetaboliteVolcno(met = metabolomicsData,
                     DataSelect = input$MetabolomicsDataset,
                     FCcut = input$MetabolomicsFC,
                     Pcut = input$MetabolomicsPval,
                     sigColSelect = input$MetabolomicsHypothesisTestDE,
                     Return = "DT") })

  #########################
  #### Methylation tab ####
  #########################
  output$MethylationDETable <- DT::renderDT({
    methylationDTload(dataset= input$MethylationDataset)
  })
}


