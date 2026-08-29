
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

#### loading tables ####
########################
#### Load genes ####
homedir <- "<path to data base>"
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
    tempDT$SYMBOL <- toupper(row.names(tempDT))
    tempDT <- as.data.table(tempDT)
    FDT <- tempDT[complete.cases(tempDT),]
    if(return == "Genes"){
      return(FDT$SYMBOL)
    } else { return(as.data.frame(FDT)) } } }

###############################
#### Volcano plot function ####
###############################
VolcanoPlotR <- function(DESE,Dataset,Comparison,sigColSelect,FCcut,Pcut,filterBy,Gene){
  if(FCcut >= 1){
    FCcut <- log2(FCcut)
    if(!Dataset == "Select Dataset"){
      Fname <- paste(Dataset, Comparison, sep = "_")
      assaySelect <- assays(DESE)[grepl(Fname, gsub("-", "_", names(assays(DESE))))]
      tempDT <- assaySelect[[1]]
      tempDT$SYMBOL <- toupper(row.names(tempDT))
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
                       aes(x=logFC,y=-log10(Pvalue) ),
                       shape=21, size=2, colour="gray", fill = "gray") +
            geom_point(data=subset(tempDT, SYMBOL == toupper(Gene)),
                       aes(x=logFC,y=-log10(Pvalue) ),
                       shape=21, size=3, colour="black", fill = "#E31A1C") }}
      return(p) } } }

############################
#### DEG count function ####
############################
DEGcount <- function(DESE = DESE, Dataset,Comparison,sigColSelect = "Pvalue",filterBy,Gene,FCcut,Pcut){
  if(FCcut >= 1){
    FCcut <- log2(FCcut)
    if(!Dataset == "Select Dataset"){
      Fname <- paste(Dataset, Comparison, sep = "_")
      assaySelect <- assays(DESE)[grepl(Fname, gsub("-", "_", names(assays(DESE)))   )  ]
      tempDT <- assaySelect[[1]]
      tempDT$SYMBOL <- row.names(tempDT)
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
  if(FCCutoff >= 1){
    FCCutoff <- log2(FCCutoff)
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
  samples_gse = paste(select_data, collapse = "|")
  print(samples_gse)
  df_subset = raw_df[rowData(raw_df)$SYMBOL%in%select_gene, grepl(samples_gse, colnames(raw_df))]
  print(df_subset)
  return(any(is.na(df_subset%>%assay())))
}

#############################################
#### Violin plot of FPKM Single data set ####
#############################################
box_violin_plot_v2.single = function(input_df, select_gene, select_proj, log10){
  # remove the NA ones, if no NA, ignore this step
  project_meta = colData(input_df)[,c("disease", "organism", "platform", "dataset", "rawcolumnnames", "cell_type", "tissue_type", "cell_line", "condition")]
  project_meta = project_meta[!is.na(project_meta$dataset),]
  raw_df_select_sample = project_meta%>%as.data.frame()%>%dplyr::filter(dataset%in%select_proj)%>%dplyr::select(rawcolumnnames)%>%unlist()%>%as.vector()
  # subset the data
  raw_df_summarized = input_df[rowData(input_df)$SYMBOL%in%select_gene, raw_df_select_sample]
  raw_df_melt = meltAssay(raw_df_summarized, assay_name="Expression", add_row_data = "SYMBOL", add_col_data = c("dataset", "rawcolumnnames", "condition"))%>%as.data.frame() # assay.type = "FPKM",
  # replace NAs with 0
  raw_df_melt = na.omit(raw_df_melt)
  # color palette
  box_count = length(unique(raw_df_melt$condition))
  coul = brewer.pal(9, "Set3")
  coul = colorRampPalette(coul)(box_count)
  # log transformation
  if(log10){
    raw_df_melt$Expression = log10(raw_df_melt$Expression+1)
    ylab_t = expression(~log[10]~(NormalizedExpression))
  } else{ ylab_t = "NormalizedExpression" }
  p=ggplot(data = raw_df_melt, aes(x=condition, y=Expression, fill=condition))+
    geom_violin(scale = "width", trim = T, alpha = 0.7)+
    geom_boxplot(outlier.shape = NA,coef = 0, fill="gray", width=0.3)+
    coord_flip()+
    ylab(ylab_t)+
    xlab("")+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(colour = 'black', size = 10),
          axis.text.x = element_text(colour = 'black', size = 10, vjust = 1, hjust=1),
          legend.position = "none",
          strip.text.y = element_text(angle = 0, size = 10))+
    stat_boxplot(geom='errorbar', linetype=1, width=0.2, position = "dodge2")+
    scale_fill_manual(values= coul) +
    facet_grid(SYMBOL~., scales = "free", space = "free", drop = F)+
    scale_x_discrete(drop = FALSE)
  print(p)
  return(p)
}

############################################
#### Violin plot of FPKM Multi data set ####
############################################
box_violin_plot_v2.multi = function(input_df, overview_file, select_gene, select_disease, select_tech, select_tissue, log10){
  # remove the NA ones, if no NA, ignore this step
  project_meta = colData(input_df)[,c("disease", "organism", "platform", "dataset", "rawcolumnnames", "cell_type", "tissue_type", "cell_line", "condition")]
  project_meta = project_meta[!is.na(project_meta$dataset),]
  select_group = overview_file%>%dplyr::filter(Disease%in%select_disease & Technology%in%select_tech & Tissue%in%select_tissue)%>%dplyr::select(ID)%>%unlist%>%as.vector()%>%unique()
  raw_df_select_sample = project_meta%>%as.data.frame()%>%dplyr::filter(dataset%in%select_group)%>%dplyr::select(rawcolumnnames)%>%unlist()%>%as.vector()
  raw_df_summarized = input_df[rowData(input_df)$SYMBOL%in%select_gene,raw_df_select_sample]
  raw_df_melt = meltAssay(raw_df_summarized, assay_name="Expression", add_row_data = "SYMBOL", add_col_data = c("dataset", "rawcolumnnames", "condition", "facet_name"))%>%as.data.frame()
  raw_df_melt = na.omit(raw_df_melt)
  # color palette
  box_count = length(unique(raw_df_melt$condition))
  coul = brewer.pal(9, "Set3")
  coul = colorRampPalette(coul)(box_count)
  # log transformation
  if(log10){
    raw_df_melt$Expression = log10(raw_df_melt$Expression+1)
    ylab_t = expression(~log[10]~(NormalizedExpression))
  } else{ ylab_t = "NormalizedExpression" }
  p=ggplot(data = raw_df_melt, aes(x=condition, y=Expression, fill=condition))+
    geom_violin(scale = "width", trim = T, alpha = 0.7)+
    geom_boxplot(outlier.shape = NA,coef = 0, fill="gray", width=0.3)+
    ylab(ylab_t)+
    xlab("")+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(colour = 'black', size = 10),
          axis.text.x = element_text(colour = 'black', size = 10),
          legend.position = "none",
          strip.text.y = element_text(angle = 0, size = 10))+
    stat_boxplot(geom='errorbar', linetype=1, width=0.2, position = "dodge2")+
    scale_fill_manual(values= coul) +
    facet_wrap(facet_name~., scales = "free", dir = "v", ncol= 1, strip.position="right", drop = F)
  return(p)
}

#################################
#### Load proteomics dataset ####
#################################
ProteomicsDataLoad <- function(Path, Fname, alpha, lfc, sigCol){
  lfc <- log2(lfc)
  #### Load proteomics summarized experiment object from DB
  SelectedList <- list()
  SelectedList[[1]] <- readRDS(file.path(Path, paste(Fname, ".rds", sep = "")))
  names(SelectedList) <- gsub(".rds", "", Fname)
  #### Proteomics Denote DE samples ####
  dataDiff <- SelectedList
  dep <- list()
  for(i in seq_along(SelectedList)){ dep[[i]] <- add_rejections2(SelectedList[[i]], alpha = alpha, lfc = lfc, sigCol=sigCol) }
  names(dep) <- names(SelectedList)
  return(dep)
}
add_rejections2 <- function (diff=dataDiff[[i]], alpha = 0.05, lfc = 1, sigCol="p.adj")  {
  if (is.integer(alpha))
    alpha <- as.numeric(alpha)
  if (is.integer(lfc))
    lfc <- as.numeric(lfc)
  assertthat::assert_that(inherits(diff, "SummarizedExperiment"),
                          is.numeric(alpha), length(alpha) == 1, is.numeric(lfc),
                          length(lfc) == 1)
  row_data <- rowData(diff, use.names = FALSE) %>% as.data.frame()
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(diff)), "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE) }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(diff)), "'\nRun test_diff() to obtain the required columns",
         call. = FALSE) }
  cols_p <- grep(sigCol, colnames(row_data))
  cols_diff <- grep("_diff", colnames(row_data))
  if (length(cols_p) == 1) {
    rowData(diff)$significant <- row_data[, cols_p] <= alpha & abs(row_data[, cols_diff]) >= lfc
    rowData(diff)$contrast_significant <- rowData(diff, use.names = FALSE)$significant
    colnames(rowData(diff))[ncol(rowData(diff, use.names = FALSE))] <- gsub(sigCol, "significant", colnames(row_data)[cols_p]) }
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
                           sign_df, by = "name") }
  return(diff)
}

##########################################
#### Proteomics Volcano plot function ####
##########################################
VolcanoPlot <- function(DEPList = DEPList, addNames = TRUE, sigCol){
  if(sigCol == "p.adj"){adj <- TRUE; BHadj <- FALSE
  } else if(sigCol == "p.val"){ adj <- FALSE; BHadj <- FALSE
  } else if(sigCol == "BHCorrection"){ adj <- FALSE; BHadj <- TRUE}
  contrasts <- colnames(rowData(DEPList[[1]]))[!colnames(rowData(DEPList[[1]])) %in% c("name","Protein.IDs","Gene.names","Description","Numberofpeptides","ID","significant")]
  contrasts <- unique(mgsub(contrasts, c("_CI.L", "_CI.R", "_diff", "_p.adj", "_p.val", "_BHCorrection", "_significant"), c("", "", "", "", "", "", "")))
  contrasts <- contrasts[grepl("_vs_", contrasts)]
  PlotList <- list()
  for(i in seq_along(DEPList)){
    comps <- contrasts
    for(b in seq_along(comps)){
      tempLit <- list()
      nam <- paste(names(DEPList)[i], comps[b], sep = "-")
      tempLit[[1]] <- plot_volcano2(DEPList[[i]], contrast = comps[b], label_size = 3, add_names = addNames, adjusted = adj, BHadjusted = BHadj) + ggtitle(nam)
      names(tempLit) <- nam
      PlotList <- c(PlotList, tempLit) } }
  return(PlotList)
}
plot_volcano2 <- function (dep, contrast, label_size = 3, add_names = TRUE, adjusted = FALSE, plot = TRUE, BHadjusted) {
  if (is.integer(label_size)){
    label_size <- as.numeric(label_size)
    assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                            is.character(contrast), length(contrast) == 1, is.numeric(label_size),
                            length(label_size) == 1, is.logical(add_names), length(add_names) ==
                              1, is.logical(adjusted), length(adjusted) == 1, is.logical(plot),
                            length(plot) == 1) }
  row_data <- rowData(dep, use.names = FALSE)
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '", deparse(substitute(dep)), "'.\nRun make_unique() to obtain required columns."), call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '", deparse(substitute(dep)), "'.\nRun test_diff() to obtain the required columns."), call. = FALSE)
  }
  if (length(grep("_significant", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_significant' columns are not present in '", deparse(substitute(dep)), "'.\nRun add_rejections() to obtain the required columns."), call. = FALSE)
  }
  if (length(grep(paste(contrast, "_diff", sep = ""), colnames(row_data))) == 0) {
    valid_cntrsts <- row_data %>% data.frame() %>% select(ends_with("_diff")) %>%
      colnames(.) %>% gsub("_diff", "", .)
    valid_cntrsts_msg <- paste0("Valid contrasts are: '", paste0(valid_cntrsts, collapse = "', '"), "'")
    stop("Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n", valid_cntrsts_msg, call. = FALSE)
  }
  diff <- grep(paste(contrast, "_diff", sep = ""), colnames(row_data))
  if (adjusted) { p_values <- grep(paste(contrast, "_p.adj", sep = ""), colnames(row_data))
  } else if(BHadjusted){ p_values <- grep(paste(contrast, "_BHCorrection", sep = ""), colnames(row_data))
  } else { p_values <- grep(paste(contrast, "_p.val", sep = ""), colnames(row_data)) }
  signif <- grep(paste(contrast, "_significant", sep = ""), colnames(row_data))
  df <- data.frame(x = row_data[, diff], y = -log10(row_data[,p_values]), significant = row_data[, signif], name = row_data$name) %>%
    filter(!is.na(significant)) %>% arrange(significant)
  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)
  p <- ggplot(df, aes(x, y)) + geom_vline(xintercept = 0) +
    geom_point(aes(col = significant)) + geom_text(data = data.frame(), aes(x = c(Inf, -Inf), y = c(-Inf, -Inf), hjust = c(1, 0), vjust = c(-1, -1), label = c(name1, name2), size = 5,
                                                                            fontface = "bold")) + labs(title = contrast, x = expression(log[2] ~ "Fold change")) + theme_DEP1() + theme(legend.position = "none") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "grey"))
  if (add_names) { p <- p + ggrepel::geom_text_repel(data = filter(df, significant), aes(label = name), size = label_size, box.padding = unit(0.1, "lines"), point.padding = unit(0.1, "lines"), segment.size = 0.5) }
  if (adjusted) { p <- p + labs(y = expression(-log[10] ~ "Adjusted p-value"))
  } else if(BHadjusted){  p <- p + labs(y = expression(-log[10] ~ "BH Adjusted p-value"))
  } else { p <- p + labs(y = expression(-log[10] ~ "P-value")) }
  if (plot) { return(p)
  } else {
    df <- df %>% select(name, x, y, significant) %>% arrange(desc(x))
    colnames(df)[c(1, 2, 3)] <- c("protein", "log2_fold_change", "p_value_-log10")
    if (adjusted) { colnames(df)[3] <- "adjusted_p_value_-log10" }
    if (BHadjusted) { colnames(df)[3] <- "BH_adjusted_p_value_-log10" }
    return(df)
  }
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

  violin_single=reactive({box_violin_plot_v2.single(raw_df, select_gene(), input$Dataset_single, input$log_single)})

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

  violin_multi=reactive({box_violin_plot_v2.multi(raw_df, overview, input$Gene_multi,  select_dis.multi(), select_tech.multi(),  select_tissue.multi(), input$log_multi)})
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
    VolcanoPlot(DEPList = ProtData(), addNames = TRUE, sigCol = isolate(input$ProteomicHypothesisTestDE))  })

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


