
###################################
#### MultiOmics User Interface ####
###################################
library(shiny)
library(plotly)
library(shinyWidgets)

#### loading tables ####
########################
homedir <- "<path to data base>"
#### Load genes ####
DESE <- readRDS(file.path(homedir, "ProcessFiles", "SumarizedExp_DB.rds"))
AvailComps <- data.table(Dataset = gsub("_.+", "", names(assays(DESE))), Comparison = gsub("^.+?_", "", names(assays(DESE))) )
# Raw RNA expression data
raw_df = readRDS(file.path(homedir, "ProcessFiles", "expression_norm.v2.RDS"))
genesAll = rowData(raw_df)$SYMBOL
dataAll = unique(colData(raw_df)$dataset)%>%na.omit()
overviewfile = fread(file.path(homedir, "OverviewFiles", "GEODataOverview3.csv"))%>%as.data.frame()
diseaseAll = unique(overviewfile$Disease)
techAll = unique(overviewfile$Technology)
tissueAll = unique(overviewfile$Tissue)
#### Proteomic data ####
ProteomicAvail <- gsub(".rds", "", list.files(file.path(homedir, "Proteomic_3")))
DEProtSE <- readRDS(file.path(homedir, "ProcessFiles", "SumarizedProtExp_DB.rds"))
#### Metabolomics data ####
metabolomicsData <- readRDS(file.path(homedir, "Metabolomics", "mtbls298.de.RDS"))
names(metabolomicsData) <- c("Basal_Artery_Basal_Vein", "Insulin_Artery_Insulin_Vein")
#### Methylation data ####
Metfiles <- list.files(file.path(homedir, "Methylation"))
Metfiles <- gsub(".rds", "", Metfiles[grepl(".rds", Metfiles)])

#### App Interface ####
#######################
ui <- fluidPage(
  tabsetPanel(
    tabPanel("Homepage", fluid = TRUE,
             sidebarPanel(markdown('<br> <br>'),
                          selectizeInput("attributes", "Select attribute:",choices = c("Tissue","Disease","Species","Technology", "Year"), selected = "Disease"),
                          markdown('<br> <br>') ),
             mainPanel(plotlyOutput("DiseaseBar")), #,
             DT::dataTableOutput("OverviewTable") ),

    tabPanel("Transcriptomics", fluid = TRUE,
             tabsetPanel(
               tabPanel("Expression Exploration - Single-Dataset", fluid = TRUE,
                        titlePanel("Expression Exploration"),
                        fluidRow(
                          column(3, wellPanel(selectizeInput("Dataset_single", "Select Dataset:", choices = dataAll, selected = "GSE102485", multiple = FALSE),
                                              selectizeInput("Gene_single", "Select genes of interest:", choices = NULL, multiple = TRUE, width = '100%',  size = 6),
                                              materialSwitch(inputId = "log_single",label = "Log scale"),
                                              actionButton("Expressionsubmit", "Submit")) ),
                          column(9, verbatimTextOutput("text_single"),
                                 uiOutput('single_data_plot'),
                                 wellPanel( radioButtons("extension_single", "Save As:", choices = c("png", "pdf", "svg"), inline = TRUE),
                                            downloadButton("download_single", "Save Plot") ) ))),

               tabPanel("Expression Exploration - Multi-Dataset", fluid = TRUE,
                        titlePanel("Expression Exploration"),
                        fluidRow(
                          column(3, wellPanel( selectizeInput("Disease_multi", "Select Disease:", choices = diseaseAll, selected = "Heart Failure", multiple = TRUE ),
                                               selectizeInput("Tissue_multi", "Select Tissue:", choices = tissueAll, selected = "Heart", multiple = TRUE ),
                                               selectizeInput("Tech_multi", "Select Technology:", choices = c("Array", "RNAseq"), selected = techAll, multiple = TRUE ),
                                               selectizeInput("Gene_multi", "Select gene of interest:", choices = NULL, multiple = FALSE, width = '100%',  size = 6),
                                               materialSwitch(inputId = "log_multi",label = "Log scale"),
                                               actionButton("Expressionsubmit", "Submit") )  ),
                          column(9, verbatimTextOutput("text_multi"),
                                 uiOutput('multi_data_plot'),
                                 wellPanel( radioButtons("extension_multi", "Save As:", choices = c("png", "pdf", "svg"), inline = TRUE),
                                            downloadButton("download_multi", "Save Plot") )))   ),

               tabPanel("Differential Expression Analysis", fluid = TRUE,
                        sidebarPanel(
                          selectizeInput("Dataset", "Select Dataset:", choices = c(unique(AvailComps$Dataset), "Select Dataset"), selected = "Select Dataset"),
                          selectizeInput("Comparison", "Select Comparison:", ""),
                          selectInput("Filter", "Highlight points by:", choices = c("Significance", "Gene")),
                          conditionalPanel(condition = "input.Filter == 'Significance'",
                                           numericInput("FC", "Enter log2(Fold Change) cutoff", 1, min = 0, max = 70),
                                           selectInput("HypothesisTestDE", "Hypothesis Test", choices = c("Pvalue", "AdjPValue"), multiple = FALSE, selected = "AdjPValue"),
                                           numericInput("Pval", "Enter p-value cutoff", 0.05, min = 0, max = 1),
                                           uiOutput("text") ),
                          conditionalPanel(condition = "input.Filter == 'Gene'",
                                           textInput("Genes", "Enter gene to highlight")),#,
                          actionButton("submit", "Plot Selection") ),
                        mainPanel(markdown(' <br> <br> <br>'),
                                  plotOutput("volcano") ),
                        fluidRow(verbatimTextOutput("DEMessage"),  align="center"),
                        DT::dataTableOutput("DETable") ),

               tabPanel("Differntial expression - Multi-Dataset", fluid = TRUE,
                        markdown('##### **Overview**
                                 - Tools in this section can be used to explore the level of gene expression across a custom selection of genes and datasets. <br>
                                 - Genes can be selected by typing in the first few letters of the gene name in the `Select genes` menu and then selecting the gene of interest from the drop down menu. Alternatively,
                                 A gene can be removed my clicking on its name and pressing the backspace button. Several genes are loaded at startup as an example. <br>
                                 - Labels under the `Select Multiple Datasets` menu follow a naming convention with underscore separaters in which the first value indicates the
                                 dataset number, the second value indicates the treatment sample, and the third value indicates the reference sample used in a differential
                                 gene expression analysis (e.g.: Dataset_Treatment_Reference).
                                 - The first dataset can be selected by clicking on a dataset of interest under the `Select Multiple Datasets` menu. Additional
                                 datasets can be selected by holding the control key and left-clicking on the additional datasets of interest. Selected comparisons are
                                 highlighted in grey. <br>
                                 - Clicking on the `Plot Selection` button will generate a heatmap of the fold changes for each selected gene across selected differential gene expression
                                 comparisons. In addition, a barplot of the average up- and down- regulated expression level of each selected gene is plotted as a bar plot below
                                 the heatmap. This barplot is ordered in descending order by the number of selected comparisons a gene is up-regulated in. The numbers
                                 above the bars indicate the number of comparisons the average expression level was calculated from. <br>
                                 - Significance and fold change cutoff values can be applied. First, select the hypothesis test to use for the significance cutoff.
                                 Next, enter the Fold Change cutoff and Significance cutoff to apply. The absolute value of the fold change is applied to the cutoff.
                                 Clicking the Plot Selection button will adjust the bargraph to only include comparison values for all genes that meet the respective cutoff values.
                                 '),
                        fluidRow(column(6, selectizeInput("GenesHeat", "Select genes:", choices = NULL, multiple = TRUE, width = '100%',  size = 6)),
                                 column(6, selectInput("DEdatasetHeat", "Select Multiple Datasets:",choices = names(assays(DESE)), multiple = TRUE, selected = NULL,width = '100%', selectize = FALSE, size = 6 ) ) ),
                        fluidRow(column(2, selectInput("ScaleData", "Row Scale Data:", choices = c(TRUE, FALSE), multiple = FALSE, selected = FALSE) ),
                                 column(2, selectInput("HypothesisTest", "Hypothesis Test", choices = c("Pvalue", "AdjPValue"), multiple = FALSE, selected = "AdjPValue") ),
                                 column(3, numericInput("FCMultiBar", "Fold Change cutoff", 1, min = 1, max = 70)   ),
                                 column(2, numericInput("PvalMultiBar", "Significance cutoff", 1, min = 0, max = 1) ),
                                 fluidRow(column(1, align="center", actionButton("Heatsubmit", "Plot Selection") )) ),
                        fluidRow(verbatimTextOutput("heatMessage"),  align="center"),
                        fluidRow(column(12, plotOutput("heat", inline = TRUE)),  align="center"), br(), br(),
                        fluidRow(column(12, plotOutput("heatBar", inline = TRUE)),  align="center"), br(), br(),
                        fluidRow(column(12, DT::dataTableOutput("heatText")),  align="center" ) ),
             ) ,width = 10),

    tabPanel("Proteomics", fluid = TRUE,
             tabsetPanel(
               tabPanel("Differential Expression Analysis", fluid = TRUE,
                        sidebarPanel(
                          selectInput("ProteomicDataset", "Select Dataset:", choices = c(ProteomicAvail, "Select Dataset"), selected = "Select Dataset"),
                          numericInput("ProteomicFC", "Enter log2(Fold Change) cutoff", 1, min = 0, max = 70),
                          selectInput("ProteomicHypothesisTestDE", "Hypothesis Test", choices = c("p.adj", "p.val", "BHCorrection"), multiple = FALSE, selected = "p.adj"),
                          numericInput("ProteomicPval", "Enter significance cutoff", 0.05, min = 0, max = 1),
                          actionButton("ProteomicSubmit", "Plot Selection"),
                          uiOutput("Proteomictext")  ),
                        mainPanel(markdown(' <br> <br> <br>'),
                                  plotOutput("Proteomicvolcano")  ),
                        DT::dataTableOutput("ProteomicDETable") ),

               tabPanel("Differntial expression - Multi-Dataset", fluid = TRUE,
                        markdown('##### **Overview**
                                 - Tools in this section can be used to explore the level of gene expression across a custom selection of genes and datasets. <br>
                                 - Genes can be selected by typing in the first few letters of the gene name in the `Select genes` menu and then selecting the gene of interest from the drop down menu. Alternatively,
                                 A gene can be removed my clicking on its name and pressing the backspace button. Several genes are loaded at startup as an example. <br>
                                 - Labels under the `Select Multiple Datasets` menu follow a naming convention with underscore separaters in which the first value indicates the
                                 dataset number, the second value indicates the treatment sample, and the third value indicates the reference sample used in a differential
                                 gene expression analysis (e.g.: Dataset_Treatment_Reference).
                                 - The first dataset can be selected by clicking on a dataset of interest under the `Select Multiple Datasets` menu. Additional
                                 datasets can be selected by holding the control key and left-clicking on the additional datasets of interest. Selected comparisons are
                                 highlighted in grey. <br>
                                 - Clicking on the `Plot Selection` button will generate a heatmap of the fold changes for each selected gene across selected differential gene expression
                                 comparisons. In addition, a barplot of the average up- and down- regulated expression level of each selected gene is plotted as a bar plot below
                                 the heatmap. This barplot is ordered in descending order by the number of selected comparisons a gene is up-regulated in. The numbers
                                 above the bars indicate the number of comparisons the average expression level was calculated from. <br>
                                 - Significance and fold change cutoff values can be applied. First, select the hypothesis test to use for the significance cutoff.
                                 Next, enter the Fold Change cutoff and Significance cutoff to apply. The absolute value of the fold change is applied to the cutoff.
                                 Clicking the Plot Selection button will adjust the bargraph to only include comparison values for all genes that meet the respective cutoff values.
                                 '),
                        fluidRow(column(6, selectizeInput("ProteomicGenesHeat", "Select genes:", choices = NULL, multiple = TRUE, width = '100%',  size = 6)),
                                 column(6, selectInput("ProteomicDEdatasetHeat", "Select Multiple Datasets:",choices = names(assays(DEProtSE)), multiple = TRUE, selected = NULL,width = '100%', selectize = FALSE, size = 6 ) ) ),
                        fluidRow(column(2, selectInput("ProteomicScaleData", "Row Scale Data:", choices = c(TRUE, FALSE), multiple = FALSE, selected = FALSE) ),
                                 column(2, selectInput("ProteomicHypothesisTest", "Hypothesis Test", choices = c("Pvalue", "AdjPValue", "BHCorrection"), multiple = FALSE, selected = "AdjPValue") ),
                                 column(3, numericInput("ProteomicFCMultiBar", "log2 Fold Change cutoff", 0, min = 0, max = 70)   ),
                                 column(2, numericInput("ProteomicPvalMultiBar", "Significance cutoff", 1, min = 0, max = 1) ),
                                 fluidRow(column(1, align="center", actionButton("ProteomicHeatsubmit", "Plot Selection") )) ),
                        fluidRow(column(12, plotOutput("Proteomicheat", inline = TRUE)),  align="center"), br(), br(),
                        fluidRow(column(12, plotOutput("ProteomicheatBar", inline = TRUE)),  align="center"), br(), br(),
                        fluidRow(column(12, DT::dataTableOutput("ProteomicheatText")),  align="center" )
               ),
             ) ,width = 10),

    tabPanel("Methylation", fluid = TRUE,
             sidebarPanel(
               selectInput("MethylationDataset", "Select Dataset:", choices = c(Metfiles, "Select Dataset"), selected = "Select Dataset")
             ),
             DT::dataTableOutput("MethylationDETable") ),

    tabPanel("Metabolomics", fluid = TRUE,
             sidebarPanel(
               selectInput("MetabolomicsDataset", "Select Comparison:", choices = c(names(metabolomicsData), "Select Dataset"), selected = "Select Dataset"),
               numericInput("MetabolomicsFC", "Enter log2(Fold Change) cutoff", 1, min = 0, max = 70),
               selectInput("MetabolomicsHypothesisTestDE", "Hypothesis Test", choices = c("pval", "padj"), multiple = FALSE, selected = "pval"),
               numericInput("MetabolomicsPval", "Enter significance cutoff", 0.05, min = 0, max = 1),
               actionButton("MetabolomicsSubmit", "Plot Selection"),
               uiOutput("Metabolomicstext")  ),
             mainPanel(markdown(' <br> <br> <br>'),
                       plotOutput("Metabolomicsvolcano")  ),
             DT::dataTableOutput("MetabolomicsDETable")
    ),
  ) )

