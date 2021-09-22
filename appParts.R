#   ____________________________________________________________________________
#   Protocol                                                               ####

protocolMenu <- function() {
    tagList(
        div(class = "container",
            h1("Dissociate Protoplast Protocol", class = "title fit-h1"),
            h2("Introduction"),
            fluidRow(column(12, img(src='images/ProtoplastProtocol.png', align = "center", width = "100%"))),
            hr(),
            p(),
            hr(),
            h2("Materials"),
            hr(),
            h3("Biological materials"),
            p("Young plant."),
            h3("Reagents"),
            tags$ol("MES (Sigma, cat. no. M8250)"),
            tags$ol("Mannitol (Sigma, cat.no.M4125)"),
            tags$ol("CaCl2 (Sigma, cat.no. C7902)"),
            tags$ol("KCl (Sigma, cat.no. P3911)"),
            tags$ol("MgCl2 (Sigma, cat.no. M9272)"),
            tags$ol("NaCl (BBI, cat.no.7647-14-5)"),
            tags$ol("CellulaseR10 (YakultPharmaceuticalInd. Co., Ltd., Japan)"),
            tags$ol("MacerozymeR10 (YakultPharmaceuticalInd. Co., Ltd., Japan)"),
            tags$ol("Bovine Serum Albumin (BSA) (Cat. No. #A8020, Solarbio, Beijing, China)"),
            tags$ol("HemiCellulase (Sigma)"),
            h3("Equipment"),
            tags$ol("Cell Strainer (40 Micrometre, cat.no.F613461-9001)"),
            tags$ol("Scalpel"),
            tags$ol("Petri dish"),
            tags$ol("10mL or 50mL centrifuge tube and 2mL centrifuge tube"),
            tags$ol("Centrifuge with adjustable acceleration"),
            h3("Reagent setup"),
            h4("Enzymolysis Solution"),
            fluidRow(column(12, tableOutput('enzymolysis'))),
            p("Water bath at 55 degree celsius for 10 minutes, then add ddH2O (0.2ml), CaCl2 (100 Micrometre) and BSA (0.01g) in sequence."),
            hr(),
            h2("Procedure"),
            hr(),
            h3("Preparation of protoplast suspension"),
            h4("1. Enzymatic hydrolysis of tissue samples"),
            p("Plant tissue were cut into 1-2 mm cuttings. Pour the previously prepared enzymolysis solution into the petri dish, add the cut hypocotyls to the enzymolysis solution (ensure that the tissue is completely immersed in the enzymolysis solution), shake at room temperature at a low speed, and avoid light for 4 -6 hours."),
            h4("2. Filtration suspension"),
            p("The single cell suspension after enzymolysis was filtered with a 40 Micrometre Cell Strainer. Take the filtered cell suspension into a centrifuge tube, add twice the volume of PBS (with 0.04% BSA) and mix it in an ultra-low speed centrifuge, centrifuge at 70g for 3 minutes, discard the supernatant, repeat this step, and finally leave a certain amount of precipitate, Mix well with a cut pipette tip, filter into a new centrifuge tube with a 40μm cell strainer."),
            h4("3. Microscopy for cell state"),
            p("A small amount of suspension was dripped onto the slide and observe the cell number and cell viability (tested by Fluorescein Diacetate) under a microscope."),
            h4("4. Reverse transcription and 10X genomics library construction"),
            p("Approximately 15,000 counted cells were used for reverse transcription and subsequently library preparation according to instruction."),
            hr(),
            h2("Troubleshooting"),
            hr(),
            fluidRow(column(12, tableOutput('troubleshooting'))),
            hr(),
            h2("Anticipated results"),
            hr(),
            img(src='images/ProtoplastProtocol.jpg', align = "center", width = "50%")
            )
        )

}
        
#   ____________________________________________________________________________
#   Workflow                                                               ####

workflowMenu <- function() {
    tagList(
        div(class = "container",
            h1("Workflow For Analysis plants scRNA Data", class = "title fit-h1"),
            h4("Customize your own workflow for the analysis of plant's scRNA by customizing", strong("Cell Ranger Count Information"),", ",strong("Quality Control threshold"),", ",strong("Normalization and ScaleData Method"),", ",strong("Dimensionality Reduction Algorithm"),", ",strong("Clustering Resolutions")," and ", strong("Marker Genes Identification threshold"),"."),
            h4(icon("link", lib = "font-awesome"), tags$a(onclick="$('li:eq(11) a').tab('show');", "Read more in About")),
            fluidRow(column(12, img(src='images/Workflow.png', align = "center", width = "100%"))),
            hr(),
            h2("Loading Cell Ranger count output", class = "h2"),
            p("Output directory from ",code("Cell Ranger count"), " and metadata for this scRNA-seq dataset. It's must conform to fixed file format."),
            fluidRow(
                column(4, textInput("projectname", label = "Project Name:", placeholder = "Cotton")),
                column(4, textInput("samplesmetadata", label = "Samples Metadata:", placeholder = "metadata")),
                column(4, textInput("genesannotation", label = "Genes Annotation File:", placeholder = "annotation"))
            ),
            fluidRow(
                column(4, textInput("mitochondria", label = "Mitochondria Genes:", placeholder = "mitochondria")),
                column(4, textInput("chloroplast", label = "Chloroplast Genes:", placeholder = "chloroplast")),
                column(4, textInput("markergenes", label = "Marker Genes Of Known Cell-Type:", placeholder = "marker"))
            ),
            fluidRow(
                column(12, textInput("filtered_feature_matrix_path", label = "filtered_feature_bc_matrix:", placeholder = "Matrix Path", width = '95%'))
            ),
            hr(),
            h2("Quality Control", class = "h2"),
            p("Set follows parameters to removing low-quality cells based on the median absolute deviation (",strong("MADs"),") from the median value of each metric across all cells, number of UMIs per cell (",strong("raw_nCount_RNA"),"), number of genes detected per cell (",strong("raw_nFeature_RNA"),"), number of genes per UMI for each cell (",strong("log10GenesPerUMI"),"), mitochondrial counts ratio (",strong("mitoRatio"),") and genes which are expressed in cells (",strong("Gene Expressiom At Least Cell"),")."),
            fluidRow(
                column(4, sliderInput("mads", label = "MADs", min=1, max=10, value=3, step=1)),
                column(4, sliderInput("rawnCount", label = "raw_nCount_RNA", min = 100, max = 50000, value = c(1000, 20000), step=100)),
                column(4, sliderInput("rawnFeature", label = "raw_nFeature_RNA", min=500, max=5000, value = c(500, 4100), step=100))
            ),
            fluidRow(
                column(4, sliderInput("log10Genes", label = "log10GenesPerUMI", min=0.1, max=1, value=0.8, step=0.1)),
                column(4, sliderInput("mitoRatio", label = "mitoRatio", min = 0.1, max = 1, value = 0.2, step=0.1)),
                column(4, sliderInput("genecell", label = "Gene Expressiom At Least Cell", min=1, max=100, value = 10, step=1))
            ),
            fluidRow(
                column(12, radioButtons("removegenes", label = "Remove Genes That Expression Valus Equal to Zero In All Cell", c("Yes", "No"), inline = TRUE))
            ),
            hr(),
            h2("Normalization and ScaleData", class = "h2"),
            p("Read literature ",icon("link", lib = "font-awesome"), tags$a(href="http://dx.doi.org/10.20517/jtgg.2020.48", target="_blank", "(Schneider et al., 2021)"), " to assess different clustering solutions of different normalization parameters."),
            fluidRow(
                column(6, selectInput("method", label = "Normalization Method:", c("CLR" = "CLR", "LogNormalize" = "LogNormalize"), width = '95%')),
                column(6, sliderInput("maxnfeatures", label = "Number Of Most Variant Genes", min=1000, max=5000, value=3000, step=100, width = '95%'))
            ),
            hr(),
            h2("Clustering", class = "h2"),
            p("The resolution is an important argument that sets the granularity of the downstream clustering and will need to be optimized for every individual experiment."),
            p("For datasets of 3,000 - 5,000 cells, the resolution set between 0.4-1.4 generally yields good clustering."),
            p("Increased resolution values lead to a greater number of clusters, which is often required for larger datasets."),
            fluidRow(
                column(12, 
                       selectInput("resolutions", label = "Resolutions For Cluster:",
                                   multiple = TRUE,
                                   c("0.2" = "0.2",
                                     "0.4" = "0.4",
                                     "0.6" = "0.6",
                                     "0.8" = "0.8",
                                     "1.0" = "1.0",
                                     "1.4" = "1.4"), width = '95%'))
            ),
            hr(),
            h2("Dimensionality Reduction", class = "h2"),
            p("Selected uses PCA, tSNE and UMAP to plot of cells in each cluster."),
            fluidRow(
                column(12,
                       checkboxGroupInput("clusteringalgorithm", 
                                          label = "clustering algorithm:",
                                          c("PCA" = "PCA", "tSNE" = "tSNE", "UMAP" = "UMAP"),
                                          selected = c("PCA","tSNE","UMAP"), inline = TRUE, width = '95%'))
            ),
            fluidRow(
                column(4, sliderInput("PCAmaxdims", label = "PCA maxdims", min=10, max=100, value=30, step=5)),
                column(4, sliderInput("tSNEperplexitys", label = "tSNE perplexitys", min=30, max=100, value=60, step=10)),
                column(4, sliderInput("UMAPmindists", label = "UMAP mindists", min=0.01, max=1, value=0.2, step=0.1))
            ),
            tags$blockquote(strong("tSNE perplexitys:"), " how to balance attention between local and global aspects of your data. See more in ",  icon("link", lib = "font-awesome"), tags$a(href="https://distill.pub/2016/misread-tsne/", target="_blank", "How to Use t-SNE Effectively"), "."),
            tags$blockquote(strong("UMAP mindists:"), " this controls how tightly the embedding is allowed compress points together. Larger values ensure embedded points are moreevenly distributed, while smaller values allow the algorithm to optimise more accurately with regard to local structure."),
            hr(),
            h2("Marker Genes Identification", class = "h2"),
            p("Identification of all markers for each cluster and annotation the cell types of the different clusters via our ", icon("link", lib = "font-awesome"), tags$a(onclick="$('li:eq(5) a').tab('show');", target="_blank", "MarkerGeneDB (Reviewed)"),"."),
            fluidRow(
                column(4, sliderInput("numberClusters", label = "Expected Number Of Cell Clusters:", min=1, max=100, value=15, step=1)),
                column(6, sliderInput("topGenes", label = "Top Gene To Plot DoHeatmap, FeaturePlot, VlnPlot and DotPlot:", min=1, max=20, value=10, step=1, width = '95%')),
                column(2, radioButtons("postype", label = "only.pos", c("Yes", "No"), inline = TRUE))
            ),
            fluidRow(
                column(4, sliderInput("minpct", label = "min.pct", min=0.1, max=1, value=0.25, step=0.05)),
                column(4, sliderInput("threshuse", label = "thresh.use", min = 0.1, max = 1, value = 0.25, step=0.05)),
                column(4, sliderInput("logfc", label = "logfc.threshold", min=0.1, max=1, value = 0.25, step=0.05))
            ),
            hr(),
            h2("Trajectory inference", class = "h2"),
            p("We are suggesting using ", icon("link", lib = "font-awesome"), tags$a(href="http://cole-trapnell-lab.github.io/monocle-release/docs/", target="_blank", "Monocle")," and ", icon("link", lib = "font-awesome"), tags$a(href="https://dynverse.org/", target="_blank", " dyno")," package to performed this 
trajectory inference for whole cluster or subcluster."),
            hr(),
            h2("Download", class = "h2"),
            p(icon("hand-o-right"), "Celebrate! Your custom workflow is almost finished. Just click the button below to download it and used to analysis your scRNA-seq data."),
            br(),
            fluidRow(
                column(4, br()),
                column(4, downloadButton('gRscript',label = "Generate R Analysis Script", class="butt"),
                       br(), 
                       tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}"))),
                column(4, br()))
        )
    )
}

#   ____________________________________________________________________________
#   SingleCellDB                                                            ####

singleCellDBMenu <- function() {
    tagList(
        div(class = "container",
            h1("Plant Single Cell Sequening Raw Data", class = "title fit-h1"),
            p("Plant Single Cell Database (SingleCellDB), is a publicly available repository of high throughput sequencing data of plant single cell. The archive accepts data from all plant."),
            p("Now (", modifiedDay, "), we have collected",length(unlist(strsplit(unique(rmgddf$DataAvailability), "_"))),"raw sequencing data of ", length(unique(rmgddf$Species)), " species."),
            h2(tags$em("Arabidopsis thaliana"), class = "h2"),
            fluidRow(
                column(6, wellPanel(
                    h2("Root", class = "h2"),
                    p("Expression profiling by high throughput sequencing."),
                    fluidRow(
                        column(6, downloadButton('atr1',"GSE158761",class="butt", href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158761"),
                            br(), 
                            tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}"))),
                        column(6, downloadButton('atr2',"GSE123818",class="butt", href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123818"),
                            br(), 
                            tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}")))
                ))),
                column(6, wellPanel(
                    h2("Flower", class = "h2"),
                    p("Expression profiling by high throughput sequencing."),
                    fluidRow(
                        column(6, downloadButton('atf1',"E-MTAB-9174",class="butt", href = "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9174/"),
                               br(), 
                               tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}")))
                )))
            ),
            fluidRow(
                column(6, wellPanel(
                    h2("Leaf", class = "h2"),
                    p("Expression profiling by high throughput sequencing."),
                    fluidRow(
                        column(6, downloadButton('atl1',"GSE161482",class="butt", href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161482"),
                               br(), 
                               tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}")))
                    ))),
                column(6, wellPanel(
                    h2("Female Gametophytes", class = "h2"),
                    p("Expression profiling by high throughput sequencing."),
                    fluidRow(
                        column(6, downloadButton('atfg1',"SRP160651",class="butt", href = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP160651"),
                               br(), 
                               tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}")))
                    )))
            ),
            hr(),
            h2(tags$em("Zea mays"), class = "h2"),
            fluidRow(
                column(6, wellPanel(
                    h2("Leaf", class = "h2"),
                    p("Expression profiling by high throughput sequencing."),
                    fluidRow(
                        column(6, downloadButton('zml1',"GSE157759",class="butt", href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157759"),
                               br(), 
                               tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}")))
                    ))),
                column(6, wellPanel(
                    h2("Staminate Primordia", class = "h2"),
                    p("Expression profiling by high throughput sequencing."),
                    fluidRow(
                        column(6, downloadButton('zmsp1',"GSE155178",class="butt", href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155178"),
                               br(), 
                               tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}")))
                    )))
            ),
            fluidRow(
                column(6, wellPanel(
                    h2("Stem", class = "h2"),
                    p("Expression profiling by high throughput sequencing."),
                    fluidRow(
                        column(6, downloadButton('zmsam1',"SRP266308",class="butt", href = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP266308"),
                               br(), 
                               tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}")))
                    ))),
                column(6, wellPanel(
                    h2("Developing Maize Ears", class = "h2"),
                    p("Expression profiling by high throughput sequencing."),
                    fluidRow(
                        column(4, downloadButton('zmdme1',"SRP272723",class="butt", href = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP272723"),
                               br(), 
                               tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}"))),
                        column(4, downloadButton('zmdme2',"SRP272726",class="butt", href = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP272726"),
                               br(), 
                               tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}"))),
                        column(4, downloadButton('zmdme3',"SRP272727",class="butt", href = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP272727"),
                               br(), 
                               tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}")))
                    )))
            ),
            hr(),
            h2(tags$em("Oryza sativa"), class = "h2"),
            fluidRow(
                column(12, wellPanel(
                    h2("Root", class = "h2"),
                    p("Expression profiling by high throughput sequencing."),
                    fluidRow(
                        column(4, downloadButton('osr1',"GSE146035",class="butt", href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146035"),
                               br(), 
                               tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}"))),
                        column(4, downloadButton('osr2',"SRP308960",class="butt", href = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP308960"),
                               br(), 
                               tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}"))),
                        column(4, downloadButton('osr3',"SRP309176",class="butt", href = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP309176"),
                               br(), 
                               tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}")))
                    )))
            ),
            hr(),
            h2(tags$em("Solanum lycopersicum"), class = "h2"),
            fluidRow(
                column(12, wellPanel(
                    h2("Shoot Apex", class = "h2"),
                    p("Expression profiling by high throughput sequencing."),
                    fluidRow(
                        column(6, downloadButton('slsa1',"SRP280242",class="butt", href = "https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP280242"),
                               br(), 
                               tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}")))
                    )))
            ),
            hr(),
        )
    )
}

#   ____________________________________________________________________________
#   SingleCellLiteratures                                                   ####

singleCellLiteratures <- function() {
    tagList(
        div(class = "container",
            h1("Plant Single Cell Literatures", class = "title fit-h1"),
            p("You can do searches of plant single cell literature using keywords by submitting your query in the search box ..."),
            p("A text mining was performed that allow us to highlight the most important keywords in searched plant scRNA literature and illustrated them as a word cloud."),
            fluidRow(column(12, textInput("searchbox", "", "arabidopsis", width = '100%'))),
            fluidRow(column(5, br()),
                     column(2, submitButton("Search", icon("search"))),
                     column(5, br())
            ),
            hr(),
            h2("Summary of Relevant Single Cell Literatures", class = "h2"),
            fluidRow(column(12, dataTableOutput('paperTable'))),
            hr(),
            h2("Statistical Graph of Keywords", class = "h2"),
            #fluidRow(column(12, wordcloud2Output(outputId = "workCloud", height = 550, width = "100%"))),
            fluidRow(column(12, plotOutput(outputId = "workCloud", height = 500, width = "100%"))),
            hr(),
            h2("Explore Frequent Terms and Associations", class = "h2"),
            fluidRow(column(12, plotlyOutput(outputId = "workFrequencies", width = 'auto'))),
            br()
        )
    )
}

#   ____________________________________________________________________________
#   Submit                                                                  ####
# https://stackoverflow.com/questions/20857068/sending-email-from-shiny
# https://github.com/senthilthyagarajan/shinyemail/blob/master/app.R
# https://www.geeksforgeeks.org/how-to-send-an-email-from-javascript/

submitMenu <- function() {
    tagList(
        div(class = "container",
            h1("Submit Data To Our Database or Feedback Question With Us", class = "title fit-h1"),
            pageWithSidebar(
                titlePanel(""),
                sidebarPanel(
                    textInput("fromEmail", "From:", value = "from@gmail.com"),
                    selectInput("emailType", "Contact Us For", c("Submit Data" = "SubmitData", "Faceback Question" = "FacebackQuestion")),
                    textInput("emailSubject", "Subject:", placeholder = "A subject"),
                    actionButton("sendEmail", icon = icon("envelope", lib = "font-awesome"), "Send mail"),
                    a(actionButton(inputId = "email1", label = "Contact Admin", 
                                   icon = icon("envelope", lib = "font-awesome")),
                      href="mailto:my_awesome_email_address.com")
                ),
                mainPanel(aceEditor("emailMessage", placeholder = "write message here",  height = "374px"))
            )
        )
    )
}

#   ____________________________________________________________________________
#   Download                                                                ####

downloadMenu <- function() {
    tagList(
        div(class = "container",
            h1("Download Data", class = "title fit-h1"),
            h2("Plant scRNA-seq Analysis Workflow Environment", class = "h2"),
            fluidRow(
                column(6, wellPanel(
                    h2("Environment", class = "h2"),
                    p(code("SingleCellCondaEnvironment.yml"), " file used to create a single cell analysis environment and ", code("RsessionInfo.txt"), " used to check the R packages version."),
                    fluidRow(
                        column(6, downloadButton('dn1',"SingleCellCondaEnvironment.yml",class="butt", href = "data/SingleCellCondaEnvironment.yml"),
                            br(), 
                            tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}"))),
                        column(6, downloadButton('dn2',"RsessionInfo.txt",class="butt", href = "data/RsessionInfo.txt"),
                            br(), 
                            tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}")))
                ))),
                column(6, wellPanel(
                    h2("Test Data", class = "h2"),
                    p("Include: Samples Metadata, Genes Annotation File, Mitochondria Genes, Chloroplast Genes, Marker Genes Of Known Cell-Type and filtered_feature_bc_matrix."),
                    fluidRow(
                        column(6, downloadButton('dn3',"TestData.tar.gz",class="butt", href = "data/TestData.tar.gz"),
                            br(), 
                            tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}")))
                )))
            ),
            hr(),
            h2("Reviewed Marker Genes", class = "h2"),
            fluidRow(
                column(12, wellPanel(
                    p("All markers have been evidenced via RNA in situ hybridization and GFP reporter. Besides, all results were verified manually."),
                    fluidRow(
                        column(12, downloadButton('dn4',"ReviewedMarkerGeneDatabase.csv",class="butt", href = "data/ReviewedMarkerGeneDatabase.csv"),
                            br(), 
                            tags$head(tags$style(".butt{background-color:#230682;} .butt{color: #e6ebef;}")))
                )))
            )
        )
    )
}
