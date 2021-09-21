#   ____________________________________________________________________________
#   Server       R 3.6.1                                                           ####

library(shiny)
library(plotly)
library(dplyr)
library(tidyr)
library(shinyjs)
library(DT)
library(rgl)
library(shinyBS)


shinyServer(function(input, output, session) {
    ##################################################################################
    #################################### Updated Day #################################
    ##################################################################################
    output$footerModifiedDay <- renderText(paste0("Updated ", modifiedDay))
    output$aboutModifiedDay <- renderText(paste0("Last modified at ", modifiedDay))
    ##################################################################################
    ####################################### Home #####################################
    ##################################################################################
    ## Tissues
    species_tissues <- rmgd %>% count(Species, Tissue) %>% select(-n) %>% count(Species)
    output$totalTissues <- renderText(paste0("We currently track ", length(unique(rmgd %>% count(Species, Tissue) %>% select(Tissue))$Tissue), " tissues ..."))
    output$speciesTissuesStat <- renderPlotly({
        tissues_count_fig <- plot_ly(species_tissues, x = ~Species, y = ~n, color = ~Species, type = "bar") %>% 
            layout(autosize = F, 
                   yaxis = list(title = 'Number'),
                   xaxis = list(tickangle = 45, title = NA),
                   plot_bgcolor  = "rgba(242, 252, 255, 0)",
                   paper_bgcolor = "rgba(242, 252, 255, 0)") %>%
            config(displaylogo = FALSE)
    })
    
    ## cell types
    cellTypes <- rmgd %>% count(Species, CellType) %>% select(-n) %>% count(Species)
    output$totalCellTypes <- renderText(paste0("... in over ", length(unique(rmgd %>% drop_na() %>% select(CellType))$CellType), " cell types ..."))
    output$cellTypesStat <- renderPlotly({
        cell_count_fig <- plot_ly(cellTypes, x = ~Species, y = ~n, color = ~Species, type = "bar") %>% 
            layout(autosize = F, 
                   yaxis = list(title = 'Number'),
                   xaxis = list(tickangle = 45, title = NA),
                   plot_bgcolor  = "rgba(242, 252, 255, 0)",
                   paper_bgcolor = "rgba(242, 252, 255, 0)") %>%
            config(displaylogo = FALSE)
    })
    
    ## Marker Genes
    tissues_count <- rmgd %>% count(Species, Tissue, xlab)
    output$totalMarkerGenes <- renderText(paste0("... in over ", sum(tissues_count$n), " marker genes that have evidenced via RNA in situ hybridization and GFP reporter"))
    
    output$markerGenesStat <- renderPlotly({
        marker_genes_count_fig <- plot_ly(tissues_count, x = ~xlab, y = ~n, color = ~Tissue, type = "bar", height = 800) %>% 
            layout(autosize = T, 
                   yaxis = list(title = 'Number'),
                   xaxis = list(tickangle = 45, title = NA),
                   plot_bgcolor  = "rgba(242, 252, 255, 0)",
                   paper_bgcolor = "rgba(242, 252, 255, 0)") %>%
            config(displaylogo = FALSE)
    })
    ##################################################################################
    ##################################### Protocol ###################################
    ##################################################################################
    output$troubleshooting <- renderTable(troubleshooting)
    output$enzymolysis <- renderTable(enzymolysisSolution)
    ##################################################################################
    ##################################### Workflow ###################################
    ##################################################################################
    # Reactive value for selected parameters ----
    generateSingleCellWorkflow <- reactive({
        scw <- readLines("F:/Dash/PlantSingleCellDataBase/www/data/SingleCellWorkflow.R")
        ## Loading Cell Ranger count output
        scwd <- gsub(pattern = "ProjectName", replace = input$projectname, x = scw)
        scwd <- gsub(pattern = "SamplesMetadata", replace = input$samplesmetadata, x = scwd)
        scwd <- gsub(pattern = "GenesAnnotationFile", replace = input$genesannotation, x = scwd)
        scwd <- gsub(pattern = "MitochondriaGenes", replace = input$mitochondria, x = scwd)
        scwd <- gsub(pattern = "ChloroplastGenes", replace = input$chloroplast, x = scwd)
        scwd <- gsub(pattern = "MarkerGenesOfKnownCellType", replace = input$markergenes, x = scwd)
        scwd <- gsub(pattern = "FilteredFeatureBcMatrix", replace = input$filtered_feature_matrix_path, x = scwd)
        ## Quality Control
        scwd <- gsub(pattern = "QCMADs", replace = input$mads, x = scwd)
        raw_nCount_RNA <- unlist(strsplit(as.character(c(paste(input$rawnCount, collapse = ","))),","))
        scwd <- gsub(pattern = "RawCountRNAMin", replace = raw_nCount_RNA[1], x = scwd)
        scwd <- gsub(pattern = "RawCountRNAMax", replace = raw_nCount_RNA[2], x = scwd)
        raw_nFeature_RNA <- unlist(strsplit(as.character(c(paste(input$rawnFeature, collapse = ","))),","))
        scwd <- gsub(pattern = "RawFeatureRNAMin", replace = raw_nFeature_RNA[1], x = scwd)
        scwd <- gsub(pattern = "RawFeatureRNAMax", replace = raw_nFeature_RNA[2], x = scwd)
        scwd <- gsub(pattern = "Log10GenesPerUMIValue", replace = input$log10Genes, x = scwd)
        scwd <- gsub(pattern = "MitoRatioValue", replace = input$mitoRatio, x = scwd)
        scwd <- gsub(pattern = "GeneExpressiomAtLeastCell", replace = input$genecell, x = scwd)
        scwd <- gsub(pattern = "RemoveGenesThatExpressionValusEqualtoZeroInAllCell", replace = input$removegenes, x = scwd)
        ## Normalization and ScaleData
        scwd <- gsub(pattern = "NormalizationMethod", replace = input$method, x = scwd)
        scwd <- gsub(pattern = "NumberOfMostVariantGenes", replace = input$maxnfeatures, x = scwd)
        ## Dimensionality Reduction
        clusteringalgorithm <- as.character(c(paste(input$clusteringalgorithm, collapse = ",")))
        scwd <- gsub(pattern = "ClusteringAlgorithm", replace = clusteringalgorithm, x = scwd)
        scwd <- gsub(pattern = "PCAMaxdims", replace = input$PCAmaxdims, x = scwd)
        scwd <- gsub(pattern = "tSNEPerplexitys", replace = input$tSNEperplexitys, x = scwd)
        scwd <- gsub(pattern = "UMAPMindists", replace = input$UMAPmindists, x = scwd)
        ## Clustering
        resolutions <- as.character(c(paste(input$resolutions, collapse = ",")))
        scwd <- gsub(pattern = "ResolutionsForCluster", replace = resolutions, x = scwd)
        ## Marker Genes Identification
        scwd <- gsub(pattern = "ExpectedNumberOfCellClusters", replace = input$numberClusters, x = scwd)
        scwd <- gsub(pattern = "TopGeneToPlotDoHeatmapFeaturePlotVlnPlotAndDotPlot", replace = input$topGenes, x = scwd)
        scwd <- gsub(pattern = "OnlyPosMarkers", replace = input$postype, x = scwd)
        scwd <- gsub(pattern = "MinPctMarkers", replace = input$minpct, x = scwd)
        scwd <- gsub(pattern = "ThreshUseMarkers", replace = input$threshuse, x = scwd)
        scwd <- gsub(pattern = "LogfcThresholdMarkers", replace = input$logfc, x = scwd)
        return(scwd)
    })
    
    ## Downloadable R script of selected dataset ----
    output$gRscript <- downloadHandler(
        filename = function() {
            paste0("F:/Dash/PlantSingleCellDataBase/www/data/", input$projectname, "_SingleCellWorkflow.R")
        },
        content = function(file) {
            writeLines(generateSingleCellWorkflow(), con = file)
        }
    )
    ##################################################################################
    ################################### MarkerGeneDB #################################
    ##################################################################################
    output$markerGenesTableData <- renderDataTable({
        rmgddf$GeneIdUrl <- with(rmgddf, 
                             ifelse(grepl("^AT", rmgddf$GeneId), paste0("<a href=\"https://www.arabidopsis.org/servlets/TairObject?name=", rmgddf$GeneId, "&type=locus\" target=\"_blank\">", rmgddf$GeneId, "</a>"),
                             ifelse(grepl("^LOC_Os", rmgddf$GeneId), paste0("<a href=\"http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=", rmgddf$GeneId, "\" target=\"_blank\">", rmgddf$GeneId, "</a>"),
                             ifelse(grepl("^Os", rmgddf$GeneId), paste0("<a href=\"https://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1?name=", rmgddf$GeneId, "\" target=\"_blank\">", rmgddf$GeneId, "</a>"), 
                             ifelse(grepl("^Zm", rmgddf$GeneId), paste0("<a href=\"https://ensembl.gramene.org/Zea_mays/Gene/Summary?g=", rmgddf$GeneId, "\" target=\"_blank\">", rmgddf$GeneId, "</a>"),
                             ifelse(grepl("^GRMZM", rmgddf$GeneId), paste0("<a href=\"https://maizegdb.org/gene_center/gene/", rmgddf$GeneId, "\" target=\"_blank\">", rmgddf$GeneId, "</a>"),
                             ifelse(grepl("^Solyc", rmgddf$GeneId), paste0("<a href=\"http://223.31.159.9/tomato2/getGene.php?trans_fac_id=", rmgddf$GeneId, "\" target=\"_blank\">", rmgddf$GeneId, "</a>"),
                             ifelse(grepl("^EF", rmgddf$GeneId), paste0("<a href=\"https://www.maizegdb.org/gene_center/gene/", rmgddf$GeneId, "\" target=\"_blank\">", rmgddf$GeneId, "</a>"),
                             ifelse(grepl("^AC", rmgddf$GeneId), paste0("<a href=\"https://www.maizegdb.org/gene_center/gene/", rmgddf$GeneId, "\" target=\"_blank\">", rmgddf$GeneId, "</a>"),
                             rmgddf$GeneId)))))))))
        
        rmgddftb <- rmgddf %>% select(Species, Tissue, CellType, SubcellularType, GeneIdUrl, RNAsituHybridization, ReferenceUrl, RefTitle) %>% 
                    mutate("RNA in situ hybridization || GFP reporter" = paste0("<img src=\"", RNAsituHybridization, "\" width=\"20%\" />")) %>% 
                    mutate(Reference = paste0("<a href=\"", ReferenceUrl, "\" target=\"_blank\">", RefTitle, "</a>")) %>% 
                    select(-RNAsituHybridization, -ReferenceUrl, -RefTitle) %>% 
                    rename("Cell Type" = CellType, "Subcellular Type" = SubcellularType, "Gene Id" = GeneIdUrl)

    }, options =list(pageLength=10, columnDefs = list(list(className = 'dt-center', targets="_all"))), escape = FALSE)
    
    output$markerGenesTableDataUnR <- renderText({ RJSONIO::toJSON(rmgddf, byrow=T, colNames=T) })
    ##################################################################################
    ############################## Single Cell Literatures ###########################
    ##################################################################################
    library(tidyverse)
    library(RISmed)
    library(openxlsx)
    library(lubridate)
    library(readxl)
    
    generateSingleCellLiteratures <- reactive({
        input$Search
        purpose <- "search"
        #key_words <- c('rice roots resolution')
        key_words <- c(paste0("scRNA single cell ", input$searchbox))
        print(key_words)
        date_range_min <- 2019
        latest <- year(ymd(Sys.Date()))
        date_range_max <- latest
        
        withProgress(message = 'Search Single Cell Literatures ...', {
            paperData <- EUtilsSummary(query = key_words, db = "pubmed", datetype = "pdat", mindate = date_range_min, maxdate = date_range_max, retmax = 100)
        
            journal <- ISOAbbreviation(EUtilsGet(paperData))
            authors <- Author(EUtilsGet(paperData))
            titles <- ArticleTitle(EUtilsGet(paperData))
            abstract <- AbstractText(EUtilsGet(paperData))
            #affiliations <- Affiliation(EUtilsGet(paperData))
            year_pub <- YearPubmed(EUtilsGet(paperData))
            article_id <- ArticleId(EUtilsGet(paperData))
        })
        
        paperdf <- data.frame(article_id, titles, abstract, journal, year_pub) %>%
                              mutate(PMID = paste0("<a href=\"https://pubmed.ncbi.nlm.nih.gov/", article_id, "\">", article_id, "</a>"))
        
        return(paperdf)
    })
    
    output$paperTable <- renderDataTable({
        paperdf <- generateSingleCellLiteratures()
        paperdfTable <- paperdf %>% rename(Title = titles, Abstract = abstract, Journal = journal, Year = year_pub) %>%
                                    select(PMID, Title, Journal, Year, Abstract)
    }, options =list(pageLength=10, columnDefs = list(list(className = 'dt-center', targets="_all"))), escape = FALSE)
    
    ## Text mining and word cloud
    ## http://www.sthda.com/english/wiki/text-mining-and-word-cloud-fundamentals-in-r-5-simple-steps-you-should-know
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### word cloud
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    library("tm")
    library("SnowballC")
    #library("wordcloud2")
    library("wordcloud")
    library("RColorBrewer")
    
    generateSingleCellLiteraturesDTM <- reactive({
        paperdf <- generateSingleCellLiteratures()
        docs <- Corpus(VectorSource(paste0(paperdf$abstract, paperdf$titles)))
        toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
        docs <- tm_map(docs, toSpace, "/")
        docs <- tm_map(docs, toSpace, "@")
        docs <- tm_map(docs, toSpace, "\\|")
        # Convert the text to lower case
        docs <- tm_map(docs, content_transformer(tolower))
        # Remove numbers
        docs <- tm_map(docs, removeNumbers)
        # Remove english common stopwords
        docs <- tm_map(docs, removeWords, stopwords("english"))
        # Remove your own stop word
        # specify your stopwords as a character vector
        docs <- tm_map(docs, removeWords, c("the", "and", "but", "single", "cell", "plant", "use", "true", "thus", "can", "gene", "type")) 
        # Remove punctuations
        docs <- tm_map(docs, removePunctuation)
        # Eliminate extra white spaces
        docs <- tm_map(docs, stripWhitespace)
        # Text stemming
        docs <- tm_map(docs, stemDocument)
        dtm <- TermDocumentMatrix(docs)
        return(dtm)
    })
    
    generateSingleCellLiteraturesDF <- reactive({
        dtm <- generateSingleCellLiteraturesDTM()
        m <- as.matrix(dtm)
        v <- sort(rowSums(m),decreasing=TRUE)
        docsdf <- data.frame(word = names(v),freq=v)
        return(docsdf)
    })
    
    #output$workCloud <- renderWordcloud2({
    output$workCloud <- renderPlot({
        docsdf <- generateSingleCellLiteraturesDF()
        set.seed(1234)
        par(mar = rep(0, 4))
        wordcloud(words = docsdf$word, freq = docsdf$freq, min.freq = 1, max.words=200, 
                  random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))
        #wordcloud2(docsdf, backgroundColor = "#F2FCFF", shape = 'circle')
    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### Plot word frequencies
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    generateSingleCellLiteraturesDFterms <- reactive({
        input$update
        
        dtm <- generateSingleCellLiteraturesDTM()
        m <- as.matrix(dtm)
        v <- sort(rowSums(m),decreasing=TRUE)
        docsdfTerms <- data.frame(word = names(v),freq=v)
        
        findFreqTerms(dtm, lowfreq = 4)
        findAssocs(dtm, terms = input$searchbox, corlimit = 0.3)
        return(docsdfTerms)
    })
    
    output$workFrequencies <- renderPlotly({
        docsdfTerms <- generateSingleCellLiteraturesDFterms()
        
        plot_ly(docsdfTerms, x = ~docsdfTerms[1:20,]$word, y = ~docsdfTerms[1:20,]$freq, color = ~docsdfTerms[1:20,]$word, type = "bar") %>% 
                layout(autosize = F,
                       width = 1140,
                       showlegend = FALSE,
                       yaxis = list(title = 'Number'),
                       xaxis = list(tickangle = 45, title = NA),
                       plot_bgcolor  = "rgba(242, 252, 255, 0)",
                       paper_bgcolor = "rgba(242, 252, 255, 0)") %>%
                config(displaylogo = FALSE)
    })
    ##################################################################################
    ###################################### Submit ####################################
    ##################################################################################
    observe({
        if (is.null(input$attachmentsPath)) return()
        file.copy(from = input$attachmentsPath$datapath, to = paste0("F:/Dash/PlantSingleCellDataBase/www/data/Submit/", input$attachmentsPath$name), recursive = TRUE)
        print(input$attachmentsPath$name)
    })
    
    output$uploadFiles <- renderText({
        HTML(paste0("The upload files have save to <a href=\"data/Submit/", input$attachmentsPath$name,"\" >",input$attachmentsPath$name,"</a>."))
    })
    ##################################################################################
    ###################################### About #####################################
    ##################################################################################
    output$contactUs <- renderText({
        HTML(paste0("Just <i class=\"far fa-envelope\"></i> <a href=\"mailto:", adminEmail, "?subject=这是邮件的主题&body=这是邮件的内容\" rel=\"nofollow\">send Email (", adminEmail, ")</a> to us."))
    })
    ##################################################################################
    ####################################### END ######################################
    ##################################################################################
})