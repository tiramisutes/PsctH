#   ____________________________________________________________________________
#   UI                                                                      ####

library(shiny)
library(plotly)
library(shinyjs)
library(shinyBS)
library(wordcloud2)

source("appParts.R")
source("readData.R")


shinyUI(fluidPage(
             navbarPage(title = div(tags$a(img(src="images/logo_horizontal.png", height=35), onclick="$('li:eq(0) a').tab('show');"),
                                    style = "position: relative; top: -5px;"), # Navigation bar
                        windowTitle = "Plant Single Cell Hub", #title for browser tab
                   theme = "style/style.css",
                   id = 'menus',
                   footer = htmlTemplate("html/footer.html"),
                   fluid = TRUE,
                   position = "fixed-top",
                   collapsible = TRUE,
                   
                   # ----------------------------------
                   # tab panel 1 - Home
                   tabPanel("Home", icon=icon("home"),
                            htmlTemplate("html/home.html"),
                            tags$script(src = "plugins/scripts.js"),
                            tags$head(
                              tags$link(rel = "stylesheet", 
                                        type = "text/css", 
                                        href = "plugins/font-awesome-4.7.0/css/font-awesome.min.css"),
                              tags$link(rel = "icon", 
                                        type = "image/png", 
                                        href = "images/logo_icon.png"),
                              tags$style(type="text/css", 
                                         "body {padding-top: 70px;}")
                            )
                   ),
                   
                   # ----------------------------------
                   # tab panel 2 - Workflow
                   navbarMenu("Guides", icon=icon("tasks"),
                              tabPanel("Dissociate Protoplast Protocol", protocolMenu()),
                              tabPanel("scRNA Analysis Pipelines", workflowMenu())
                   ),
                   
                   # ----------------------------------
                   # tab panel 3 - MarkerGeneDB
                   navbarMenu("MarkerGeneDB", icon=icon("database"),
                              tabPanel("MarkerGeneDB (Reviewed)", htmlTemplate("html/MarkerGeneDB.html"), icon=icon("battery-full")),
                              tabPanel("MarkerGeneDB (Unreviewed)", htmlTemplate("html/MarkerGeneDBUnR.html"), icon=icon("battery-half"))
                   ),
                   
                   # ----------------------------------
                   # tab panel 4 - SingleCellDB
                   tabPanel("SingleCellDB", icon=icon("university"),
                            singleCellDBMenu()
                   ),
                   
                   # ----------------------------------
                   # tab panel 5 - Single Cell Literatures
                   tabPanel("SingleCellLiteratures", icon=icon("graduation-cap"),
                            singleCellLiteratures()
                   ),
                   
                   # ----------------------------------
                   # tab panel 6 - Submit
                   #tabPanel("Submit", icon=icon("upload"), submitMenu()),
                   
                   # https://www.geeksforgeeks.org/how-to-send-an-email-from-javascript/
                   tabPanel("Submit", htmlTemplate("html/email.html", filedd = fluidRow(column(12, fileInput("attachmentsPath", "Attachments Path", multiple = FALSE, width = "100%")))), icon=icon("upload"),
                            tags$head(
                              tags$link(rel = "stylesheet",
                                        type = "text/css",
                                        href = "plugins/email_page.css"),
                              tags$script(src = "plugins/my_email_info.js")
                              )
                            ),
                   
                   # ----------------------------------
                   # tab panel 7 - Downlod
                   tabPanel("Downlod", icon=icon("download"),
                            downloadMenu()
                   ),
                   
                   # ----------------------------------
                   # tab panel 8 - About
                   tabPanel("About", icon=icon("question-circle"),
                            htmlTemplate("html/about.html"),
                            shinyjs::useShinyjs(),
                            tags$head(
                                tags$link(rel = "stylesheet", 
                                          type = "text/css", 
                                          href = "plugins/carousel.css"),
                                tags$script(src = "plugins/holder.js")
                            ),
                            tags$style(type="text/css",
                                       ".shiny-output-error { visibility: hidden; }",
                                       ".shiny-output-error:before { visibility: hidden; }"
                            )
                   )
                   
             ),
             # https://stackoverflow.com/questions/35584644/r-shiny-navbarpage-right-aligned-tabs
             HTML("<script>var parent = document.getElementsByClassName('navbar-nav');
                   parent[0].insertAdjacentHTML( 'afterend', '<ul class=\"nav navbar-nav navbar-right\"><li class=\"twitter\"><a href=\"https://twitter.com/hopetogy\" target=\"_blank\"><i class=\"fa fa-twitter\"></i></a></li><li class=\"weixin\"><a href=\"images/xzp_Wechat.jpg\" target=\"_blank\"><i class=\"fa fa-weixin\"></i></a></li></ul>' );
                  </script>")
))