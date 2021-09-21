library(tidyr)
library(dplyr)
library(emayili)
library(magrittr)
library(readr)

#   ____________________________________________________________________________
#   Database modified Day                                                   ####

modifiedDay <- "2021-05-23"

#   ____________________________________________________________________________
#   Send Email                                                              ####

adminEmail <- "1099808298@qq.com"
adminEmailPassword <- "otaudkznvbgthihf"

smtp <- server(host = "smtp.qq.com",
               port = 587,
               username = adminEmail,
               password = adminEmailPassword)

#   ____________________________________________________________________________
#   Reviewed Marker Gene Database                                           ####

rmgddf <- read.csv("www/data/ReviewedMarkerGeneDatabase.csv", stringsAsFactors = FALSE, check.names = FALSE)
rmgd <- unite(rmgddf, xlab, c(Tissue, CellType), remove=FALSE)

print("Reviewed Marker Gene Database Complete")

#   ____________________________________________________________________________
#   Enzymolysis Solution                                                    ####

enzymolysisSolution <- read_delim("www/data/enzymolysisSolution.csv",delim = "\t", col_types = cols(.default = "c"))

print("Protocol Enzymolysis Solution Data Complete")

#   ____________________________________________________________________________
#   Protocol Troubleshooting                                                ####

troubleshooting <- read_delim("www/data/ProtocolTroubleshooting.csv",delim = "\t", col_types = cols(.default = "c"))

print("Protocol Troubleshooting Data Complete")
