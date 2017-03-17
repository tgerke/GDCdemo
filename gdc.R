#######################################################################################
# https://github.com/Bioconductor/GenomicDataCommons
# to install:
# library(devtools)
# install_github("Bioconductor/GenomicDataCommons")
# a note to Mac users: don't try this with XCode 8 

library(GenomicDataCommons)
options(stringsAsFactors=FALSE)

# convenient to set a working directory where you'd like to download
setwd("~/Documents/tcga/gdc/")

#######################################################################################
# check status (should get "OK")
GenomicDataCommons::status()

#######################################################################################
# one way to filter: projects

# starts with assigning a blank query that we will later filter
projquery <- projects()
# note could also have created a query from one of the 3 other subclasses:
#caseq <- cases()
#fileq <- files()
#annq <- annotations()

projquery
# fields are the fields returned when we retrieve data
# filters gets populated with the filters() function
# facets is for aggregating data
# note that legacy can be swapped to "legacy" (maybe TRUE) to access the legacy archive

# this is going to retreive some metadata 
# the count of records satisfying filter criteria, in this case everything 
count(projquery)

# retrieve some metadata
projresults <- results(projquery)
projresults
# this looks wrong! Reason: default of results() is to return 10 records

# can get everything with results_all, but use with caution. 
# advised: test with counts() and results() before calling results_all()
projresults <- results_all(projquery)

# convenient to convert this to a dataframe
projresults <- as.data.frame(projresults)
projresults

#######################################################################################
# have a look at fields; these vary based on the query type

# the set of default fields is a reasonable subset of those available
default_fields("projects")
default_fields("files")
length(available_fields("files"))
default_fields("cases")
length(available_fields("cases"))

#######################################################################################
# one way to filter: cases

caseq <- cases()
caseq

# use select() to grab the fields you want
caseq <- GenomicDataCommons::select(cases(), available_fields('cases'))
length(caseq$fields)

# this is a useful function
grep_fields("cases", "infiltration")
grep_fields("annotations", "id")
# so is this (loads a shiny app)
field_picker("cases")

#######################################################################################
# within fields, we likely need to know the subfields available

default_fields("files")
types <- aggregations(facet(files(), c("access", "type")))
types

#######################################################################################
# query some metadata to find the files we want to download

# how many total files
fileq <- files()
count(fileq)

# only gene expression files
aggregations(facet(files(), "type"))
fileq <- filter(files(), ~ type=="gene_expression")

#######################################################################################
# now that we have all of the building blocks, let's build a query
# want project TCGA PRAD gene expression from tumor tissue
# here's a way to do it using a file query

# where might we subset files by PRAD project?
grep_fields("files", "proj")
# try a couple of these
aggregations(facet(files(), c("cases.project.primary_site", "cases.project.project_id")))
# looks like either would work, but let's go with the TCGA-PRAD one

fileq <- filter(files(), ~ cases.project.project_id=="TCGA-PRAD" & 
                   type=="gene_expression")
# here's our manifest
manifest(fileq)
# oops, got duplicate quantities here; supposed I only want the FPKM
# (this actually takes some hunting to find the right key, which is "HTSeq - FPKM")
aggregations(facet(files(), "analysis.workflow_type"))

# try again
fileq <- filter(files(), ~ cases.project.project_id=="TCGA-PRAD" & 
                  type=="gene_expression" &
                  analysis.workflow_type=="HTSeq - FPKM")
# here's our manifest
manifest(fileq)
# still not right -- I happen to know there are <500 tumor tissue profiles + some normals

grep_fields("files", "type")
aggregations(facet(files(), "cases.samples.sample_type"))

# try again
fileq <- filter(files(), ~ cases.project.project_id=="TCGA-PRAD" & 
                   type=="gene_expression" &
                   analysis.workflow_type=="HTSeq - FPKM" &
                   cases.samples.sample_type=="Primary Tumor")
# here's our manifest
manifest(fileq)
count(fileq)

#######################################################################################
# finally, let's do the download

manifestdata <- manifest(fileq)
head(manifestdata)
# only download the first 2 files
fileids <- gdcdata(manifestdata$id[1:2], destination_dir=getwd(), progress=FALSE)

# did it work?
list.files(getwd())
library(R.utils)
zipfiles <- list.files(getwd())
for (i in 1:length(zipfiles)) {gunzip(zipfiles[i])}
readfiles <- list.files(getwd())
gep <- vector("list", length=length(readfiles))
for (i in 1:length(readfiles)) {
   gep[[i]] <- read.table(readfiles[i], sep="\t", header=FALSE)
}
head(gep[[1]])

#######################################################################################
# find the linker for the files we just downloaded

expands <- c("diagnoses","annotations", "demographic","exposures")
clinResults = cases() %>% 
  GenomicDataCommons::select(NULL) %>%
  GenomicDataCommons::expand(expands) %>% 
  results(size=50)
clinDF = as.data.frame(clinResults)

# the first case ID (not found because we only downloaded the first 50 cases!)
searchid <- substr(readfiles[1], 1, 36)
grep(searchid, clinDF$case_id)
# should work if you filter the above clinResults query

#######################################################################################
# working with the token

# gdc_token() searches for token resolved in this order: 
# As a string stored in the environment variable, GDC_TOKEN
# As a file, stored in the file named by the environment variable, GDC_TOKEN_FILE
# In a file in the user home directory, called .gdc_token

GDC_TOKEN <- read.table("gdcToken.txt")$V1
gdc_token()
gdc_token
Sys.getenv("GDC_TOKEN")
Sys.setenv(GDC_TOKEN=GDC_TOKEN)

