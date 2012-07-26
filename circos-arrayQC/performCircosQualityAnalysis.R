
#/**
#* Copyright (C) [2012] [martin koch] This program is free software; you can
#* redistribute it and/or modify it under the terms of the GNU General
#* Public License as published by the Free Software Foundation; either
#* version 3 of the License, or (at your option) any later version. This
#* program is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
#* more details. You should have received a copy of the GNU General Public
#* License along with this program; if not, see
#* <http://www.gnu.org/licenses/>.
#**/

source("C:/Users/mak/git/circos-arrayQC/circos-arrayQC/makeCircosFiles.R")
setwd("Arrays/projectGEO_Quality/")


## preload all necessary libraries:	
	require.libs()

### prepare data:

	#### use GEO series:
	series <- "GSE18965"
	series <- "GSE37200"
	series <- "GSE9936"

	# getfiles:
	getGEOSuppFiles(as.character(series))

	# create dirs:
	dataDir <- paste(series,"data", sep="_")
	dir.create(dataDir)
	untar(paste(paste(series,series,sep="/"),"RAW.tar", sep="_"), exdir=dataDir)
	setwd(dataDir)
		
	cel_files <- dir(pattern="*.gz") 
	sapply(cel_files, gunzip)
	
  
	# get data from cel files
	celNames <- dir(pattern="*.CEL") 
	data <- read.affybatch(filenames=celNames)

	## get rid of the ".CEL" filesuffix:
	sampleNames(data) <- sub('.CEL','',sampleNames(data))
	a <- sapply(strsplit(sampleNames(data),'_'),unlist) # a<-a[1,]

	sampleNames(data)  <- sub('GSM', '', unlist(a))
	print(sampleNames(data))	
	
## write circos files::
	workdir <- getwd()
	pathToCircos <- "C:/Users/mak/Documents/circos-0.6"
	#pathToCircos <- "/home/mak/mak/circos-0.60"

	fileName <- paste("circosQC_", series, sep='')

	writeCircos.files(data, celNames, workdir, fileName, pathToCircos)

	# clean up 
	rm(data, cel)

	
	
	expressionQCPipeline
	
	##########################################################
	# use Circos for illumina arrays
	library("GEOquery")
	library("beadarray")	
	library("beadarrayExampleData")	
	data(exampleBLData)
	library("lumi")
	library("lumiHumanIDMapping")
	library("lumiHumanAll.db")
	  source("http://bioconductor.org/biocLite.R")
    biocLite("illuminaHumanv3.db")

	library("illuminaHumanv3.db")
	
	
	series <- "GSE10220"
	series <- "GSE13290"
	series <- "GSE8686"

	
	setwd("Arrays/projectGEO_Quality/")

	
	# getfiles:
	getGEOSuppFiles(as.character(series))

	# create dirs:
	dataDir <- paste(series,"data", sep="_")
	dir.create(dataDir)
	untar(paste(paste(series,series,sep="/"),"RAW.tar", sep="_"), exdir=dataDir)
	setwd(dataDir)
		
	files <- dir(pattern="*.gz") 
	sapply(files, gunzip)
	
	# get data from gpr files
	gprNames <- dir(pattern="*.") 

	eset <- lumiExpresso(lumiR.batch(gprNames, na.rm = TRUE, lib.mapping="lumiHumanIDMapping",  checkDupId = TRUE, QC="FALSE"))

	BLData <- readIllumina(useImages = FALSE, illuminaAnnotation = "lumiHumanIDMapping")
	suggestAnnotation(BLData, verbose = TRUE)
	annotation(exampleBLData) <- "Humanv3"
	BLData <- 
	
	p <- combinedControlPlot(exampleBLData)

	annoName <- annotation(exampleBLData)
    cProf <- makeControlProfile(annoName)
	array = gprNames[1]

	inten <- getBeadData(exampleBLData, array = array, what = "Grn")
    pIDs <- getBeadData(exampleBLData, array = array, what = "ProbeID")
	wts <- getBeadData(data, array, what = wtsName)
	
	 bsv <- identifyControlBeads(exampleBLData, array = 1, "Humanv3")
	
	## get rid of the ".CEL" filesuffix:
	sampleNames(data) <- sub('.CEL','',sampleNames(data))
	a <- sapply(strsplit(sampleNames(data),'_'),unlist) # a<-a[1,]

	sampleNames(data)  <- sub('GSM', '', unlist(a))
	print(sampleNames(data))	
	
	
	
	
	
	
	
	##########################################################
	# use Circos for illumina arrays
	#source("http://bioconductor.org/biocLite.R")
    #biocLite("affyPLM")
	library("affyPLM")	
	library("arrayQualityMetrics")	
	library("lumi")
	
	setwd("Arrays/projectGEO_Quality/")

	
	series <- "GSE6968"
	dataDir <- paste(series,"data", sep="_")
	setwd(dataDir)
		
	gprNames <- dir(pattern="*.") 

	eset <- lumiExpresso(lumiR.batch(gprNames, na.rm = TRUE, lib.mapping="lumiHumanIDMapping",  checkDupId = TRUE, QC="FALSE"))

	qa <- arrayQualityMetrics(expressionset = eset, outdir = series, force = TRUE, do.logtransform = TRUE)
	
	outliers(exprs(eset)[1:10,1:10], method = c("KS", "sum", "upperquartile"))
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
