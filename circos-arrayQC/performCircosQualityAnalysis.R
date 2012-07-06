
## preload all necessary libraries:	
	require.libs()

### prepare data:

	#### use GEO series:
	series <- "GSE18965"
	
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
	a <- sapply(strsplit(sampleNames(data),'_'),unlist)
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

	
