
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



writeCircosConfig <- function(path, name, sepp){

	fName<- "circosQCconfig.txt"
		sink(fName, append=TRUE)
			cat(paste("<image>", "<<include etc/image.conf>>",
			paste("file*=", name, sep=''),
			"angle_offset* = -88",
			"</image>", sep='\n'))
			cat('\n')
			cat('\n')
			cat(paste("show_pca=yes","show_degradation=yes","show_outliers_general=yes","show_outliers_bio=yes","show_outliers_internal.a=yes","show_outliers_internal.b=yes",paste("karyotype=", paste(path, "Karyotype.txt", sep=sepp)),"chromosomes_order_by_karyotype=yes","chromosomes_units=1","chromosomes_display_default=yes",sep='\n'))
			cat('\n')
			cat('\n')
			cat(paste("<plots>", "<plot>", "show=conf(show_pca)", "type=scatter", paste("file=", paste(path, "PCA.txt", sep=sepp)),
			"r1=0.99r","r0=0.8r","min=-35","max=20", "fill_color=blue","stroke_color=vdgrey","stroke_thickness=1","glyph=circle", "glyph_size=18","<<include background.conf>>",
			# "<axes>","color     = grey_a1","thickness = 1","<axis>","position = 0","color    = black_a1","</axis>","<axis>","spacing = 5","</axis>","</axes>",
			"axis = yes","axis_color     = black_a1","axis_thickness = 1","axis_spacing   = 15",
			"<rules>","<rule>", "condition = var(id) eq 'A'","fill_color=red", "stroke_color=black","glyph_size=24","</rule>","</rules>","</plot>",sep='\n')) 
			cat('\n')
			cat('\n')
			cat(paste("<plot>","show = conf(show_degradation)", "type = heatmap", paste("file=", paste(path, "RNAdegradation.txt", sep=sepp)),"r1 = 0.77r","r0 = 0.72r","stroke_thickness = 1", "stroke_color = grey_a3", "color = rdbu-11-div","<<include background.conf>>","</plot>", sep='\n'))
			cat('\n')
			cat('\n')
			cat(paste("<plot>","show = conf(show_outliers_general)", "type  = scatter", paste("file=", paste(path, "Outliers_General.txt", sep=sepp)),"r1 = 0.69r","r0 = 0.60r","fill_color = red","<<include outliers.glyphs.conf>>","<<include background.conf>>","</plot>", sep='\n'))
			cat('\n')
			cat('\n')
			cat(paste("<plot>", "show = conf(show_outliers_bio)", "type  = scatter", paste("file=", paste(path, "Outliers_Bio.txt", sep=sepp)),"r1 = 0.59r", "r0 = 0.50r", "fill_color = blue", "<<include outliers.glyphs.conf>>", "<<include background.conf>>","<rules>", "<rule>",  "condition  = var(id) eq 'A'", "fill_color = red", "stroke_color = black", "glyph_size = 22", "</rule>", "</rules>", "</plot>", sep='\n'))
			cat('\n')
			cat('\n')
			cat(paste("<plot>", "show = conf(show_outliers_internal.a)", "type  = scatter", paste("file=", paste(path, "Outliers_Internal.A.txt", sep=sepp)),"r1 = 0.49r", "r0 = 0.40r", "fill_color = red", "<<include outliers.glyphs.conf>>", "<<include background.conf>>","</plot>"	, sep='\n'))
			cat('\n')
			cat('\n')
			cat(paste("<plot>", "show = conf(show_outliers_internal.b)", "type  = scatter", paste("file=", paste(path, "Outliers_Internal.B.txt", sep=sepp)),"r1 = 0.39r", "r0 = 0.30r", "fill_color = red", "<<include outliers.glyphs.conf>>", "<<include background.conf>>", "</plot>", "</plots>", sep='\n'))
			cat('\n')
			cat('\n')
			cat(paste("<<include ideogram.conf>>", "<<include ticks.conf>>", "<<include etc/colors_fonts_patterns.conf>>", "<<include etc/housekeeping.conf>>", sep='\n'))
			cat('\n')
		sink()
	}
	
	
    write.ticks <- function(){
		fName<- "ticks.conf"
		sink(fName, append=TRUE)
		  cat(paste("show_ticks       = yes","show_tick_labels = yes","show_grid        = yes","<ticks>","radius           = 1r","color            = black","thickness        = 2p","format           = %d","grid_end       = 1r","grid_start     = 0.3r","<tick>","spacing        = 5u","size           = 5p","</tick>","<tick>","spacing        = 20u","size           = 12p","show_label     = yes","label_size     = 20p","label_offset   = 10p","format         = %d","grid = yes","grid_color = grey_a3","grid_thickness = 1p","</tick>","</ticks>", sep='\n'))
		sink()
	}
	
    write.ideogram <- function(){
	fName<- "ideogram.conf"
		sink(fName, append=TRUE)
		  cat(paste("<ideogram>","<spacing>","default = 0.002r","<pairwise hs12 hs2>","spacing = 10r","</pairwise>","</spacing>","thickness = 1p","fill = yes","fill_color = black","radius = 0.90r","show_label = yes","label_parallel = yes","label_size = 36p","label_font = default","label_radius = dims(image,radius)-50p","show_bands = no","</ideogram>", sep='\n'))
		sink()
	}
	 write.background <- function(){
		fName<- "background.conf"
		sink(fName, append=TRUE)
		  cat(paste("background       = yes","background_color            = grey_a5","background_stroke_color     = grey_a3","background_stroke_thickness = 1", sep='\n'))
		sink()
	}
	write.outliers.glyphs <- function(){
		fName<- "outliers.glyphs.conf"
		sink(fName, append=TRUE)
		  cat(paste("stroke_color = black","stroke_thickness = 1","glyph = circle","glyph_size = 20", sep='\n'))
		sink()
	}

writeHighlights <- function(cel, rown){

	fName<- "Highlights.txt"

 	for(i in 1:length(cel)){
		sink(fName, append=TRUE)
			cat(paste(paste("hs", i, sep=''), paste(1, rown+1, sep='\t'), sep='\t'), sep='\n')
		sink()
	}
}



# chr - hs9 GSM251362 1 21  chr9
# chr - hs10 GSM251363 1 21  chr10
# band hs1 1007_s_at 1007_s_at 1 2 gneg

oddEven <- function(i){
	res <- ifelse(i %%2, "gneg", "gpos") 
	
	if(res == "gpos"){
		ifelse(i %%5, "gpos25", "gpos50")
	}
	else
		res
}


writeKaryotype <- function(cel, probes, d){

fName <- "Karyotype.txt"
	rown <- length(probes)
	
	## KaryotypeHeader:
 	for(i in 1:length(cel)){	
		sink(fName, append=TRUE)
			cat(paste('chr', '-', paste("hs", d[i,], sep=''), sub('.CEL', '', cel[d[i,]]), 1, rown+1, 'white',  sep=' '))
			cat('\n')
		sink()
	}

	
	## KaryotypeBand:
	for(i in 1:length(cel)){
		sink(fName, append=TRUE)
			for(j in 1:length(probes)){
				cat(paste('band', paste("hs", i, sep=''), probes[j], probes[j], j, j+1, oddEven(j), sep=' '))
				cat('\n')			
			}
		sink()
	}
	
}

writePCA <- function(eset, size){

	fName <- "PCA.txt"
	pcaResult <- pca(eset, method="svd", center=TRUE, nPcs=2)

	allsamples <- dim(scores(pcaResult))[1]
	size <- (size - allsamples)/allsamples
		
	  
	sink(fName, append=TRUE)
		for(i in 1:dim(scores(pcaResult))[1]){
			for(j in 1:dim(scores(pcaResult))[1]){
				score <- scores(pcaResult)[j,1]
				cat(paste('hs',i,sep=''), round(j*size, digits = 0),round((j+1)*size, digits = 0), score, getlabel(i,j) , sep='\t')
				cat('\n')	
			}
		}
	sink()
}


getlabel <- function(i,j){
	ifelse(i == j, 'id=A','id=B' )
}

writeRNAdegradation <- function(data, size){

  fName <- "RNAdegradation.txt"
    size <- (size - 11)/11
	result <- list()
	
	if(file.exists(fName)){
		print(paste(paste('The File: ', fName, sep=" "), 'is already there.', sep=" "))
		print('Maybe you are trying to overwrite existsting results?')
		print('Aborted..')
	}
	else {
		  for(i in 1:length(sampleNames(data))){
			rnaDeg <-  AffyRNAdeg(data[,i])
			
			sink(fName, append=TRUE)
			#res <- shift.scale(rnaDeg)
			res <- rnaDeg$means.by.number
				for(j in 1: length(res)){
					cat(paste('hs',i,sep=''),round(j*size, digits = 0),round((j+1)*size, digits = 0), res[j], sep='\t')
					cat('\n')
				}
			sink()
			result[i] <- sum(res)
		  }
	 }
  return(result)
}

shift.scale <- function(rnaDeg){
		#(transform == "shift.scale")
		sds <- rnaDeg$ses
		mns <- rnaDeg$means.by.number
		 
 mn <- mns[, 1]
 mns <- sweep(mns, 1, mn)
 mns <- mns/(sds)
 sweep(mns, 1, 1:(dim(mns)[1]), "+")
}

#hs7 0 4999999 9522.000000
#hs7 5000000 9999999 30013.000000
#hs7 10000000 14999999 14083.000000
#hs7 15000000 19999999 10505.000000
#hs7 20000000 24999999 36482.000000

writeHeatmap <- function(probes, fName){
	mprobes <- mean(probes)
	
  	for(i in 1:length(colnames(probes))){
		sink(fName, append=TRUE)
			for(j in 1:length(rownames(probes))){
				cat(paste(paste("hs", i, sep=''), j,(j+1), getneg(probes[j,i], mprobes), sep=' '))
				cat('\n')			
			}
		sink()
	}
}


getneg <- function(value, mdata){
	ifelse(value > mdata, value, -value)
}
	


 writeOutliers.General<- function(yqc, snames, size){
 
	Out.SFS <- getOutliers(yqc, "sfs") 		# scale factor
	Out.AVBG <- getOutliers(yqc, "avbg")	# average background
	Out.AVNS <- getOutliers(yqc, "avns") 	# average noise
	#Out.PP <- getOutliers(yqc, "pp")		# percentage present
	 
	fName <- "Outliers_General.txt"
	 
	allsamples <- length(snames)
	size <- (size - allsamples)/allsamples
	
	sink(fName, append=TRUE)
		for(i in 1:length(snames)){
			if(length(Out.SFS) != 0){
				for(j in 1:length(Out.SFS)){
					if(snames[i] == names(Out.SFS[j])){
						cat(paste(paste("hs", i, sep=''), round(i*size, digits = 0),round((i+1)*size, digits = 0), Out.SFS[j], sep=' '))
						cat('\n')			
					}
				}
			}
			if(length(Out.AVBG) != 0){
				for(j in 1:length(Out.AVBG)){
					if(snames[i] == names(Out.AVBG[j])){
						cat(paste(paste("hs", i, sep=''), round(i*size, digits = 0),round((i+1)*size, digits = 0), Out.AVBG[j], sep=' '))
						cat('\n')			
					}
				}
			}
			if(length(Out.AVNS) != 0){
				for(j in 1:length(Out.AVNS)){
					if(snames[i] == names(Out.AVNS[j])){
						cat(paste(paste("hs", i, sep=''), round(i*size, digits = 0),round((i+1)*size, digits = 0), Out.AVNS[j], sep=' '))
						cat('\n')			
					}
				}
			}
		}
	sink()
}

 writeOutliers.Bio<- function(yqc, snames, size){
 
	yqcRatios <- getQCRatios(yqc)
	Out.Actin <- getOutliers(yqc, "actin")	# beta-actin 3/5 ratio
	Out.GADPH <- getOutliers(yqc, "gapdh")	# GAPDH 3/5 ratio
	 
	fName <- "Outliers_Bio.txt"
	 
	allsamples <- length(snames)
	size <- (size - allsamples)/allsamples
	 
	sink(fName, append=TRUE)
	# write actin:
	for(i in 1:length(snames)){
		cat(paste(paste("hs", i, sep=''), round(i*size, digits = 0),round((i+1)*size, digits = 0), yqcRatios[1,i], 'id=B', sep=' '))
		cat('\n')			
	}	
	# write gapdh:
	for(i in 1:length(snames)){
		cat(paste(paste("hs", i, sep=''), round(i*size, digits = 0),round((i+1)*size, digits = 0), yqcRatios[2,i], 'id=B', sep=' '))
		cat('\n')			
	}
	# write outliers:
	if(length(Out.Actin) != 0){
		for(j in 1:length(Out.Actin)){
			for(i in 1:length(snames)){
				if(snames[i] == names(Out.Actin[j])){
					cat(paste(paste("hs", i, sep=''), round(i*size, digits = 0),round((i+1)*size, digits = 0), Out.Actin[j], 'id=A', sep=' '))
					cat('\n')			
				}
			}
		}		
	}
	if(length(Out.GADPH) != 0){
		for(j in 1:length(Out.GADPH)){
			for(i in 1:length(snames)){
				if(snames[i] == names(Out.GADPH[j])){
					cat(paste(paste("hs", i, sep=''), round(i*size, digits = 0),round((i+1)*size, digits = 0), Out.GADPH[j], 'id=A', sep=' '))
					cat('\n')			
				}
			}
		}
	}
	sink()
}

 writeOutliers.Internal.A <- function(yqc, snames, size){
 
	 Out.bioB <- getOutliers(yqc, "biob")	# internal bioB control
	 Out.bioC <- getOutliers(yqc, "bioc")	# internal bioC control
	 Out.bioD <- getOutliers(yqc, "biod")	# internal bioD control
	 
	fName <- "Outliers_Internal.A.txt"
	 
	allsamples <- length(snames)
	size <- (size - allsamples)/allsamples
	 
	sink(fName, append=TRUE)
		for(i in 1:length(snames)){
			if(length(Out.bioB) != 0){
				for(j in 1:length(Out.bioB)){
					if(snames[i] == names(Out.bioB[j])){
						cat(paste(paste("hs", i, sep=''), round(i*size, digits = 0),round((i+1)*size, digits = 0), Out.bioB[j], sep=' '))
						cat('\n')			
					}
				}
			}
			if(length(Out.bioC) != 0){
				for(j in 1:length(Out.bioC)){
					if(snames[i] == names(Out.bioC[j])){
						cat(paste(paste("hs", i, sep=''), round(i*size, digits = 0),round((i+1)*size, digits = 0), Out.bioC[j], sep=' '))
						cat('\n')			
					}
				}
			}
			if(length(Out.bioD) != 0){
				for(j in 1:length(Out.bioD)){
					if(snames[i] == names(Out.bioD[j])){
						cat(paste(paste("hs", i, sep=''), round(i*size, digits = 0),round((i+1)*size, digits = 0), Out.bioD[j], sep=' '))
						cat('\n')			
					}
				}
			}
		}
	sink()
}

 writeOutliers.Internal.B <- function(yqc, snames, size){
 
	 Out.dap <- getOutliers(yqc, "dap")		# Dap spike control
	 Out.thr <- getOutliers(yqc, "thr")		# thr spike control
	 Out.phe <- getOutliers(yqc, "phe")		# phe spike control
	 Out.lys <- getOutliers(yqc, "lys")		# lys spike control
 
	fName <- "Outliers_Internal.B.txt"
	 
	allsamples <- length(snames)
	size <- (size - allsamples)/allsamples
	 
	sink(fName, append=TRUE)
		for(i in 1:length(snames)){
			if(length(Out.dap) != 0){
				for(j in 1:length(Out.dap)){
					if(snames[i] == names(Out.dap[j])){
						cat(paste(paste("hs", i, sep=''), round(i*size, digits = 0),round((i+1)*size, digits = 0), Out.dap[j], sep=' '))
						cat('\n')			
					}
				}
			}
			if(length(Out.thr) != 0){
				for(j in 1:length(Out.thr)){
					if(snames[i] == names(Out.thr[j])){
						cat(paste(paste("hs", i, sep=''), round(i*size, digits = 0),round((i+1)*size, digits = 0), Out.thr[j], sep=' '))
						cat('\n')			
					}
				}
			}
			if(length(Out.phe) != 0){
				for(j in 1:length(Out.phe)){
					if(snames[i] == names(Out.phe[j])){
						cat(paste(paste("hs", i, sep=''), round(i*size, digits = 0),round((i+1)*size, digits = 0), Out.phe[j], sep=' '))
						cat('\n')			
					}
				}
			}
			if(length(Out.lys) != 0){
				for(j in 1:length(Out.lys)){
					if(snames[i] == names(Out.lys[j])){
						cat(paste(paste("hs", i, sep=''), round(i*size, digits = 0),round((i+1)*size, digits = 0), Out.lys[j], sep=' '))
						cat('\n')			
					}
				}
			}
		}
	sink()
}

get.path <- function(workdir, pathToCircos){

	a<-unlist(strsplit(pathToCircos,'/'))
	b<-unlist(strsplit(workdir,'/'))
	
	path <- NULL
	i <- NULL
	for(i in length(a)){
	  if(b[i]!=a[i])
	    break
	}
	for(i in length(b):i){
	    path <-paste(b[i], path,  sep="/")

	}
	path <- paste("..", path, sep="/")
	path <- paste(path, "circosFiles", sep="")
	path <- gsub("/","\\",path,fixed=TRUE) 
	
  return(path)
}


require.libs <- function(){

	require(Biobase)
	require(GEOquery)
	require(affy)
	require(pcaMethods)
	require(genefilter)
	require(IRanges)
	require(yaqcaffy)
	
	library(Biobase)
	library(GEOquery)
	library(affy)
	library(pcaMethods)
	library(genefilter)
	library(IRanges)
	library(yaqcaffy)
	# source("http://www.bioconductor.org/biocLite.R")
	
}


writeCircos.files <- function(data, cel, workdir, fileName, pathToCircos ){

      ## get libraries:
	require.libs()

	## normalize the data:
	eset.norm <- rma(data, destructive=TRUE, normalize=TRUE, background=FALSE)

	dir.create("circosFiles")
	setwd("circosFiles") 


	n_probes <- nrow(eset.norm)
	noa <- c(1:(length(sampleNames(eset.norm))))	# number if arrays
	nop <- c(1:200)	# number of probes

	
	### RNA degradation:
	degrad <- writeRNAdegradation(data[,noa],  length(nop))
	d <- data.frame(noa, unlist(degrad), sort(unlist(degrad)))
	d <- d[order(d[,2]),]
	d <-  d[1]

	### Karyotype & Highlites:
	writeKaryotype(sampleNames(data)[noa], featureNames(eset.norm)[nop], d)
	writeHighlights(sampleNames(data)[noa],  length(nop))
	 
	### PCA:
	writePCA(eset.norm[,noa],  length(nop))
	 
	###   yaqcaffy
	yqc <-yaqc(data[,noa])

	snames <- sampleNames(data)[noa]

	 ## get OUtliers
	 writeOutliers.General(yqc, snames, length(nop))
	 writeOutliers.Bio(yqc, snames, length(nop))
	 writeOutliers.Internal.A(yqc, snames, length(nop))
	 writeOutliers.Internal.B(yqc, snames, length(nop))
	 
	 
	### Circos configuration files:
	if(.Platform$OS.typ == "windows"){
	  path <- get.path(workdir, pathToCircos)
	  writeCircosConfig(path, fileName, sepp="\\")
	}
	else{
	path <- paste(workdir, "circosFiles", sep="/")
	  writeCircosConfig(path, fileName, sepp="/")
	}
	
	write.ticks()
	write.ideogram()
	write.background()
	write.outliers.glyphs()
	 
	print("Circos files are generated! please go to circos home dir and execute:")
	print("")
	print(paste("~/circos-0.60>   perl  bin/circos -conf",paste(path, "circosQCconfig.txt", sep='/'), sep=' '))
	 
	 
}




