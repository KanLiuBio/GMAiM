library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)
library(dplyr)
library(edgeR)
library(ggplot2)




### All the inputs and options to run the package. ### 

# Input miRNA.dat file from miRBase
datFile <- "/work5/liukan/Imprecision/test/miRNA.dat"
# Give the absolute directory path of all the sorted and indexed bam files.
bamDir <- "/work5/liukan/Imprecision/test"


# Output csv file for the detailed miRNA imprecision read count table
outputFile1 <- "/work5/liukan/Imprecision/test/Detailed_imprecision.csv"
# Output csv file for the final differential imprecision comparision table.
outputFile2 <- "/work5/liukan/Imprecision/test/DE_Comparison.csv"
# Output csv file for the imprecision plot.
outputFile3 <- "/work5/liukan/Imprecision/test/deviation.csv"
# Output csv file for the imprecision plot 2.
outputFile4 <- "/work5/liukan/Imprecision/test/deviation_detail.csv"
# Output directory to store all the miRNA imprecision plot.
dirOut= '/work5/liukan/Imprecision/test/Plot/'

# Species three-letter abbreviations in miRBase.
# Please find this information from http://www.mirbase.org/cgi-bin/browse.pl
myid='ath'
#Specify the detailed sample names for comparison
SampleList=c('col','smx')

# How many total mismatches tolerate for imprecision cutoff
total_bias=2

# How many replicates for each sample.
repNum=2
if (repNum <2) {
  print ("Error: cannot parse whether replicates or not.")
  exit()
}
# Which type of imprecision the user wants to compare with the precision miRNA.
# Please choose one of the four types (case sensitive): mismatch, longer,shorter, shifted.
my_flag="mismatch"

# Flanking nucleotide size of upstream and downstream mature miRNA for plot
boundary_size=5











#### DE comparison function ####
mTest.edgeR <- function(counts, groups,samples=NULL, genes=rownames(genes),names.counts = c("precision","imprecision"),thr.count=5,thr.sample=round(ncol(counts)/2),verbose=T,counts.attached=T){
  n = length(groups)  
  ugp = unique(groups)
  ng = length(ugp)
  
  if(is.null(samples)){
    samples=groups
    for(gg in ugp){
      ixgg = (groups==gg)
      ngg = sum(ixgg)
      samples[ixgg] = paste(samples[ixgg],(1:ngg),sep='_')
    }
  }
  ###
  y <- DGEList(counts=counts, group=rep(groups,each=2), genes=genes)
  
  ### filter genes with low total counts in most samples
  counts.total <- t(rowsum(t(y$counts), group=gl(n,2)))
  keep <- rowSums(counts.total >= thr.count) >=thr.sample
  if(verbose){
    print('The numbers of genes filtered and kept are')
    print(table(keep))
  }
  ngene=sum(keep)
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  ### contrcut the design matrix
  sam <- factor(rep(samples, each=2),levels=samples)
  stts <- factor(rep(names.counts,n), levels=names.counts)
  design <- model.matrix(~ sam + stts)
  colnames(design) <- gsub("sam","",colnames(design))
  colnames(design) <- gsub("stts","",colnames(design))
  colnames(design)[1] <- "Int"
  noff = ncol(design)
  design.diff = matrix(0,2*n,ng-1)
  colnames(design.diff) = paste0(ugp[-1],'-',ugp[1])
  for(ii in 1:(ng-1)){
    idgg = which(groups==ugp[ii+1])
    design.diff[2*idgg-1,ii] = 1
  }
  design <- cbind(design, design.diff)
  if(verbose){
    print('The design matrix is')
    print(design)
  }
  ### estimate the over dispersion
  y <- estimateDisp(y, design=design, trend="none")
  #  y$common.dispersion
  #  summary(y$prior.df)
  
  ### start the analysis 
  fit <- glmFit(y, design)
  
  output = list()
  for(ii in 2:ng){
    lrt.ii <- glmLRT(fit, coef=noff+ii-1)
    tbl.lrt.ii = topTags(lrt.ii,n=ngene,sort.by = 'none')
    output.ii = tbl.lrt.ii@.Data[[1]]
    if(counts.attached) output.ii = cbind(output.ii,y$counts)
    output.ii = output.ii[order(output.ii$PValue),]
    name.ii = colnames(design)[noff+ii-1]
    output[[name.ii]] = output.ii
  }
  if(ng>2){
    for(ii in 2:(ng-1)){
      for(jj in (ii+1):ng){
        contrast.iijj = rep(0,ncol(design))
        contrast.iijj[noff+jj-1]=1
        contrast.iijj[noff+ii-1]=-1
        lrt.iijj <- glmLRT(fit, contrast =contrast.iijj)
        tbl.lrt.iijj = topTags(lrt.iijj,n=ngene,sort.by = 'none')
        output.iijj = tbl.lrt.iijj@.Data[[1]]
        if(counts.attached) output.iijj = cbind(output.iijj,y$counts)
        output.iijj = output.iijj[order(output.iijj$PValue),]
        name.iijj = paste0(ugp[jj],'-',ugp[ii])
        output[[name.iijj]] = output.iijj
      }
    }
  }
  if(verbose){
    print('All pair-wise contrasts considered include:')
    print(names(output))
  }
  return(output)  ### a list of 
}






# Read and extract information from "miRNA.dat" file. 
con=file(datFile, open="r")

lines=readLines(con) 
lnum=length(lines)
allm=data.frame()
myidin=0
for (i in 1:lnum){
  first2=substr(lines[i], start = 1, stop = 2)
  first610=substr(lines[i], start = 6, stop = 10)
  ff=strsplit(lines[i], "\\s+")
  if (identical(first2,"ID")) {
    mirnaID=sapply(ff, `[`, 2)
    mirnaID3=substr(mirnaID, 1,3)
    if (identical(mirnaID3,myid)){
        myidin=1;
    } else {
        myidin=0;
    }
  }
  if(myidin ==0) next
  if (identical(first2,"FT") & identical(first610,"miRNA")) {
    corrd1=sapply(ff, `[`, 3)
    nf=strsplit(corrd1,"\\.\\.")
    nstart=sapply(nf,`[`, 1)
    nstop=sapply(nf,`[`, 2)
    prd=sapply(strsplit(lines[i+2],"\\s+"), `[`, 2)
    pf=strsplit(prd, "\"")
    pname=sapply(pf,`[`, 2)
    allm=rbind(allm, data.frame("miRNA_ID"=mirnaID,"miRNA_ID_sub"=pname, "Start"=as.numeric(nstart), "End"=as.numeric(nstop)))
  }
}







#Initialize the data frame for imprecision plot.
boundary_size=5
deviated_Table <- data.frame(matrix(ncol = (boundary_size *4 + 3) , nrow = 0))

boundaryName <- vector(length=0)
for (x_n in boundary_size:1) {
  x_x <-paste('upstream', x_n, sep = "-")
  boundaryName=c(boundaryName,x_x)
}
boundaryName=c(boundaryName,'upstream-0')
for (x_n in 1:boundary_size) {
  x_x <-paste('upstream', x_n, sep = "+")
  boundaryName=c(boundaryName,x_x)
}
for (x_n in boundary_size:1) {
  x_x <-paste('downstream', x_n, sep = "-")
  boundaryName=c(boundaryName,x_x)
}
boundaryName=c(boundaryName,'downstream-0')
for (x_n in 1:boundary_size) {
  x_x <-paste('downstream', x_n, sep = "+")
  boundaryName=c(boundaryName,x_x)
}

colnames(deviated_Table) <- c('miRNA',boundaryName)

df_x <- data.frame(matrix(0, ncol = (boundary_size *4 +2) , nrow = 1))
colnames(df_x) <- boundaryName

miRNA_ID=unique(allm$miRNA_ID_sub)
for (x_m in miRNA_ID) {
  deviated_Table=rbind(deviated_Table,cbind('miRNA'= x_m,df_x))
}

# Build another detailed data frame
boundary_size=5
deviated_Table_2 <- data.frame(matrix(ncol = (boundary_size *2 + 1) * (boundary_size *2 + 1) , nrow = 0))
boundaryName2 <- vector(length=0)
for (x_k in -boundary_size:boundary_size) {
  for (x_l in -boundary_size:boundary_size) {
    boundaryName2=c(boundaryName2,paste(x_k,x_l,sep=','))
  }
}
colnames(deviated_Table_2) <- c('miRNA',boundaryName2)

df_y <- data.frame(matrix(0, ncol = (boundary_size *2 + 1) * (boundary_size *2 + 1) , nrow = 1))
colnames(df_y) <- boundaryName2

for (x_m in miRNA_ID) {
  deviated_Table_2=rbind(deviated_Table_2,cbind('miRNA'= x_m,df_y))
}







#Initialize the result as data frame.
countResult=data.frame()

allHairpinID=unique(allm$miRNA_ID)




# Read all the alignments of sorted Bam files under the given directory.
my_fls <- list.files(bamDir, recursive=TRUE, pattern="*bam$", full=TRUE)
my_BamName =  vector()
for (i in my_fls){

  # File name for list in output file1
  bamName <- basename(i)
  my_BamName <- c(my_BamName, bamName) 
  my_bam <- readGAlignments(file = i)
  reads.gr <- granges(my_bam)

  for (hairpin in allHairpinID){
    tempDF=data.frame()
    # Subset all the mature miRNA records of each hairpin from miRNA.dat file.
    tempDF=allm[allm$miRNA_ID==hairpin,]

    # Subset all the alignments of each hairpin based on its ID.
    candidate= reads.gr[reads.gr@seqnames == hairpin, ]
    test <- as.data.frame(candidate@ranges)

    # Check if the hairpin has no alignments or this hairpin has more than 2 mature miRNAs.
    if (nrow(test)==0 || nrow(tempDF)>2){
      for (kk in 1:nrow(tempDF)){
        countResult=rbind(countResult, cbind("SampleID"=bamName, tempDF[kk,], data.frame("match"=0, "mismatch"=0, "longer"=0, "shorter"=0, "shifted"=0,"other"=0 )))
      }
      next
    }

    # Remove duplicated read alignments of each hairpin ID
    test.uniq <- test %>% group_by(start,end) %>% summarise(freq=n())

    for (jj in 1:nrow(tempDF)) {

      #Initialize all the variable vectors
      tempTotalCount <- as.vector(integer(nrow(tempDF)))
      tempMatchCount <- as.vector(integer(nrow(tempDF)))
      tempMisatchCount <- as.vector(integer(nrow(tempDF)))
      tempLongerCount <- as.vector(integer(nrow(tempDF)))
      tempShorterCount <- as.vector(integer(nrow(tempDF)))
      tempShiftCount <- as.vector(integer(nrow(tempDF)))
      tempOther <- as.vector(integer(nrow(tempDF)))

      matureStart=tempDF[jj,]$Start
      matureStop=tempDF[jj,]$End
      matureID=tempDF[jj,]$miRNA_ID_sub

      # iterate the whole uniquely filtered read mapping data frame. 
      for (ii in 1:nrow(test.uniq)) {
        readStart=test.uniq[ii,]$start
        readEnd=test.uniq[ii,]$end

        newName <- ''
        newName2 <- ''
        # count deviation details for upstream
        if (abs(readStart-matureStart) <= boundary_size ) {
	  if (readStart-matureStart >0 ){
	    newName=as.character(paste('upstream',as.character(readStart-matureStart),sep='+')) ## such as upstream+2 
	  } else {
	    newName=as.character(paste('upstream',as.character(matureStart-readStart),sep='-')) ## such as upstream-2 or upstream-0
	  }
	  deviated_Table[deviated_Table$miRNA==matureID,newName] = as.numeric(deviated_Table[deviated_Table$miRNA==matureID,newName]) + test.uniq[ii,]$freq
	}

        if (abs(readEnd-matureStop) <= boundary_size) {
	  if (readEnd-matureStop >0 ){
	    newName=as.character(paste('downstream',as.character(readEnd-matureStop),sep='+')) ## such as downstream+2 
	  } else {
	    newName=as.character(paste('downstream',as.character(matureStop-readEnd),sep='-')) ## such as downstream-2 or downstream-0
	  }
	  deviated_Table[deviated_Table$miRNA==matureID,newName] = as.numeric(deviated_Table[deviated_Table$miRNA==matureID,newName]) + test.uniq[ii,]$freq
	}


        if (abs(readStart-matureStart) <= boundary_size & abs(readEnd-matureStop) <= boundary_size) {
	  newName2=as.character(paste(as.character(readStart-matureStart),as.character(readEnd-matureStop),sep=','))
	  deviated_Table_2[deviated_Table_2$miRNA==matureID,newName2] = as.numeric(deviated_Table_2[deviated_Table_2$miRNA==matureID,newName2]) + test.uniq[ii,]$freq
	}


        # Check if the unique read alignment intersect with mature miRNA
        if ((readStart >= matureStart & readStart <= matureStop) | (readEnd >= matureStart & readEnd <= matureStop)) {
          tempTotalCount[jj] = tempTotalCount[jj]+test.uniq[ii,]$freq
          left_bias <- abs(readStart - matureStart)
          right_bias <- abs(readEnd -matureStop)
          if (left_bias + right_bias <= total_bias){  ### Considered as exact match reads
            tempMatchCount[jj] <-tempMatchCount[jj] + test.uniq[ii,]$freq
          }
          else{
            tempMisatchCount[jj] <-tempMisatchCount[jj] + test.uniq[ii,]$freq
            if (readStart <= matureStart & readEnd >= matureStop) {
              tempLongerCount[jj] <- tempLongerCount[jj] + test.uniq[ii,]$freq
            } else if (readStart >= matureStart & readEnd <= matureStop) {
              tempShorterCount[jj] <- tempShorterCount[jj] + test.uniq[ii,]$freq
            } else if (readStart > matureStart & readEnd > matureStop) {
              tempShiftCount[jj] <- tempShiftCount[jj] + test.uniq[ii,]$freq
            } else if (readStart < matureStart & readEnd < matureStop) {
              tempShiftCount[jj] <- tempShiftCount[jj] + test.uniq[ii,]$freq
            } else{
              tempOther[jj] <- tempOther[jj] + test.uniq[ii,]$freq
            }
          }
        }
      }
      countResult=rbind(countResult, cbind("SampleID"=bamName, tempDF[jj,], data.frame(  "match"=as.numeric(tempMatchCount[jj]), "mismatch"=as.numeric(tempMisatchCount[jj]), "longer"=as.numeric(tempLongerCount[jj]), "shorter"=as.numeric(tempShorterCount[jj]), "shifted"=as.numeric(tempShiftCount[jj]) , "other"=as.numeric(tempOther[jj]) )))
    }
  }
}





#Tranform into DE analysis format
new_table=data.frame()
for (miRNA in countResult[countResult$SampleID==my_BamName[1],]$miRNA_ID_sub){
  new_table <- rbind( new_table, data.frame("miRNA"=miRNA))
}

for (bamName in my_BamName){
  subset <- countResult[countResult$SampleID==bamName,]
  subset2 <-subset[,c("match",my_flag)]
  precision=paste(bamName,"match",sep='_')
  imprecision=paste(bamName,my_flag,sep='_')
  colnames(subset2) <- c(precision,imprecision)
  new_table=cbind(new_table,subset2)
}

nmgn2 = new_table[,1]
ctbl2=new_table[-1]

DE_result = mTest.edgeR(ctbl2,groups = rep(SampleList,each=repNum),genes = nmgn2)


#Output results
write.table(countResult, file = outputFile1, sep = ",", col.names = NA)
write.table(DE_result, file = outputFile2, sep = ",", col.names = NA)
write.table(deviated_Table, file = outputFile3, sep = ",", col.names = NA)
write.table(deviated_Table_2, file = outputFile4, sep = ",", col.names = NA)



###Generate imprecision plot figures for all miRNAs.

df_detail = read.csv(outputFile4, header = T,check.names = F)
df_detail2=df_detail[-1]
lineCount=length(df_detail2)

for (i in 1:lineCount){

  indi_df <- data.frame()
  test <- df_detail2[i,]
  miRNA <- df_detail2[i,]$miRNA
  test <- test[-1]
  total_count <- rowSums(test)
  if (total_count == 0){
    next
  }
  for (j in colnames(test)){
    pair<-strsplit(j, "\\,+")
    up<-sapply(pair, `[`, 1)
    down<-sapply(pair, `[`, 2)
    freq<-as.numeric(test[,j])
    if (freq >=1) {
      indi_df=rbind(indi_df,data.frame("up"=up,"down"=down,"freq"=freq))
    }
  }

  pdf_file=paste(dirOut,miRNA,sep='')
  pdf_file=paste(pdf_file,"pdf",sep='.')
  pdf(pdf_file)

  p6 <- ggplot( indi_df, aes(x = up, y = down, size = freq)) +
        geom_point(shape = 21, colour = "mediumvioletred",fill = "springgreen") +
        ggtitle("miRNA imprecision frequency") +
        labs(x = "Upstream deviation", y = "Downstream deviation")
  plot(p6)

  dev.off()

}













