###################################
####    RNA-seq Pipeline       ####
#######     Min Heo          ######
#####  2021 1-14 ~ 2021 - 2-2  ####
###################################

setRepositories(ind=c(1:8))

library(HelpersMG)
library(bfork)
### Directory Inputs
Args <- commandArgs(trailingOnly = T)

idxRNA_seq_Directory <- which(Args=="-Drna")
idxPutYoutFastqDir <- which(Args=="-Dfq")
idxPutYourTrimProgramDir_JAVA <- which(Args=="-Dtri")  ## "/program/Trimmomatic/trimmomatic-0.39.jar"
idxPutYourFastQC_Program_Directory <- which(Args=="-Dqc")
idxPutYourPutYourHISAT2_Pro_DIR <- which(Args=="-Dhi")
idxPutYouSamtools_Pro_Dir <- which(Args=="-Dsm")
idxPutYourfeatureCounts_Pro_Dir <- which(Args=="-Dfc")
idxScientificName <- which(Args=="-sn") 
idxReferenceID <- which(Args=="-rid")
idxPutYourThreadsNumber <- which(Args=="-t")
idxPutYourTrimAdapter_Dir_ID <- which(Args=="-a")

idxEnsembl_ReleaseVersion <- which(Args=="--r")
idxTophred <- which(Args=="--p")
idxMismatches <- which(Args=="--m")
idxPalinClipThresh <- which(Args=="--pthr")
idxSimpleClipThresh <- which(Args=="--sthr")
idxLEADING <- which(Args=="--l")
idxTRAILING <- which(Args=="--t")
idxWindowsize <- which(Args=="--w")
idxReQuiredQuality <- which(Args=="--q")
idxMINLEN <- which(Args=="--ml")
idxRNA_Strandness <- which(Args=="--rna-strandness")
idxStranded_Featurecounts <- which(Args=="--s")

##file 

RNA_seq_Directory <- Args[idxRNA_seq_Directory+1]
PutYoutFastqDir <- Args[idxPutYoutFastqDir+1]
PutYourTrimProgramDir_JAVA <- Args[idxPutYourTrimProgramDir_JAVA+1]  ## "/program/Trimmomatic/trimmomatic-0.39.jar"
PutYourFastQC_Program_Directory <- Args[idxPutYourFastQC_Program_Directory+1]
PutYourPutYourHISAT2_Pro_DIR <- Args[idxPutYourPutYourHISAT2_Pro_DIR+1]
PutYouSamtools_Pro_Dir <- Args[idxPutYouSamtools_Pro_Dir+1]
PutYourfeatureCounts_Pro_Dir <- Args[idxPutYourfeatureCounts_Pro_Dir+1]
ScientificName <- Args[idxScientificName+1] 
ReferenceID <- Args[idxReferenceID+1]
PutYourThreadsNumber <- Args[idxPutYourThreadsNumber+1]
PutYourTrimAdapter_Dir_ID <- Args[idxPutYourTrimAdapter_Dir_ID+1]

Ensembl_ReleaseVersion <- Args[idxEnsembl_ReleaseVersion+1]
PutYourTophredNumber <- Args[idxTophred+1]
PutSeedMismatches <- Args[idxMismatches+1] 
PutPalindromClipThreshold <- Args[idxPalinClipThresh+1] 
PutSimpleClipThreshold <- Args[idxSimpleClipThresh+1]
PutLEADING <- Args[idxLEADING+1] 
PutTRAILING <- Args[idxTRAILING+1]
PutWindowsize <- Args[idxWindowsize+1] 
PutReQuiredQuality <- Args[idxReQuiredQuality+1] 
PutMINLEN <- Args[idxMINLEN+1]
PutRNA_Strandness <- Args[idxRNA_Strandness+1]
PutStranded_FeatureCounts <- Args[idxStranded_Featurecounts+1]

## default 

Ensembl_ReleaseVersion <- 102
PutYourTophredNumber <- 33
PutSeedMismatches <- 2
PutPalindromClipThreshold <- 30
PutSimpleClipThreshold <- 10
PutLEADING <- 3
PutTRAILING <- 3
PutWindowsize <- 4
PutReQuiredQuality <- 15
PutMINLEN <- 36
PutRNA_Strandness <- ""
PutStranded_FeatureCounts <- 0 ## 0:unstranded 1:Stranded 2:Reverse Stranded 

if(length(idxEnsembl_ReleaseVersion)>0){
  Ensembl_ReleaseVersion <- Args[idxEnsembl_ReleaseVersion+1]
}
if(length(idxTophred)>0){
  PutYourTophredNumber <- Args[idxTophred+1]
}
if(length(idxMismatches)>0){
  PutSeedMismatches <- Args[idxMismatches+1] 
}
if(length(idxPalinClipThresh)>0){
  PutPalindromClipThreshold <- Args[idxPalinClipThresh+1] 
}
if(length(idxSimpleClipThresh)>0){
  PutSimpleClipThreshold <- Args[idxSimpleClipThresh+1]
}
if(length(idxLEADING)>0){
  PutLEADING <- Args[idxLEADING+1] 
}
if(length(idxTRAILING)>0){
  PutTRAILING <- Args[idxTRAILING+1]
}
if(length(idxWindowsize)>0){
  PutWindowsize <- Args[idxWindowsize+1] 
}
if(length(idxReQuiredQuality)>0){
  PutReQuiredQuality <- Args[idxReQuiredQuality+1] 
}
if(length(idxMINLEN)>0){
  PutMINLEN <- Args[idxMINLEN+1]
}
if(length(idxRNA_Strandness)>0){
  PutRNA_Strandness <- paste0(" --rna-strandness ",Args[idxRNA_Strandness+1])
}
if(length(idxStranded_Featurecounts)>0){
  PutStranded_FeatureCounts <- Args[idxStranded_Featurecounts+1]
}
sample_ID_Number <- length(list.files(path=PutYoutFastqDir,pattern="fastq.gz"))%/%2

CPUCMD <- paste0("top -n 1 -b| grep -i cpu\\(s\\)| awk \'{print $8}\'")

Current_CPU_remain <- system(CPUCMD,intern=T)
Current_CPU_remain_n <- as.numeric(Current_CPU_remain)

## Make Directory
MKdir <- c()
MKdir[1] <- paste0("mkdir ",RNA_seq_Directory,"/1.Reference_annotation")
MKdir[2] <- paste0("mkdir ",RNA_seq_Directory,"/2.Trimming_Trimmomatic")
MKdir[3] <- paste0("mkdir ",RNA_seq_Directory,"/3.QC_fastQC")
MKdir[4] <- paste0("mkdir ",RNA_seq_Directory,"/4.Mapping_HISAT2")
MKdir[5] <- paste0("mkdir ",RNA_seq_Directory,"/5.Sorting_SAMtools")
MKdir[6] <- paste0("mkdir ",RNA_seq_Directory,"/6.Quantification_featureCounts")
MKdir[7] <- paste0("mkdir ",RNA_seq_Directory,"/2.Trimming_Trimmomatic/Paired")
MKdir[8] <- paste0("mkdir ",RNA_seq_Directory,"/2.Trimming_Trimmomatic/UnPaired")
MKdir[9] <- paste0("mkdir ",RNA_seq_Directory,"/4.Mapping_HISAT2/INDEX")
MKdir[10] <- paste0("mkdir ",RNA_seq_Directory,"/5.Sorting_SAMtools/Sorted")

for(u in 1:length(MKdir)){
  system(MKdir[u])
}

Download_directory <- paste0(RNA_seq_Directory,"/1.Reference_annotation")
Trimming_directory <- paste0(RNA_seq_Directory,"/2.Trimming_Trimmomatic")
QC_directory <- paste0(RNA_seq_Directory,"/3.QC_fastQC")
Mapping_directory <- paste0(RNA_seq_Directory,"/4.Mapping_HISAT2")
Sorting_directory <- paste0(RNA_seq_Directory,"/5.Sorting_SAMtools")
Quantification_directory <- paste0(RNA_seq_Directory,"/6.Quantification_featureCounts")
PutYourTrimmoOutputDir_Pair <- paste0(RNA_seq_Directory,"/2.Trimming_Trimmomatic/Paired")
PutYourTrimmoOutputDir_UnPair <- paste0(RNA_seq_Directory,"/2.Trimming_Trimmomatic/UnPaired") 
Mk_INDEX_File_DIR <- paste0(RNA_seq_Directory,"/4.Mapping_HISAT2/INDEX")
PutYourSamtools_Sort_Dir <- paste0(RNA_seq_Directory,"/5.Sorting_SAMtools/Sorted")

###################################################################################
## 2. Trimmomatic Process  "Paried Ends", CROP, HEADCROP, MAXINFO are fixed   #####
###################################################################################
adapter <- system(paste0("find ",PutYourTrimAdapter_Dir_ID),intern=T)

setwd(Trimming_directory)

if(Current_CPU_remain_n < 50){
  PutYourThreadsNumber <- 16%/%sample_ID_Number
  print("Threads number is changed : 16")
}else{
  PutYourThreadsNumber <- 32%/%sample_ID_Number
  print("Threads number is unchanged : 32")
}

library(bracer)

if(length(list.files(path=PutYourTrimmoOutputDir_Pair,pattern=".paired.fastq.gz"))!=0){
  print("Data is existed ... Skip Trimming process")
}else{
  if(length(list.files(path=PutYoutFastqDir,pattern=".fastq.gz"))!=0 && length(adapter)!= 0){
    
    FastqFile <- list.files(path=PutYoutFastqDir,pattern="_1.fastq.gz")
    fastq_file_length <- length(list.files(path=PutYoutFastqDir,pattern="_1.fastq.gz"))
    
    id_r <- strsplit(FastqFile,"_")
    id_unlist <- as.data.frame(id_r)
    id <- c()
    id_left <- c()
    id_right<- c()
    id_left_output_paired<- c()
    id_right_output_paired<- c()
    id_left_output_unpaired<- c()
    id_right_output_unpaired<- c()
    TrimmomaticCMD <- c()
    
    for(i in 1:fastq_file_length){
      id[i] <- as.character(id_unlist[1,i])
      
      id_left[i] <-  paste0(id[i],"_1.fastq.gz")
      id_right[i] <-  paste0(id[i],"_2.fastq.gz") 
      
      id_left_output_paired[i] <-  paste0(id[i],"_1.paired.fastq.gz")  
      id_right_output_paired[i] <- paste0(id[i],"_2.paired.fastq.gz") 
      
      id_left_output_unpaired[i] <-  paste0(id[i],"_1.unpaired.fastq.gz")  
      id_right_output_unpaired[i] <-  paste0(id[i],"_2.unpaired.fastq.gz") 
    }
    
    for(k in 1:length(id)){
      TrimmomaticCMD[k] <- paste0("java -jar ",PutYourTrimProgramDir_JAVA," PE -threads ",PutYourThreadsNumber," -phred",PutYourTophredNumber," ",PutYoutFastqDir,"/",id_left[k]," ",PutYoutFastqDir,"/",id_right[k]," ",PutYourTrimmoOutputDir_Pair,"/",id_left_output_paired[k]," ",PutYourTrimmoOutputDir_UnPair,"/",id_left_output_unpaired[k]," ",PutYourTrimmoOutputDir_Pair,"/",id_right_output_paired[k]," ",PutYourTrimmoOutputDir_UnPair,"/",id_right_output_unpaired[k]," ","ILLUMINACLIP:",PutYourTrimAdapter_Dir_ID,":",PutSeedMismatches,":",PutPalindromClipThreshold,":",PutSimpleClipThreshold," ","LEADING:",PutLEADING," ","TRAILING:",PutTRAILING," ","SLIDINGWINDOW:",PutWindowsize,":",PutReQuiredQuality," ","MINLEN:",PutMINLEN," ","2> ",Trimming_directory,"/",id[k],".log &")
      TrimmomaticCMD_r <- as.data.frame(TrimmomaticCMD)
      TrimProcess <- as.character(TrimmomaticCMD_r[k,1])
      system(TrimProcess)
    }
    while(TRUE){
      if(length(system("ps -ef | grep 'java' | grep -v 'grep' | grep -v 'Darknight' | awk '{print $2}'",intern=T))==0)
        break
    }
    print("End Trimming processing")
    print("Start fastQC")
  }else{
    print("Error in your inputs for processing Trimmomatic")
    print("Please check your inputs & directories & files related to trimming")
    stop()
  }
}

## 3. Quality Control by fastQC 

if(Current_CPU_remain_n < 50){
  PutYourThreadsNumber <- 16%/%sample_ID_Number
  print("Threads number is changed : 16")
}else{
  PutYourThreadsNumber <- 32%/%sample_ID_Number
  print("Threads number is unchanged : 32")
}

if(length(list.files(path=QC_directory,pattern=".html"))!=0){
  print("Data is existed ... Skip fastQC process")
  
}else{
  if(length(list.files(path=PutYourTrimmoOutputDir_Pair,pattern=".paired.fastq.gz"))!=0){
    Fastq_QC_File <- list.files(path=PutYourTrimmoOutputDir_Pair,pattern="_1.paired.fastq.gz")
    Fastq_QC_file_length <- length(list.files(path=PutYourTrimmoOutputDir_Pair,pattern="_1.paired.fastq.gz"))
    
    QC_id_r <- strsplit(Fastq_QC_File,"_")
    QC_id_unlist <- as.data.frame(QC_id_r)
    QC_id <- c()
    QC_id_left <- c()
    QC_id_right <- c()
    
    for(a in 1:Fastq_QC_file_length){
      QC_id[a] <- as.character(QC_id_unlist[1,a])
      
      QC_id_left[a] <-  paste0(QC_id[a],"_1.paired.fastq.gz")
      QC_id_right[a] <-  paste0(QC_id[a],"_2.paired.fastq.gz") 
    }
    FastQC_CMD <- c()
    for(b in 1:length(QC_id)){
      FastQC_CMD[b] <- paste0(PutYourFastQC_Program_Directory," ","-o ",QC_directory," ","--noextract"," -f fastq -t ",PutYourThreadsNumber," ",PutYourTrimmoOutputDir_Pair,"/",QC_id_left[b]," ",PutYourTrimmoOutputDir_Pair,"/",QC_id_right[b]," &")
      system(FastQC_CMD[b])
    }
    while(TRUE){
      if(length(system("ps -ef | grep 'java' | grep -v 'grep' | grep -v 'Darknight' | awk '{print $2}'",intern=T))==0)
        break
    }
    print("End FastQC process")
    print("Start HISAT2 ")
  }else{
    print("Error in your inputs for processing FastQC")
    print("Please check your inputs & directories & files")
    stop()
  }
}

### 4. HISAT2 & Samtools

# (2) mapping & sam -> Sorted.bam

if(Current_CPU_remain_n < 50){
  PutYourThreadsNumber <- 16%/%sample_ID_Number
  print("Threads number is changed : 16")
}else{
  PutYourThreadsNumber <- 32%/%sample_ID_Number
  print("Threads number is unchanged : 32")
}

if(length(list.files(path=PutYourSamtools_Sort_Dir,pattern =".Sorted.bam"))!=0){
  print("Data is existed ... Skip mapping & sorting")
}else{
  if(length(list.files(path=PutYourTrimmoOutputDir_Pair,pattern=".paired.fastq.gz"))!=0 && length(list.files(path=Mk_INDEX_File_DIR,pattern=".ht2"))!=0 ){
    setwd(Mapping_directory)
    
    Hisat2_File <- list.files(path=PutYourTrimmoOutputDir_Pair,pattern="_1.paired.fastq.gz")
    Hisat2_file_length <- length(list.files(path=PutYourTrimmoOutputDir_Pair,pattern="_1.paired.fastq.gz"))
    Hisat2_id_r <- strsplit(Hisat2_File,"_")
    Hisat2_id_unlist <- as.data.frame(Hisat2_id_r)
    Hisat2_id <- c()
    Hisat2_id_left <- c()
    Hisat2_id_right <- c()
    
    for(c in 1:Hisat2_file_length){
      Hisat2_id[c] <- as.character(Hisat2_id_unlist[1,c])
      Hisat2_id_left[c] <-  paste0(Hisat2_id[c],"_1.paired.fastq.gz")
      Hisat2_id_right[c] <- paste0(Hisat2_id[c],"_2.paired.fastq.gz") 
    }
    Mapping_SAM_BAM_CMD <- c()
    for(d in 1:length(Hisat2_id)){
      Mapping_SAM_BAM_CMD[d] <- paste0(PutYourPutYourHISAT2_Pro_DIR,"/hisat2 -p ",PutYourThreadsNumber,PutRNA_Strandness," -x ",Mapping_directory,"/INDEX/IndexName -1 ",PutYourTrimmoOutputDir_Pair,"/",Hisat2_id_left[d]," -2 ",PutYourTrimmoOutputDir_Pair,"/",Hisat2_id_right[d]," 2> ",Hisat2_id[d],".log | ",PutYouSamtools_Pro_Dir," sort -@ ",PutYourThreadsNumber," -o ",PutYourSamtools_Sort_Dir,"/",Hisat2_id[d],".Sorted.bam &")
      system(Mapping_SAM_BAM_CMD[d])
    }
    hisat2_flag=0
    samtools_flag=0
    while(TRUE){
      if(length(system("ps -ef | grep 'hisat2-align-s' | grep -v 'grep' | grep -v 'practice' | awk '{print $2}'",intern=T))==0){
        hisat2_flag=1
      }
      if(length(system("ps -ef | grep 'samtools' | grep -v 'grep' | grep -v 'Sbo' | grep -v 'SCREEN' | grep -v usr | awk '{print $2}'",intern=T))==0){
        samtools_flag=1
      }
      if(hisat2_flag==1 && samtools_flag ==1){
        break
      }
    }
    print("End mapping process")
    print("Start featureCounts")
  }else{
    print("Error in your inputs for mapping & sam -> Sorted.bam")
    print("Please check your inputs & directories & files")
    stop()
  }
}

### 5. featureCounts

if(Current_CPU_remain_n < 50){
  PutYourThreadsNumber <- 16%/%sample_ID_Number
  print("Threads number is changed : 16")
}else{
  PutYourThreadsNumber <- 32%/%sample_ID_Number
  print("Threads number is unchanged : 32")
}

if(length(list.files(path=Quantification_directory, pattern =".txt.summary"))!=0){
  print("Data is existed ... Skip FeatureCounts process")
  print("Darknight is over ..")
}else{
  if(length(list.files(path=PutYourSamtools_Sort_Dir, pattern =".Sorted.bam"))!=0){
    setwd(Quantification_directory)
    featureCounts_File <- list.files(path=PutYourSamtools_Sort_Dir,pattern=".Sorted.bam")
    featureCounts_file_length <- length(list.files(path=PutYourSamtools_Sort_Dir,pattern=".Sorted.bam"))
    featureCounts_id_r <- strsplit(featureCounts_File,"[.]")
    featureCounts_id_unlist <- as.data.frame(featureCounts_id_r)
    featureCounts_id <- c()
    featureCounts_CMD <- c()
    
    for(q in 1:featureCounts_file_length){
      featureCounts_id[q] <- as.character(featureCounts_id_unlist[1,q])
    }
    for(w in 1:length(featureCounts_id)){
      featureCounts_CMD[w]<- paste0(PutYourfeatureCounts_Pro_Dir," -T ",PutYourThreadsNumber," -s ",PutStranded_FeatureCounts," -p -t exon -g gene_id -a ",Download_directory,"/",ScientificName,".",ReferenceID,".",Ensembl_ReleaseVersion,".gtf"," -o ",Quantification_directory,"/", featureCounts_id[w],"_counts.txt ",PutYourSamtools_Sort_Dir,"/",featureCounts_id[w],".Sorted.bam &")
      system(featureCounts_CMD[w])
    }
    while(TRUE){
      if(length(system("ps -ef | grep 'featureCounts' | grep -v 'grep' | grep -v 'usr' | awk '{print $2}'",intern=T))==0){
        break
      }
    }
    print("End featureCounts process")
    print(" End all process, check your output files")
  }else{
    print("Error in your inputs for featureCounts")
    print("Please check your featureCounts")
    stop()
  }
}