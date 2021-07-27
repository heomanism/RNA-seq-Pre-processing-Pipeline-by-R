# RNA-seq-Pre-processing Pipeline by R language

Manual for Darknight made by Min Heo ^^ 

Inputs' index, if index is "-a", it is essential index, and if index is "--a", it is not essential, it has default value. If you want to change these values, you can change them.

## Index 

### Necessary Index 
-Drna : Put your RNA-seq pipeline directory which you want to store ouputs. 

-Dfq: Put your directory that includes fastq files(sample files).

-Dtri: Put your directory that includes Trimmomatic program.

-Dqc: Put your directory that includes fastQC program.

-Dhi: Put your directory that includes Hisat2 program.

-Dsm: Put your directory that includes Samtools program.

-Dfc: Put your directory that includes featureCounts program.

-sn: Put Scientific name of your data(fastqc) for your Annotation & Reference. ex)Homo_sapiens  <- first letter should be capitalized.

-rid: Put Reference ID of your data. Please check ID in Ensembl site.

-t: Put your Threads numbers, but if your CPU usage is more than 50%, the threads numbers will be changed. (32 -> 16, 32 is maximum number.)

-a: Put your Adapter Directory & ID. 


### Not Necessary Index (Put options if you need these. () is default value.)
--r: Ensembl Release version (102)

-p: Tophred number (33)

--m: Trimmomatic's mismatch numbers (2)

--pthr: Trimmomatic's Palindrome Clip Threshold (30)

--sthr: Trimmomatic's Simple Clip Threshold (10)

--l: Trimmomatic's Leading (3)

--t: Trimmomatic's Trailing (3)

--w: Trimmomatic's Window Sliding Size (4)

--q: Trimmomatic's requiredQuality (15)

--ml: Trimmomatic's MINLENGTH (36)

--rna-strandness: Hisat2's rna-strandness option (unstranded) others: RF FR

--s: FeatureCount's stranded option (0:unstranded) others: 1:stranded, 2:reverse stranded 


### Example for inputs 

RNA_seq_Directory <- "/disk11/2.Pipeline_test_hm/1.Test"    (Put your RNA-seq process Directory) /disk1/9.PipelineTest/3.MinHeo/V2.Test

PutYoutFastqDir <- "/disk11/2.Pipeline_test_hm/4.RNA-seq_pipeline_HM/0.fastq/sample" (Put your fastq Directory) /disk1/9.PipelineTest/3.MinHeo

PutYourTrimProgramDir_JAVA <- "/program/Trimmomatic/trimmomatic-0.39.jar"  (Put your Trimmomatic Directory)

PutYourFastQC_Program_Directory <- "/program/FastQC/fastqc" (Put your fastqc Directory)

PutYourPutYourHISAT2_Pro_DIR <- "/program/HISAT2" (Put your HISAT2 Directory)

PutYouSamtools_Pro_Dir <-"/program/samtools/bin/samtools" (Put your SAMtools Directory)

PutYourfeatureCounts_Pro_Dir <-"/program/subread/bin/featureCounts" (Put your featureCounts Directory)

ScientificName <- "Coturnix_japonica"  #First letter is must be capitalized. 

ReferenceID <- "Coturnix_japonica_2.0" # You can find in Ensembl

PutYourThreadsNumber <- 32 

PutYourTrimAdapter_Dir_ID <- "/program/Trimmomatic/adapters/TruSeq3-PE-2.fa" ##Illumina version only

### Code example 
time Rscript Darknight.v.1.1.R -Drna "/disk1/9.PipelineTest/3.MinHeo/V3.Test" -Dfq "/disk1/9.PipelineTest/data" -Dtri "/program/Trimmomatic/trimmomatic-0.39.jar" -Dqc "/program/FastQC/fastqc" -Dhi "/program/HISAT2" -Dsm "/program/samtools/bin/samtools" -Dfc "/program/subread/bin/featureCounts" -sn "Coturnix_japonica" -rid "Coturnix_japonica_2.0" -t 32 -a "/program/Trimmomatic/adapters/TruSeq3-PE-2.fa" --s 0

