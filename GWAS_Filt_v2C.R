# From 2B,


### Load Packages ###############################################################################
# install.packages("stringi")
library(stringi)

# BiocManager::install("liftOver")
library(liftOver)


# BiocManager::install("ensemblVEP")
#library(ensemblVEP)


# BiocManager::install("VariantAnnotation")
library(VariantAnnotation)

#install.packages('annovarR')
#library(annovarR)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#install.packages('vcfR')
#library(vcfR)

#BiocManager::install("TVTB")
library(TVTB)




### Specify directory/file ###############################################################################

drt0 = "I:\\Daniel\\Projects\\IGVF\\GWAS2\\"
wDrt0 = paste(drt0, "Out20220824\\", sep = "")

pThresh = 2

saveSig = 1
runId = "20220824"
pThreshAnno = log(5e-8)


### Load data / initialize ###############################################################################

dir.create(wDrt0,showWarnings = FALSE)


# For liftover
path = system.file(package = "liftOver", "extdata", "hg19ToHg38.over.chain")
ch = import.chain(path)


# Load ClinVar
drt1 = "I:\\Daniel\\Datasets\\SNPs\\"
file0 = "variant_summary_20220624.csv"
fID = paste(drt1, file0, sep = "")


# Get headers and format for read in
headers0 = read.table(fID, header = T, nrows = 2,sep=',')
headers0 = names(headers0)


nColumns = length(headers0)

cols = rep("NULL",  nColumns)
cols[c( 8,9,10,11,13)] = "character"



clinVar0 = read.csv(fID, header = T,colClasses = cols)


ii = is.element(clinVar0$Assembly,"GRCh38")
clinVar0 = clinVar0[ii,]
clinVar0$Assembly = NULL

clinVar0$pos = paste(paste("chr",clinVar0$Chromosome,sep=''),clinVar0$Start,sep="_")

clinVar0$Chromosome = NULL
clinVar0$Start = NULL

unique(clinVar0$ReviewStatus)
ii = is.element(clinVar0$ReviewStatus,c("criteria provided, single submitter","criteria provided, multiple submitters, no conflicts","practice guideline","criteria provided, conflicting interpretations","reviewed by expert panel") )
clinVar1 = clinVar0[ii,]

# Read annotations
drt2 = "I:\\Daniel\\Projects\\IGVF\\GWAS\\"
file2 = "gwasSnpAnno20220623v1.tsv"
file2 = paste(drt2, file2, sep = "")

nColumns = 2
cols = rep("character",  nColumns)

anno0 = read.table(file2, header = F, colClasses = cols,sep = "\t")





### Run (loop)--Pos + Meta P/Beta/AF with Liftover + Annotations ###############################################################################

###Get files


allFiles = list.files(drt0)
n = seq(from = 1,
        to = length(allFiles),
        by = 1)

vOut = rep(FALSE, length(allFiles))
for (j in n) {
  v0 = allFiles[j]
  
  v1 = substr(v0, nchar(v0) - 2, nchar(v0))
  
  if (stri_cmp_eq("bgz", v1)) {
    vOut[j] = TRUE
  }
}
allFiles = allFiles[vOut]

###Full loop

v00= c('coding','spliceSite','threeUTR','fiveUTR','intron','promoter','NonCoding')

nTotal = length(allFiles)
nFiles = seq(from = 1,
             to = nTotal,
             by = 1)



# rm(list = c("vAllOut","vAllOut2"))
for (i in nFiles) {
  
  file0 = allFiles[i]
 
  
  file00 = substr(file0, 1,nchar(file0) - 8)
  file1 = paste(file00, "_hg38.tsv", sep = "")
  file1 = paste(drt0, file1, sep = "")
  print(paste(i,"of",nTotal, file00))
  if (!file.exists(file1)){
    # Make file path and open
    fileFull = paste(drt0, file0, sep = "")
    fID = gzfile(fileFull)
    
    
    # Get headers and format for read in
    headers0 = read.table(fID, header = T, nrows = 1)
    headers0 = names(headers0)
    
    
    nColumns = length(headers0)
    
    cols = rep("NULL",  nColumns)
    cols[2] = "numeric"
    cols[c(1, 3, 4)] = "character"
    
    v0 = any(stri_cmp_eq("pval_meta", headers0))
    if (any(stri_cmp_eq("pval_meta", headers0))) {
      cols[stri_cmp_eq("pval_meta", headers0)] = "numeric"
      cols[stri_cmp_eq("beta_meta", headers0)] = "numeric"
      cols[stri_cmp_eq("af_cases_meta", headers0)] = "numeric"
      cols[stri_cmp_eq("af_controls_meta", headers0)] = "numeric"
      useMeta = 1
      
    } else {
      cols[stri_cmp_eq("pval_EUR", headers0)] = "numeric"
      cols[stri_cmp_eq("beta_EUR", headers0)] = "numeric"
      cols[stri_cmp_eq("af_cases_EUR", headers0)] = "numeric"
      cols[stri_cmp_eq("af_controls_EUR", headers0)] = "numeric"
      useMeta = 0
    }
    
    
    
    
    # Read file in
    fID = gzfile(fileFull)
    dat0 = read.table(fID, header = T, colClasses = cols)
    
    
    # Retain if P-value exists
    if (any(stri_cmp_eq("pval_meta", headers0))) {
      idcs =  !is.na(dat0$pval_meta)
    } else {
      idcs =  !is.na(dat0$pval_EUR)
    }
    
    
    
    
    dat0 = dat0[idcs,]
    
    
    
    
    if (useMeta == 1) {
      gr <- GRanges(
        seqnames = dat0$chr,
        ranges = IRanges(dat0$pos, width = 1),
        strand = rep("*", nrow(dat0)),
        ref = dat0$ref,
        alt = dat0$alt,
        af_cases = dat0$af_cases_meta,
        af_controls = dat0$af_controls_meta,
        beta = dat0$beta_meta,
        pval = dat0$pval_meta
      )
    } else{
      gr <- GRanges(
        seqnames = dat0$chr,
        ranges = IRanges(dat0$pos, width = 1),
        strand = rep("*", nrow(dat0)),
        ref = dat0$ref,
        alt = dat0$alt,
        af_cases = dat0$af_cases_EUR,
        af_controls = dat0$af_controls_EUR,
        beta = dat0$beta_EUR,
        pval = dat0$pval_EUR
      )
    }
    
    
    if (pThresh<1){
      idcs = gr$pval<pThresh
      gr = gr[idcs,]
    }
    
    seqlevelsStyle(gr) = "UCSC"  # necessary
    gr = liftOver(gr, ch)
    gr = unlist(gr)
    
    
    
    if (exists("grAll")){
      grAll = union(grAll, gr[,1:2])
    }else {
      grAll =  gr[,1:2]
    }
    
    
    gr = data.frame(
      chr = as.factor(gr@seqnames),
      pos = as.factor(gr@ranges),
      ref = gr$ref,
      alt = gr$alt,
      af_cases = gr$af_cases,
      af_controls = gr$af_controls,
      beta = gr$beta,
      pval = gr$pval
    )
    
    
    
    # Export file as tsv
    file1 = substr(file0, 1, nchar(file0) - 8)
    file1 = paste(file1, "_hg38.tsv", sep = "")
    file1 = paste(drt0, file1, sep = "")
    
    
    write.table(gr,
                file1,
                quote = FALSE,
                sep = '\t',
                row.names = FALSE)
    
  } else {
    # Get headers and format for read in
    headers0 = read.table(file1, header = T, nrows = 1)
    headers0 = names(headers0)
    
    nColumns = length(headers0)
    cols = rep("NULL",  nColumns)
    cols[1] = "character"
    cols[2] = "character"
    # cols[6] = "numeric"
    cols[8] = "numeric"
    
    gr = read.table(file1, header = T, colClasses = cols)
    
  }
  
  

  ### Annotate 
  
  
  ii = gr$pval<pThreshAnno 
  gr = gr[ii,]
  
  
  gr$pos = paste(gr$chr, gr$pos, sep = "_")
  gr$chr = NULL
  
  # For all variants
  ii = is.element(anno0$V1,gr$pos )
  if(saveSig==1){
    v0 = anno0[ii,]
    file1 = paste(wDrt0, "SigPos_", file00,".csv", sep = "")
    write.table(v0,
                file1,
                quote = FALSE,
                sep = ',',
                row.names = FALSE)
  }
  v0 = anno0[ii,]$V2
  v0 = c(v0,v00)
  vS = table(v0)-1
  vS = vS[c(1,6,7,2,3,5,4) ]
  
  
  v0 = data.frame(
    vS = as.factor(vS)
  )
  names(v0)[names(v0) == "vS"] <- file00
  
  # rm(vAllOut)
  if (exists("vAllOut")){
    vAllOut = cbind(vAllOut, v0)
  }else{
    vAllOut = v0
  }
  
  
  
  file1 = paste("mergedAnnoCountsAll_",runId,".csv", sep = "")
  file1 = paste(wDrt0, file1, sep = "")
  
  write.table(vAllOut,
              file1,
              quote = FALSE,
              sep = ',',
              row.names = TRUE)
  
  
  
  # For ClinVar variants
  
  ii = is.element(gr$pos,clinVar1$pos )
  dat1 = gr[ii,]
  
  ii = is.element(anno0$V1,dat1$pos )
  v0 = anno0[ii,]$V2
  v0 = c(v0,v00)
  vS = table(v0)-1
  vS = vS[c(1,6,7,2,3,5,4) ]
  
  
  v0 = data.frame(
    vS = as.factor(vS)
  )
  names(v0)[names(v0) == "vS"] <- file00
  
  # rm(vAllOut2)
  if (exists("vAllOut2")){
    vAllOut2 = cbind(vAllOut2, v0)
  }else{
    vAllOut2 = v0
  }
  
  
  
  file1 = paste("mergedAnnoCountsClinVar_",runId,".csv", sep = "")
  file1 = paste(wDrt0, file1, sep = "")
  
  write.table(vAllOut2,
              file1,
              quote = FALSE,
              sep = ',',
              row.names = TRUE)
  
  
}







