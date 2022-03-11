rm(list = ls())
library(rtracklayer) ## for liftOver
library(vcfR)        ## to import 1000 Genomes data
library(beepr)
library(biomaRt)

varlist_save <- c()

dir <- "~/Downloads/" # change the directory to where you downloaded the data

#################################### GET VARIANT AGE DATA ####################################
## let's import the data from human.genome.dating
## downloaded at: https://human.genome.dating/gene/EDAR
EDAR = read.csv(paste0(dir,"atlas.EDAR.gene.csv"),skip = 7, header = T) # skip 7 because the first lines are a preamble.
pos = EDAR[,3]
pos_unique = unique(pos)
## note: the position of these variants is in the format GRCh37
## we need to convert them to GRCh38

# let's get the age of the variants using the joint model. 
# We need to be careful because some variants are repeated 
# and for those we need to use the combined estimate
from_pos_to_unique = match(pos, pos_unique)
from_unique_to_pos = match(pos_unique, pos)
from_unique_to_pos_adj <- from_unique_to_pos

combined <- which(EDAR[,"DataSource"] == " Combined")
combined2 <- c()
for(x in 1:(length(from_unique_to_pos)-1)){
  if(from_unique_to_pos[x]+1 != from_unique_to_pos[x+1]){
    from_unique_to_pos_adj[x] = from_unique_to_pos[x]+2
    combined2 <- c(combined2, from_unique_to_pos[x]+2)
  }
}
# if the last one is not the last row, then probably (CHECK IN THE FUTURE!) there is another "Combined"
if(from_unique_to_pos[length(from_unique_to_pos)] != nrow(EDAR)){
  x = length(from_unique_to_pos)
  from_unique_to_pos_adj[x] = from_unique_to_pos[x]+2
  combined2 <- c(combined2, from_unique_to_pos[x]+2)
}
# all.equal(combined2,combined)
## in this case we have something weird, because "combined" appears after 1 instead of 2
which(combined2-combined != 0) #  248 279
combined2[c(248,279)]
which(from_unique_to_pos_adj == combined2[248]) # 682
which(from_unique_to_pos_adj == combined2[279]) # 770
from_unique_to_pos_adj[682] <- combined[248]
from_unique_to_pos_adj[770] <- combined[279]
## how can we do this better?



index_EDAR = from_unique_to_pos_adj

agemode = EDAR[index_EDAR,"AgeMode_Jnt"]
agemedian = EDAR[index_EDAR,"AgeMedian_Jnt"]
agemean = EDAR[index_EDAR,"AgeMean_Jnt"]
# hist(agemode[agemode < 10000], breaks = 100)
## use the MODE

## we can also save a subset of the original EDAR
EDAR_combined <- EDAR[index_EDAR,]  ## ORDERED by pos and pos38

varlist_save <- c(varlist_save, "EDAR", "pos", "pos_unique", 
                  "index_EDAR","agemode","EDAR_combined")



# The chain file gives the "conversion key". I got the chain file from
# https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/
ch = import.chain(paste0(dir,"hg19ToHg38.over.chain/hg19ToHg38.over.chain"))

# to use the function liftOver, we need the `positions` to be in format GRanges
gr <- GRanges(seqnames=Rle(c("chr2"), length(pos_unique)), ## gene EDAR is on chromosome 2
              ranges=IRanges(pos_unique, pos_unique),			 ## since it's SNP, the start and end coincide
              strand=rep(c("+"), length(pos_unique)))
gr
pos38 = liftOver(gr, ch)
pos38 = pos38@unlistData@ranges@start

varlist_save <- c(varlist_save, "pos38")

#################################### GET SAMPLE DATA FROM 1000 GENOME ####################################
## Now we need to import the data from 1000 Genomes
## from their data portal: https://www.internationalgenome.org/data-portal/sample 
## click on any sample and select "1000 Genomes on GRCh38" - then find the right chromosome
## Note that we had to decompress it first (using 7-zip on windows or gunzip on mac)
## note: I downloaded the file using internet explorer (on windows) or safari (on mac)
time1 = Sys.time()
chr2 <- read.vcfR(file = paste0(dir,"ALL.chr2.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf"),
                  verbose = TRUE,convertNA = TRUE, cols = 1:3)
time2 = Sys.time() - time1
beepr::beep()

# we need to look for the right position to start importing the data:
pos1 <- chr2@fix[,2]
pos1 <- as.numeric(pos1)
range(pos1) # 10103 242171544

## EDAR on Ensemble (GRCh38)
## Chromosome 13: 108894471-108989372  
## > range(pos38)
## [1] 108,894,472 108,989,311

## after some trial and error to tune skip and nrows:
time1 = Sys.time()
chr2B <- read.vcfR(file = paste0(dir,"ALL.chr2.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf"),
                   verbose = TRUE, skip = 3096000, nrows = 100000, convertNA = TRUE) # you can try with nrows = 4000
time2 = Sys.time() - time1
beepr::beep()

pos2 <- chr2B@fix[,2]
pos2 <- as.numeric(pos2)
range(pos2, na.rm = T)  #  105,409,830 142,067,265

i_min = min(which(pos2 >= 108894472))
i_max = max(which(pos2 <= 108989311))

index_chr2B <- seq(i_min,i_max)
Gen1000_fix = chr2B@fix[index_chr2B,-8]    ## NOT ORDERED by pos38 (higher dimensionality)
Gen1000_gt = chr2B@gt[index_chr2B,-1]

Gen1000_fix = as.data.frame(Gen1000_fix)

pos_Gen1000 = as.numeric(Gen1000_fix[,2])
# index_ensembl = match(pos38,pos_Gen1000)
from_dating_to_ensembl = match(pos38,pos_Gen1000)
from_ensembl_to_dating = match(pos_Gen1000,pos38)

# these are the same, except for NA's
# head(pos_Gen1000[match(pos38,pos_Gen1000)])
# head(pos38)

sum(is.na(from_dating_to_ensembl)) # 102 dated species are not found in ensembl -- NOTE  I ran it a second time and got 109 na.
head(Gen1000_fix[from_dating_to_ensembl,])

from_dating_to_ensembl_noNA <- from_dating_to_ensembl[!is.na(from_dating_to_ensembl)]

varlist_save <- c(varlist_save, "Gen1000_fix","Gen1000_gt",
                  "from_dating_to_ensembl","from_ensembl_to_dating", "from_dating_to_ensembl_noNA")

Gen1000_fix$AGE_MODE = NA
Gen1000_fix$AGE_MODE[from_dating_to_ensembl_noNA] = agemode[!is.na(from_dating_to_ensembl)]

Gen1000_fix$ID_dating <- NA
Gen1000_fix$ID_dating[from_dating_to_ensembl_noNA] = EDAR_combined[,"VariantID"][!is.na(from_dating_to_ensembl)]

Gen1000_fix$POS_dating_GRCh37 <- NA
Gen1000_fix$POS_dating_GRCh37[from_dating_to_ensembl_noNA] = EDAR_combined[,"Position"][!is.na(from_dating_to_ensembl)]

# we could potentially save a dataset, but it's NOT ordered by pos38 because the NA were removed.
EDAR_combined_pos38NA <- EDAR_combined[!is.na(from_dating_to_ensembl),]

# this is to check that what we did is right
# as.numeric(Gen1000_fix$POS[from_dating_to_ensembl])[1:10] # pos38[1:10] # visually comparing these two
# all.equal(as.numeric(Gen1000_fix$POS[from_dating_to_ensembl_noNA]),pos38[!is.na(from_dating_to_ensembl)])

#################################### CHECK THE VARIANTS ID IF THEY MATCH ####################################
## we tried to use biomart to check for the IDs but it times out. 
## even with very small position range it times out, so I have downloaded from the website
## from https://www.ensembl.org/biomart/martview/ 
## select Database: Ensembl Variation 105 - select Dataset: Human Short Variants (SNPs and indels excluding flagged variants)
## add filter: GENE ASSOCIATED VARIANT FILTERS > Gene stable ID(s) <- ENSG00000139618 
## ps. this ID is found on the BRCA2 page of Ensembl
## keep the standard attributes: 
## Variant name; Variant source; Chromosome/scaffold name; Chromosome/scaffold position start (bp); Chromosome/scaffold position end (bp)
## finally click "results" in the top left, and export all results to compressed (.gz) CSV [unique results]

biomart_data = read.csv(paste0(dir, "martquery_EDARunique.txt"), header =TRUE, sep = ",")
colnames(biomart_data) <- c("variant_name","source","chrom","start","end")
dim(biomart_data)
mean(pos38 %in% biomart_data[,"start"]) # 1.00!!

length(unique(biomart_data$variant_name)) # 24688
length(unique(biomart_data$start)) # 23980
length(unique(biomart_data$end))   # 23954

biomart_data$start_eq_end <- sapply(1:nrow(biomart_data),FUN=function(x)biomart_data$start[x]==biomart_data$end[x])
biomart_data$start_end <- NA
for(x in which(biomart_data$start_eq_end)){
  biomart_data$start_end[x] <- biomart_data$start[x]
}

# from_pos38_to_mart_old = match(pos38,biomart_data[,"start"])
from_pos38_to_mart = match(pos38,biomart_data[,"start_end"])
# all.equal(biomart_data[from_pos38_to_mart,"start"],pos38) #TRUE
# biomart_data[from_pos38_to_mart,"variant_name"]
# sum(is.na(biomart_data[from_pos38_to_mart,"variant_name"]))  # 0

# note: pos_unique has the same order of EDAR_combined and
# pos38 is matched to pos_unique. So we can do
EDAR_combined$VariantID_mart = biomart_data[from_pos38_to_mart,"variant_name"]
# let's save this also in Gen1000_fix
Gen1000_fix$ID_mart <- NA
Gen1000_fix$ID_mart[from_dating_to_ensembl_noNA] = EDAR_combined[,"VariantID_mart"][!is.na(from_dating_to_ensembl)]

Gen1000_fix$POSstart_mart <- NA
Gen1000_fix$POSend_mart <- NA
Gen1000_fix$POSstart_mart[from_dating_to_ensembl_noNA] = biomart_data[from_pos38_to_mart,"start"][!is.na(from_dating_to_ensembl)]
Gen1000_fix$POSend_mart[from_dating_to_ensembl_noNA] = biomart_data[from_pos38_to_mart,"end"][!is.na(from_dating_to_ensembl)]

EDAR_combined[1:10,c("VariantID","VariantID_mart")] # promising!
all.equal(EDAR_combined[,"VariantID"],EDAR_combined[,"VariantID_mart"]) # "22 string mismatches"
all.equal(EDAR_combined[!is.na(from_dating_to_ensembl),"VariantID"],EDAR_combined[!is.na(from_dating_to_ensembl),"VariantID_mart"]) #  "5 string mismatches"

equal_name = apply(EDAR_combined[,c("VariantID","VariantID_mart")],MARGIN = 1, FUN = function(x) x[1]==x[2])
names(equal_name) = NULL

Gen1000_fix$Equal_name = NA
Gen1000_fix$Equal_name[from_dating_to_ensembl_noNA] = equal_name[!is.na(from_dating_to_ensembl)]

EDAR_combined[!equal_name,c("VariantID","VariantID_mart")]
# most of them have weird names, like X409552936 in VariantID and a rs... in VariantID_mart.
# those with normal name have same start and end, and exist in biomart with two names. Like:
biomart_data[which(biomart_data[,"variant_name"]=="rs143471570"),] 
biomart_data[which(biomart_data[,"variant_name"]=="rs955846975"),] 
# if you check https://www.ensembl.org/Homo_sapiens/Variation/Explore?r=2:108928827-108929827;v=rs143471570;vdb=variation;vf=250895712
# you see it has a "co-located variant" with the second name
# same below:
biomart_data[which(biomart_data[,"variant_name"]=="rs147059377"),] 
biomart_data[which(biomart_data[,"variant_name"]=="rs1558814135"),] 


# biomart_data[biomart_data[,"variant_name"] %in% EDAR_combined[!equal_name,"VariantID_mart"],]

## this is difficult to understand
match( pos38[match(EDAR_combined[!equal_name,c("Position")],pos_unique)],
       Gen1000_fix[,"POS"])
# this is exactly the same as before, but maybe easier?
from_dating_to_ensembl[match(EDAR_combined[!equal_name,c("Position")],pos_unique)]
# ie. most of the different names also are missing in the ensemble data in terms of position

# the ones that are present in ensembl (Gen1000_fix)
# noequalname_inensembl = which(!is.na(from_dating_to_ensembl[match(EDAR_combined[!equal_name,c("Position")],pos_unique)]))
## note which(!equal_name) is the same as match(EDAR_combined[!equal_name,c("Position")],pos_unique)
# all.equal(match(EDAR_combined[!equal_name,c("Position")],pos_unique), which(!equal_name) )
# noequalname_inensembl = which(!is.na(from_dating_to_ensembl[ !equal_name ])) ## wrong
noequalname_inensembl = which(  (!is.na(from_dating_to_ensembl))  &  (!equal_name)  )

EDAR_combined[ noequalname_inensembl ,c("VariantID","VariantID_mart")] # different names but both rs...
# while the others are missing
biomart_data[from_pos38_to_mart[noequalname_inensembl],]

varlist_save <- c(varlist_save, "biomart_data","from_pos38_to_mart", "equal_name",
                  "noequalname_inensembl")
save(list = varlist_save,file= "data/output_cleandata_EDAR.RData")


Gen1000_fix_keep = Gen1000_fix[which(Gen1000_fix$Equal_name == TRUE),]
sum(is.na(Gen1000_fix_keep[,-c(3,6)]))
Gen1000_fix_keep = Gen1000_fix_keep[,-c(3,6)]

Gen1000_gt_keep = Gen1000_gt[which(Gen1000_fix$Equal_name == TRUE),]
Gen1000_gt_keep01 = ifelse(Gen1000_gt_keep == "0|0",0,1)
tot_count = rowSums(Gen1000_gt_keep01)

table(tot_count)

varlist_save <- c(varlist_save, "Gen1000_fix_keep","Gen1000_gt_keep", "Gen1000_gt_keep01",
                  "tot_count")
save(list = varlist_save,file= "data/output_cleandata_EDAR.RData")


tmp = load("data/output_cleandata_EDAR.RData")

head(Gen1000_fix_keep$AGE_MODE)
Gen1000_fix_keep$rank <- NA

# this would give the larges (in terms of large N) rank to the oldest (large age)
# and rank 1 to the youngest (small age)
# match(Gen1000_fix_keep$AGE_MODE[1:10],sort(Gen1000_fix_keep$AGE_MODE[1:10], decreasing = FALSE))
# cbind(1:10,Gen1000_fix_keep$AGE_MODE[1:10])
Gen1000_fix_keep$rank <- match(Gen1000_fix_keep$AGE_MODE,
                               sort(Gen1000_fix_keep$AGE_MODE, decreasing = FALSE))

which.max(Gen1000_fix_keep$rank) # 1311
Gen1000_fix_keep[1311,]       # rank 1482

tot_count[1311] # 1014 # not the largest cluster! (2548)
# there are other 139 species with more than this frequency 
# and the most frequent species has 2548 observations (i.e. all samples).
which(tot_count == max(tot_count))
Gen1000_fix_keep[which(tot_count == max(tot_count)),] 
# the most frequent are other very old species (rank 1456-1468 instead of 1482 and age 55K/63K instead of 123K).

save(list = c("Gen1000_fix_keep","Gen1000_gt_keep", "Gen1000_gt_keep01"),
     file = "data/Gen1000_cleaned_EDAR.Rdata")

