library(rtracklayer) ## for liftOver
library(vcfR)        ## to import 1000 Genomes data
library(beepr)
library(biomaRt)

varlist_save <- c()

dir <- "~/Downloads/" # change the directory to where you downloaded the data

#################################### GET VARIANT AGE DATA ####################################
## let's import the data from human.genome.dating
## downloaded at: https://human.genome.dating/gene/BRCA2
BRCA2 = read.csv(paste0(dir,"atlas.BRCA2.gene.csv"),skip = 7, header = T) # skip 7 because the first lines are a preamble.
pos = BRCA2[,3]
pos_unique = unique(pos)
## note: the position of these variants is in the format GRCh37
## we need to convert them to GRCh38

# let's get the age of the variants using the joint model. 
# We need to be careful because some variants are repeated 
# and for those we need to use the combined estimate
from_pos_to_unique = match(pos, pos_unique)
from_unique_to_pos = match(pos_unique, pos)
from_unique_to_pos_adj <- from_unique_to_pos

combined <- which(BRCA2[,"DataSource"] == " Combined")
combined2 <- c()
for(x in 2:length(from_unique_to_pos)){
  if(!(from_unique_to_pos[x] == from_unique_to_pos[x-1]+1)){
    from_unique_to_pos_adj[x-1] = from_unique_to_pos[x-1]+2
    combined2 <- c(combined2, from_unique_to_pos[x-1]+2)
  }
}
# all.equal(combined2,combined) # TRUE

index_BRCA2 = from_unique_to_pos_adj

agemode = BRCA2[index_BRCA2,"AgeMode_Jnt"]
agemedian = BRCA2[index_BRCA2,"AgeMedian_Jnt"]
agemean = BRCA2[index_BRCA2,"AgeMean_Jnt"]
# hist(agemode[agemode < 10000], breaks = 100)
## use the MODE

## we can also save a subset of the original BRCA2
BRCA2_combined <- BRCA2[index_BRCA2,]  ## ORDERED by pos and pos38

varlist_save <- c(varlist_save, "BRCA2", "pos", "pos_unique", 
                  "index_BRCA2","agemode","BRCA2_combined")



# The chain file gives the "conversion key". I got the chain file from
# https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/
ch = import.chain(paste0(dir,"hg19ToHg38.over.chain/hg19ToHg38.over.chain"))

# to use the function liftOver, we need the `positions` to be in format GRanges
gr <- GRanges(seqnames=Rle(c("chr13"), length(pos_unique)), ## gene BRCA2 is on chromosome 13
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
chr13 <- read.vcfR(file = paste0(dir,"ALL.chr13.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf"),
                   verbose = TRUE, nrows = 90000, convertNA = TRUE)
time2 = Sys.time() - time1
beepr::beep()

# we need to look for the right position to start importing the data:
pos1 <- chr13@fix[,2]
pos1 <- as.numeric(pos1)
range(pos1) # 18171423 21534858

## BRCA2 on Ensemble (GRCh38)
## Chromosome 13: 32315086-32400268
## > range(pos38)
## [1] 32315519 32399625

## after some trial and error to tune skip and nrows:
time1 = Sys.time()
chr13B <- read.vcfR(file = "ALL.chr13.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf",
                    verbose = TRUE, skip = 419060, nrows = 2400, convertNA = TRUE)
time2 = Sys.time() - time1
beepr::beep()

pos2 <- chr13B@fix[,2]
pos2 <- as.numeric(pos2)
range(pos2)  # 32315519 32404028

i_min = min(which(pos2 >= 32315519))
i_max = max(which(pos2 <= 32399625))

index_chr13B <- seq(i_min,i_max)
Gen1000_fix = chr13B@fix[index_chr13B,-8]    ## NOT ORDERED by pos38 (higher dimensionality)
Gen1000_gt = chr13B@gt[index_chr13B,-1]

Gen1000_fix = as.data.frame(Gen1000_fix)

pos_Gen1000 = as.numeric(Gen1000_fix[,2])
# index_ensembl = match(pos38,pos_Gen1000)
from_dating_to_ensembl = match(pos38,pos_Gen1000)
from_ensembl_to_dating = match(pos_Gen1000,pos38)

# these are the same, except for NA's
# head(pos_Gen1000[match(pos38,pos_Gen1000)])
# head(pos38)

sum(is.na(from_dating_to_ensembl)) # 111 dated species are not found in ensembl -- NOTE  I ran it a second time and got 109 na.
head(Gen1000_fix[from_dating_to_ensembl,])

from_dating_to_ensembl_noNA <- from_dating_to_ensembl[!is.na(from_dating_to_ensembl)]

varlist_save <- c(varlist_save, "Gen1000_fix","Gen1000_gt",
	"from_dating_to_ensembl","from_ensembl_to_dating", "from_dating_to_ensembl_noNA")

Gen1000_fix$AGE_MODE = NA
Gen1000_fix$AGE_MODE[from_dating_to_ensembl_noNA] = agemode[!is.na(from_dating_to_ensembl)]

Gen1000_fix$ID_dating <- NA
Gen1000_fix$ID_dating[from_dating_to_ensembl_noNA] = BRCA2_combined[,"VariantID"][!is.na(from_dating_to_ensembl)]

Gen1000_fix$POS_dating_GRCh37 <- NA
Gen1000_fix$POS_dating_GRCh37[from_dating_to_ensembl_noNA] = BRCA2_combined[,"Position"][!is.na(from_dating_to_ensembl)]

# we could potentially save a dataset, but it's NOT ordered by pos38 because the NA were removed.
BRCA2_combined_pos38NA <- BRCA2_combined[!is.na(from_dating_to_ensembl),]

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

biomart_data = read.csv(paste0(dir, "mart_export.txt"), header =TRUE)
colnames(biomart_data) <- c("variant_name","source","chrom","start","end")
dim(biomart_data)
mean(pos38 %in% biomart_data[,"start"]) # 1.00!!

from_pos38_to_mart = match(pos38,biomart_data[,"start"])
# all.equal(biomart_data[from_pos38_to_mart,"start"],pos38) #TRUE
# biomart_data[from_pos38_to_mart,"variant_name"]

# note: pos_unique has the same order of BRCA2_combined and
# pos38 is matched to pos_unique. So we can do
BRCA2_combined$VariantID_mart = biomart_data[from_pos38_to_mart,"variant_name"]
# let's save this also in Gen1000_fix
Gen1000_fix$ID_mart <- NA
Gen1000_fix$ID_mart[from_dating_to_ensembl_noNA] = BRCA2_combined[,"VariantID_mart"][!is.na(from_dating_to_ensembl)]

Gen1000_fix$POSstart_mart <- NA
Gen1000_fix$POSend_mart <- NA
Gen1000_fix$POSstart_mart[from_dating_to_ensembl_noNA] = biomart_data[from_pos38_to_mart,"start"][!is.na(from_dating_to_ensembl)]
Gen1000_fix$POSend_mart[from_dating_to_ensembl_noNA] = biomart_data[from_pos38_to_mart,"end"][!is.na(from_dating_to_ensembl)]

BRCA2_combined[1:10,c("VariantID","VariantID_mart")] # promising!
all.equal(BRCA2_combined[,"VariantID"],BRCA2_combined[,"VariantID_mart"]) # "25 string mismatches"
all.equal(BRCA2_combined[!is.na(from_dating_to_ensembl),"VariantID"],BRCA2_combined[!is.na(from_dating_to_ensembl),"VariantID_mart"])

equal_name = apply(BRCA2_combined[,c("VariantID","VariantID_mart")],MARGIN = 1, FUN = function(x) x[1]==x[2])
names(equal_name) = NULL

Gen1000_fix$Equal_name = NA
Gen1000_fix$Equal_name[from_dating_to_ensembl_noNA] = equal_name[!is.na(from_dating_to_ensembl)]

BRCA2_combined[!equal_name,c("VariantID","VariantID_mart")]
# the first one has a normal name: rs200268825 for dating, rs11571577 for mart. But it's not the only one with "normal names"
biomart_data[which(biomart_data[,"variant_name"]=="rs11571577"),] ## the first one has a end that's larger than the start (weird.)

biomart_data[biomart_data[,"variant_name"] %in% BRCA2_combined[!equal_name,"VariantID_mart"],]
## some of them, but not all, have weird start and end...

## this is difficult to understand
match( pos38[match(BRCA2_combined[!equal_name,c("Position")],pos_unique)],
       Gen1000_fix[,"POS"])
# this is exactly the same as before, but maybe easier?
from_dating_to_ensembl[match(BRCA2_combined[!equal_name,c("Position")],pos_unique)]
# ie. most of the different names also are missing in the ensemble data in terms of position

# the ones that are present in ensembl (Gen1000_fix)
# noequalname_inensembl = which(!is.na(from_dating_to_ensembl[match(BRCA2_combined[!equal_name,c("Position")],pos_unique)]))
## note which(!equal_name) is the same as match(BRCA2_combined[!equal_name,c("Position")],pos_unique)
# all.equal(match(BRCA2_combined[!equal_name,c("Position")],pos_unique), which(!equal_name) )
# noequalname_inensembl = which(!is.na(from_dating_to_ensembl[ !equal_name ])) ## wrong
noequalname_inensembl = which(  (!is.na(from_dating_to_ensembl))  &  (!equal_name)  )

BRCA2_combined[ noequalname_inensembl ,c("VariantID","VariantID_mart")] # different names but both rs...
# while the others are missing
biomart_data[from_pos38_to_mart[noequalname_inensembl],]

varlist_save <- c(varlist_save, "biomart_data","from_pos38_to_mart", "equal_name",
	"noequalname_inensembl")
save(list = varlist_save,file= "data/output_cleandata_BRCA.RData")


Gen1000_fix_keep = Gen1000_fix[which(Gen1000_fix$Equal_name == TRUE),]
sum(is.na(Gen1000_fix_keep[,-c(3,6)]))
Gen1000_fix_keep = Gen1000_fix_keep[,-c(3,6)]

Gen1000_gt_keep = Gen1000_gt[which(Gen1000_fix$Equal_name == TRUE),]
Gen1000_gt_keep01 = ifelse(Gen1000_gt_keep == "0|0",0,1)
tot_count = rowSums(Gen1000_gt_keep01)

table(tot_count)

varlist_save <- c(varlist_save, "Gen1000_fix_keep","Gen1000_gt_keep", "Gen1000_gt_keep01",
                  "tot_count")
save(list = varlist_save,file= "data/output_cleandata_BRCA.RData")


tmp = load("data/output_cleandata_BRCA.RData")

head(Gen1000_fix_keep$AGE_MODE)
Gen1000_fix_keep$rank <- NA

# this would give the larges (in terms of large N) rank to the oldest (large age)
# and rank 1 to the youngest (small age)
# match(Gen1000_fix_keep$AGE_MODE[1:10],sort(Gen1000_fix_keep$AGE_MODE[1:10], decreasing = FALSE))
# cbind(1:10,Gen1000_fix_keep$AGE_MODE[1:10])
Gen1000_fix_keep$rank <- match(Gen1000_fix_keep$AGE_MODE,
                               sort(Gen1000_fix_keep$AGE_MODE, decreasing = FALSE))

which.max(Gen1000_fix_keep$rank) # 42
max(Gen1000_fix_keep$rank)       # 1073

tot_count[42] # 2542 # one of the largest clusters!
              # there are other six species with this frequency (2542)
              # and the most frequent species has 2548 observations.
which(tot_count == max(tot_count))
Gen1000_fix_keep[c(42,545),] 
# the most frequent is another very old species (rank 1040 instead of 1073 and age 43K instead of 78K).

save(list = c("Gen1000_fix_keep","Gen1000_gt_keep", "Gen1000_gt_keep01"),
     file = "data/Gen1000_cleaned_BRCA.RData")

