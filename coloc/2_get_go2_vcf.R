.libPaths(c("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/R", .libPaths()))
project.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/"
data.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/"
setwd(project.path)

source("scripts/Coloc_helper_functions.R")

library(data.table)
library(VariantAnnotation)
library(gwasvcf)
library(magrittr)
library("GenomicRanges")
library(rtracklayer)
library(coloc)
library(dplyr)

variant_ann <- fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/variant_ann_hg38.txt.gz")

################################################################################
#----------------------- Get GO in VCF format ----------------------------------
################################################################################
for(trait in c("KNEE", "TKR", "ALLOA")){
  ncases <- ifelse(trait=="KNEE", 172256, ifelse(trait=="TKR", 48161, 489952))
  ncontrols <- ifelse(trait=="KNEE", 1144244, ifelse(trait=="TKR", 958463, 1472094))
  GWAS_n <- c(ncases, ncontrols)
  GWASfile=paste0("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/GO2_", trait, "_sumstats_hg38.txt.gz")
  make_vcf(GWASfile=GWASfile, 
           chrom="chr.b38", 
           pos="pos.b38", 
           nea="NEA", 
           ea="EA", 
           snp="SNP", 
           ea_af="EAF", 
           effect="BETA", 
           se="SE", 
           pval="P",
           hg="hg38", # hg build of summary stats
           lift_down=FALSE, # if lift_down=TRUE then hg38 is converted to hg19. If FALSE, hg19 is converted to hg38. Ignored if WantToLiftOver is set ot FALSE 
           WantToLiftOver=FALSE, # Whether or not you want to lift over the coordinates
           GWAS_n=GWAS_n, # a vector of one or two elements. If quantitative trait, total sample size, if case control, number of cases and controls
           variant_ann=variant_ann, # reference file to map the missing values
           output=NULL)
}