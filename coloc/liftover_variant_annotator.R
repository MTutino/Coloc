project.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/"
data.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/"
setwd(project.path)

################################################################################
#----------------------- Liftover variant_ann ----------------------------------
################################################################################
variant_ann <- fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/FunctionsAndData/variants.tsv", data.table=FALSE)
variant_ann$variant_id_chrpos <- paste(variant_ann$chr, variant_ann$pos, toupper(variant_ann$alt), sep="_")

# Create .bed 
variant_ann <- as.data.table(variant_ann)
variant_ann$start <- variant_ann$pos-1 
fwrite(variant_ann[, .(chr=as.integer(chr), start=as.integer(start), pos=as.integer(pos), minor_allele, rsid)], "variant_ann_hg19.bed", col.names = F, row.names = F, quote = F, sep = "\t")

system("CrossMap.py bed /lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/FunctionsAndData/hg19ToHg38.over.chain /lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/variant_ann_hg19.bed /lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/variant_ann_hg38.bed 2>&1")

# Import the lifted positions and keep only autosomes
liftedover = data.table::fread("variant_ann_hg38.bed")
colnames(liftedover) <- c("chr.b38", "start", "pos.b38", "ea", "rsid")
dim(liftedover)

# Create the final file
liftedover[, `:=` (start=NULL, ea=NULL)]
liftedover <- unique(liftedover)
res <- merge(variant_ann, liftedover, by="rsid")
variant_ann <- res
fwrite(res, "variant_ann_hg38.txt.gz")
rm(liftedover)
rm(res)
variant_ann$variant_id_chrpos <- paste(variant_ann$chr.b38, variant_ann$pos.b38, toupper(variant_ann$alt), sep="_")
fwrite(variant_ann, "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/variant_ann_hg38.txt.gz")