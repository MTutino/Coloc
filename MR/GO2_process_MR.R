qtl.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/"
project.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/MR/"
setwd(project.path)

source(paste0(project.path, "scripts/functions_MR.R"))

################################################################################
#-------------------------- Process cis and trans MR --------------------------#
################################################################################
# Process cis MR
dt <- data.table::fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/MR/GO2_results_cis_MR.csv")
dt$analysis <- "cis"
dt[, c("tissue", "gene"):=data.table::tstrsplit(sub("eQTL_", "", exposure), "_(?=[^_]+$)", perl=TRUE)]

processed.dt <- function.processMR(dt)
processed.dt[, c("tissue", "gene"):=data.table::tstrsplit(sub("eQTL_", "", exposure), "_(?=[^_]+$)", perl=TRUE)]
processed.dt[p_value_plot==0, .N]
processed.dt[p_value_plot==0, data.table::uniqueN(gene)]
processed.dt[p_value_plot==0, .N, by=.(tissue, outcome)]
mr.genes <- processed.dt[p_value_plot==0, unique(gene)]
data.table::fwrite(processed.dt, paste0(project.path,"GO2_processed_results_cis.csv"))
data.table::fwrite(processed.dt[p_value_plot==0], paste0(project.path,"GO2_sig_processed_results_cis.csv"))

################################################################################
#------------------------ Compare with colocalization -------------------------#
################################################################################
coloc.tissue.trait <- unique(fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/coloc_per_gene_tissue_trait.csv"))
processed.dt <- fread("GO2_sig_processed_results_cis.csv")

# Add gene names
coloc.genes.dt <- fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/all_colocalizing_genes_with_names.csv")
mr.genes <- processed.dt[, unique(gene)]
mr.genes.dt <- merge(data.table(gene.id=mr.genes), coloc.genes.dt, by.x="gene.id", by.y="ensembl_gene_id", all.x=TRUE)
data.table::fwrite(mr.genes.dt, "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/GO2_coloc_mr_genes.txt")
