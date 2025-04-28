#prepare annotations for HiC
project.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/"
setwd(project.path)

library(data.table)
library(openxlsx)
library(dplyr)

chondro_hic <- read.xlsx('/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/torus/ard-2023-224945-inline-supplementary-material-7.xlsx', startRow = 2)
chondro_hic$Chromosome <- paste0('chr',chondro_hic$Chromosome)
chondro_hic_tad <- unique(chondro_hic[,1:3])
chondro_hic_ctcf <- chondro_hic[,c(1,4:5)]
#chondro_hic_tad$Feature <- paste('Chondrocyte_TAD')
#chondro_hic_ctcf$Feature <- paste('Chondrocyte_CTCF_binding_site')

#only these are annotated to respective genes
chondro_prom_enh <- read.xlsx('/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/torus/ard-2023-224945-inline-supplementary-material-2.xlsx',startRow = 2)
chondro_prom_enh$Chromosome <- paste0('chr',chondro_prom_enh$Chromosome)
chondro_prom_enh <- unique(chondro_prom_enh[,c(1:3,5)])
#chondro_prom_enh$Feature <- paste('Chondrocyte_active_promoter')

chondro_enhancer <- read.csv('/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/torus/All_Active_Enhancers_in_PromEnhLoopsList_corrected_source_single_ATAC.csv')
chondro_enhancer <- chondro_enhancer %>% dplyr::rename('Chromosome' = chr)
chondro_enhancer <- unique(chondro_enhancer[,1:3])
#chondro_enhancer$Feature <- paste('Chondrocyte_active_enhancer')

rownames(chondro_hic_tad) <- NULL
rownames(chondro_hic_ctcf) <- NULL
rownames(chondro_prom_enh) <- NULL
rownames(chondro_enhancer) <- NULL

chondro_hic_ctcf <- chondro_hic_ctcf %>% dplyr::filter(Chromosome %in% paste0('chr',1:22))
chondro_prom_enh <- chondro_prom_enh %>% dplyr::filter(Chromosome %in% paste0('chr',1:22))
chondro_enhancer <- chondro_enhancer %>% dplyr::filter(Chromosome %in% paste0('chr',1:22))


colnames(chondro_hic_ctcf) <- c('Chromosome','Start','End','Feature')
colnames(chondro_enhancer) <- c('Chromosome','Start','End','Feature')
colnames(chondro_prom_enh) <- c('Chromosome','Start','End','Feature')


chondro_hic_list <- list(chondro_hic_ctcf, chondro_enhancer, chondro_prom_enh)
names(chondro_hic_list) <- c('Chondrocyte_CTCF_binding_site','Chondrocyte_active_enhancer','Chondrocyte_active_promoter')
annotation_list_all_cart <- names(chondro_hic_list)

#save as bed
annotations_path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/HiCannotation"

#write the bed files
lapply(1:length(chondro_hic_list),function(i){
  readr::write_tsv(chondro_hic_list[[i]], 
                   file = file.path(annotations_path, paste0(annotation_list_all_cart[i],'.bed')),
                   col_names=FALSE)
})

coloc.dt <- data.table()
for(trait in c("KneeOA", "AllOA", "TKR")){
  for(tissue in c("high_grade_cartilage", "low_grade_cartilage", "synovium", "fat_pad")){
    # Read coloc results
    coloc <- fread(file.path("coloc", paste0(trait, "_coloc_rda_files"), paste0("colocABF_df_", tissue, "_results.txt")))
    # Subset to significant results
    sig.coloc <- coloc[PP4 >= 0.8 | (PP4 > 0.6 & PP4 < 0.8 & PP4/PP3 > 2)]
    for(gene in unique(sig.coloc$cpg)){
      # Extract the credible
      credset <- unlist(strsplit(sig.coloc[cpg==gene, credible_set], ","))
      #cur.chr <- unique(sig.coloc[cpg==gene, coloc_lead_snp_chr])
      coloc.dt <- rbind(coloc.dt, data.table(rsid=credset, gene=gene, trait=trait, tissue=tissue))
    }
  }
}

# Add POS of rsids
ensembl <- biomaRt::useEnsembl("snp",dataset = "hsapiens_snp")
SNPs <- biomaRt::getBM(attributes=c("refsnp_id",
                                    "chr_name",
                                    "chrom_start",
                                    "chrom_end"),
                       filters="snp_filter", values=unique(coloc.dt$rsid), mart = ensembl, uniqueRows=TRUE)
SNPs <- as.data.table(SNPs) %>% .[chr_name %in% seq(1,22)]
coloc.dt <- merge(coloc.dt, SNPs[, .(refsnp_id, chr=chr_name, pos=chrom_start)], by.x="rsid", by.y="refsnp_id", all.x=TRUE)
fwrite(coloc.dt, "colocalization_credset_snps.csv")

library(GenomicRanges)

snpRanges <- makeGRangesFromDataFrame(coloc.dt[!is.na(pos)], seqnames.field='chr', start.field='pos', end.field='pos', keep.extra.columns = TRUE)
snpRanges <- plyranges::mutate(snpRanges, snp=rsid)

annotations <- list.files(annotations_path, full.names = TRUE)

# head(chondro_hic_list[[3]])
# names(chondro_hic_list[[3]])
# f=annotations[2]
for(f in annotations){
  cat(sprintf('Annotating SNPs with %s ...\n', basename(f)))
  annot.gr <- rtracklayer::import(f, format='bed')
  # snpRangesIn <- IRanges::subsetByOverlaps(snpRanges, annot.gr)
  snpRangesIn <- IRanges::mergeByOverlaps(snpRanges, annot.gr)
  snpRangesIn[snpRangesIn$gene==snpRangesIn$name]
  res.dt <- as.data.table(snpRangesIn)
  res.dt <- res.dt[, .(snp, chr=snpRanges.seqnames, pos=snpRanges.start, trait, tissue, coloc.gene=gene, hic.gene=name, hic.start=annot.gr.start, hic.end=annot.gr.end)]
  res.dt[coloc.gene==hic.gene, .(snp, chr, pos, trait, tissue, coloc.gene)]
  snpsIn <- unique(snpRangesIn$snp)
  unique(res.dt[coloc.gene==hic.gene, .(snp, chr, pos, trait, tissue, coloc.gene)], by=c("snp", "coloc.gene"))
  
  fwrite(res.dt[coloc.gene==hic.gene, .(snp, chr, pos, trait, tissue, coloc.gene)], paste0("HiCannotation/result_", sub(".bed", "", basename(f)), ".txt"))
  fwrite(unique(res.dt[coloc.gene==hic.gene, .(snp, chr, pos, trait, tissue, coloc.gene)], by=c("snp", "coloc.gene")), paste0("HiCannotation/unique_result_", sub(".bed", "", basename(f)), ".txt"))
}

library(data.table)
library(dplyr)

project.path <- "C:/Users/AnaLuizaArruda/OneDrive - Helmholtz Zentrum München/Projects/molQTL300"
setwd(project.path)

credset.snps <- fread("coloc/colocalization_credset_snps.csv")
enhancer.dt <- fread("ActEnh_ActPromoter_Loops_Source_plus_Gene_info.csv")
enhancer.dt <- enhancer.dt[chr!="X"]
enhancer.dt$gene <- enhancer.dt$ENSEMBL.gene.ID
enhancer.dt$chr <- as.integer(enhancer.dt$chr)

matched_pairs <- credset.snps %>%
  inner_join(enhancer.dt, by = c("chr", "gene"), relationship = "many-to-many") %>%
  filter(pos >= start_enh & pos <= end_enh)

fwrite(unique(matched_pairs[source=="fanc"]), "overlap_credset_enhancer.csv")
# View the result
