project.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/"
setwd(project.path)

library(data.table)

res <- data.table()
overlap.snps <- NULL
for(trait in c("KneeOA", "AllOA", "TKR")){
  # Read MR IVs
  ivs.dt <- fread(paste0(project.path, "MR/data/", trait, "_ivs.txt"))
  
  for(tissue in c("high_grade_cartilage", "low_grade_cartilage", "synovium", "fat_pad")){
    # Read coloc credset
    coloc <- fread(file.path("coloc", paste0(trait, "_coloc_rda_files"), paste0("colocABF_df_", tissue, "_results.txt")))
    sig.coloc <- coloc[PP4 >= 0.8 | (PP4 > 0.6 & PP4 < 0.8 & PP4/PP3 > 2)]
    for(gene in unique(sig.coloc$cpg)){
      credset <- unlist(strsplit(sig.coloc[cpg==gene, credible_set], ","))
      if(length(credset[credset %in% ivs.dt[tissue==tissue & gene.id==gene, ifelse(!is.na(proxy), proxy, iv)]])!=0){
        res <- rbind(res, data.table(tissue=tissue,
                                     trait=trait,
                                     gene=gene,
                                     overlap=paste(credset[credset %in% ivs.dt[tissue==tissue & gene.id==gene, ifelse(!is.na(proxy), proxy, iv)]], collapse=",")))
        overlap.snps <- c(overlap.snps, credset[credset %in% ivs.dt[tissue==tissue & gene.id==gene, ifelse(!is.na(proxy), proxy, iv)]])
      }
    }
  }
}

fwrite(res, "mr_coloc_overlap_per_trait_tissue.tsv")
write(unique(overlap.snps), "mr_coloc_overlap_snps.tsv")

