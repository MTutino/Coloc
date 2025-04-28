library(data.table)

setwd("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/")

local.clump <- function(data, ref.bfile.path, pval.thrs, pval.col="pval_beta", rsid.col="rsid"){
  # https://mrcieu.github.io/ieugwasr/articles/local_ld.html
  clumped.dt <- ieugwasr::ld_clump(dplyr::tibble(rsid=data[, get(rsid.col)], pval=data[, get(pval.col)]),
                                   plink_bin = genetics.binaRies::get_plink_binary(),
                                   bfile = ref.bfile.path,
                                   clump_p = pval.thrs,
                                   clump_kb = 10000,
                                   clump_r2 = 0.001)
  # return(data[data$rsid %in% clumped.dt$rsid, ])
  return(clumped.dt$rsid)
}

res.dt <- data.table()
for(t in c("high_grade_cartilage", "low_grade_cartilage", "synovium", "fat_pad")){
  peer <- ifelse(t=="high_grade_cartilage", 30, ifelse(t=="fat_pad", 15, 45))
  dt <- fread(paste0("pipeline_RNA_seq_analysis/eQTL/", t, "/matrixEQTL/peer", peer, "/matrixeQTL_all_signif_pairs_transeVariant_eGene.txt.gz"))
  dt[, `:=`(rsid=sub(".*_", "", SNP), tissue=t)]
  
  for(g in unique(dt$gene)){
    if(nrow(dt)>1){
      clumped.snp <- tryCatch({local.clump(dt[gene==g], 
                                           ref.bfile.path="/lustre/groups/itg/shared/referenceData/1kG/EAS-SAS-AMR-EUR-AFR_1kg_v3/EUR",
                                           pval.thrs=1, pval.col="p-value", rsid.col="rsid")},
                              error=function(e) e)
      if(inherits(clumped.snp, "error")) res.dt <- rbind(res.dt, data.table(gene=g,
                                                                            tissue=t,
                                                                            rsid="error"), fill=T)
      else res.dt <- rbind(res.dt, dt[gene==g & rsid %in% clumped.snp], fill=T)
    }
    else res.dt <- rbind(res.dt, dt[gene==g], fill=T)
  }
}

res.dt[rsid=="."]  # are not clumped (only one IV)
res.dt[is.na(rsid)]

fwrite(res.dt, "pipeline_RNA_seq_analysis/eQTL/trans_eqtl_clumped_all_tissues.tsv")

