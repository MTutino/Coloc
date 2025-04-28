.libPaths(c("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/R", .libPaths()))

qtl.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/pipeline_RNA_seq_analysis/eQTL/"
go.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/odysseas/GO2/pQTL_new/GWAS_data/"  # b38
project.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/MR/"
setwd(project.path)

library(data.table)
library(ieugwasr)

################################################################################
#--------------------------------- FUNCTIONS ----------------------------------#
################################################################################
get.oa.rsid <- function(trait, project.path){
  load(paste0("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/GO2_GWAS_data/", trait, ".rda"))
  oa <- as.data.table(GWAS_associations)
  oa[, chr:=as.numeric(levels(chr))[chr]]
  oa <- oa[!(rsid %in% seq(1,22))]
  snps <- oa[!(rsid %in% c("", ".")), rsid]
  data.table::fwrite(data.table::data.table(SNP=snps), paste0(project.path, trait, ".rsid.txt"))
}

get.oa.alphaid <- function(trait, project.path){
  load(paste0("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/GO2_GWAS_data/", trait, ".rda"))
  oa <- as.data.table(GWAS_associations)
  rm(GWAS_associations)
  oa[, chr:=as.numeric(levels(chr))[chr]]
  oa[, alphaid:=paste(paste(chr, position, sep=":"), pmin(ea, nea), pmax(ea, nea), sep="_")]
  snps <- oa[, alphaid]
  data.table::fwrite(data.table::data.table(SNP=snps), paste0(project.path, trait, ".alphaid.txt"))
}

get.go.data <- function(trait, path){
  outcome_dat <- TwoSampleMR::read_outcome_data(filename = paste0(path, trait, "_subset.txt"),
                                                snps = snps,
                                                chr_col = "CHR",
                                                pos_col = "POS",
                                                sep = "\t",
                                                snp_col = "SNP",
                                                beta_col = "BETA",
                                                se_col = "SE",
                                                effect_allele_col = "EA",
                                                other_allele_col = "NEA",
                                                eaf_col = "EAF",
                                                pval_col = "P",
                                                samplesize_col = "N")
}

read_qtl_data <- function(path, tissue, cis=TRUE){
  if(cis){
    filename <- ifelse(tissue=="high_grade_cartilage", "hg", ifelse(tissue=="low_grade_cartilage", "lg", tissue))
    dt <- data.table::fread(paste0(path, tissue, "/tensorqtl/output_tensorqtl/nom/sig_pairs_df_", filename, ".txt"))
    dt[, c("chr", "pos", "REF", "ALT", "rsid"):=data.table::tstrsplit(variant_id, "_")]
    dt[, maf:=ifelse(af<0.5, af, 1-af)]
    dt <- dt[rsid!="."]
  }
  
  else{
    peer <- ifelse(tissue=="high_grade_cartilage", 30, ifelse(tissue=="low_grade_cartilage", 45, ifelse(tissue=="synovium", 45, 15)))

    dt <- data.table::fread(paste0(path, tissue, "/matrixEQTL/peer", peer, "/matrixeQTL_all_signif_pairs_transeVariant_eGene.txt.gz")) 
    dt[, c("chr", "pos", "REF", "ALT", "rsid"):=data.table::tstrsplit(SNP, "_")]
    dt[, chr:=as.integer(sub("chr", "", chr))]
    dt[, pos:=as.integer(pos)]
    dt[, slope_se:=abs(beta/qnorm(`p-value`/2))]
    colnames(dt)[match("gene", names(dt))] <- "phenotype_id"
    colnames(dt)[match("p-value", names(dt))] <- "pval_beta"
    colnames(dt)[match("beta", names(dt))] <- "slope"
    dt <- dt[rsid!="."]
    print(length(unique(dt$phenotype_id)))
    dt$ma_count <- ifelse(tissue=="high_grade_cartilage", 216, ifelse(tissue=="low_grade_cartilage", 263, ifelse(tissue=="synovium", 278, 94)))
  }
  
  return(dt)
}

extract_ivs_go <- function(ivs, trait.lst, project.path, out.file){
  for (trait in trait.lst){
    load(paste0("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/GO2_GWAS_data/", trait, ".rda"))
    oa <- as.data.table(GWAS_associations)
    rm(GWAS_associations)
    oa[, chr:=as.numeric(levels(chr))[chr]]
    oa[, alphaid:=paste(paste(chr, position, sep=":"), pmin(ea, nea), pmax(ea, nea), sep="_")]
    oa$cptid <- oa$rsid
    oa$rsid <- NULL
    oa[, ncases:=ifelse(trait=="KNEE", 172256, ifelse(trait=="ALLOA", 489952, 48161))]
    data.table::fwrite(oa[alphaid %in% ivs], paste0(project.path, trait, out.file))
  }
}

extract_ivs_qtl <- function(ivs.dt, tissue.lst, qtl.path, project.path, out.file, cis=TRUE){
  for (tissue in tissue.lst){
    qtl <- read_qtl_data(path=qtl.path, tissue, cis=cis)
    qtl <- qtl[phenotype_id %in% ivs.dt$gene.id & rsid %in% unique(c(ivs.dt$iv, ivs.dt$proxy))]
    data.table::fwrite(qtl, paste0(project.path, "eQTL_", tissue, out.file))
  }
}

local.clump <- function(data, ref.bfile.path, pval.thrs, pval.col="pval_beta", rsid.col="rsid"){
  # https://mrcieu.github.io/ieugwasr/articles/local_ld.html
  clumped.dt <- ieugwasr::ld_clump(dplyr::tibble(rsid=data[, get(rsid.col)], pval=data[, get(pval.col)]),
                                   plink_bin = genetics.binaRies::get_plink_binary(),
                                   bfile = ref.bfile.path,
                                   clump_p = pval.thrs,
                                   clump_kb = 10000,
                                   clump_r2 = 0.001)
  return(clumped.dt$rsid)
}

proxy_lookup <- function(snp, api.token, population, gwas.data){
  chosen.proxy <- NA
  all.potential.proxies <- LDlinkR::LDproxy(snp = snp, 
                                            pop = population, 
                                            r2d = "r2", 
                                            token = api.token,
                                            genome_build = "grch38_high_coverage")
  if(unique(grepl("error", all.potential.proxies))) return("no proxy found")
  
  all.potential.proxies <- data.table::as.data.table(all.potential.proxies)
  all.potential.proxies <- all.potential.proxies[R2>0.7 & RS_Number!=".", RS_Number]
  
  #Check if any of the proxies is included in the oa datad
  if(nrow(subset(gwas.data[rsid %in% all.potential.proxies]))>0){
    for (potential.proxy in all.potential.proxies){
      # Check if potential.proxy is in T2D data
      if(nrow(gwas.data[rsid==potential.proxy])!=0) break
    }
    chosen.proxy <- potential.proxy
  }
  
  return(chosen.proxy)
}

get_ivs <- function(exposure, ref.bfile.path, pval.thrs, beta.col="slope", se.col="slope_se", rsid.col="rsid"){
  #Clumping
  clumped.ivs <- local.clump(data=exposure, ref.bfile.path, pval.thrs)
  
  #F-stat
  exposure[, snp_fstat:=(get(beta.col))^2/(get(se.col))^2]
  
  return(exposure[get(rsid.col) %in% clumped.ivs & snp_fstat>=10, get(rsid.col)])
}

preprocess_data <- function(tissues.lst, outcome.lst, project.path, out.file, ref.bfile.path, api.token, coloc.results, cis=TRUE){
  # Read eQTL data for hub genes
  ivs.dt <- data.table::data.table()
  
  for(gwas in outcome.lst){
    # Read OA data
    oa.snps <- data.table::fread(paste0(project.path, gwas, ".alphaid.txt")) 
    for (t in tissues.lst){
      dt <- read_qtl_data(path=qtl.path, tissue=t, cis=cis)
      dt[, alphaid:=paste(paste(sub("chr", "", chr), pos, sep=":"), pmin(REF, ALT), pmax(REF, ALT), sep="_")]
      cur.coloc.res <- coloc.results[GWAS_ID==gwas & tissue==t]
      cur.oa.snps <- merge(oa.snps, dt[,.(SNP=alphaid, rsid)], by="SNP")
      
      for (gene in unique(cur.coloc.res$gene)){  
        # Get IVs
        if(nrow(dt[phenotype_id==gene])>1){
          clumped.ivs <- tryCatch(
            get_ivs(exposure=dt[phenotype_id==gene], ref.bfile.path, pval.thrs=1),
            error=function(e) e
          )
          
          if(inherits(clumped.ivs, "error")) {
            print("No SNPs after clumping failed.")
            next 
          }
        }
        else clumped.ivs <- dt[phenotype_id==gene, rsid]
        
        cur.ivs.dt <- dt[rsid %in% clumped.ivs, .(alphaid, rsid, tissue=t, gene.id=gene, iv=clumped.ivs, proxy.needed=FALSE)]
        
        # Check if we need proxies for any IV
        need.proxy <- cur.ivs.dt[(!alphaid %in% oa.snps$SNP), rsid]
        cur.ivs.dt[rsid %in% need.proxy, proxy.needed:=TRUE]
        
        #Find proxies
        if(length(need.proxy)>0){
          # print(paste0("Gene: ", gene))
          for(i in 1:length(need.proxy)){
            chosen.proxy <- proxy_lookup(snp=need.proxy[i], api.token, population="EUR", gwas.data=cur.oa.snps)
            # create data table with columns: iv, OA trait, proxy
            cur.ivs.dt[rsid==need.proxy[i], `:=`(oa.trait=gwas, proxy=chosen.proxy)]
          }
        }
        
        ivs.dt <- rbind(ivs.dt, cur.ivs.dt, fill=T)  
      }
    }
    data.table::fwrite(ivs.dt, out.file)
  }
}

################################################################################
#------------------------------- MAIN cis MR ----------------------------------#
################################################################################
# Extract rsids for OA phenotypes of interest
get.oa.alphaid("KNEE", project.path="/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/MR/")
get.oa.alphaid("ALLOA", project.path="/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/MR/")
get.oa.alphaid("TKR", project.path="/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/MR/")

coloc.results <- fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/coloc_per_gene_tissue_trait.csv")

# Get IVs (including proxies if needed)
for(out in c("ALLOA", "KNEE", "TKR")){  # 
  preprocess_data(tissues.lst = c("high_grade_cartilage", "low_grade_cartilage", "synovium", "fat_pad"), 
                  outcome.lst = out, 
                  project.path, 
                  out.file = paste0(project.path, "data/", out, "_ivs.txt"),
                  ref.bfile.path = "/lustre/groups/itg/shared/referenceData/1kG/EAS-SAS-AMR-EUR-AFR_1kg_v3/EUR", 
                  api.token = "d1fe6afcb763",
                  coloc.results = coloc.results,
                  cis=TRUE)
}

ivs.dt <- rbind(data.table::fread(paste0(project.path, "data/KNEE_ivs.txt")),
                data.table::fread(paste0(project.path, "data/TKR_ivs.txt")),
                data.table::fread(paste0(project.path, "data/ALLOA_ivs.txt")))

# Extract IVs from GO1 data
extract_ivs_go(ivs=unique(c(ivs.dt$alphaid, ivs.dt$proxy)), 
               trait.lst=c("KNEE", "TKR", "ALLOA"), 
               project.path, 
               out.file="_subset.txt")

# Extract IVs from FunGen data
extract_ivs_qtl(ivs.dt, 
                tissue.lst = c("high_grade_cartilage", "low_grade_cartilage", "synovium", "fat_pad"), 
                qtl.path, 
                project.path, 
                out.fil="_subset.txt")

# Run MR
source(paste0(project.path, "scripts/functions_MR.R"))
wrap.MR(tissue.lst = c("high_grade_cartilage", "low_grade_cartilage", "synovium", "fat_pad"), 
        trait.lst = c("KNEE", "TKR", "ALLOA"),
        in.file = "_subset.txt",
        coloc.results = coloc.results,
        project.path = project.path, 
        output.file = paste0(project.path,"GO2_results_cis_MR.csv"))
