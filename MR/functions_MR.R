library(dplyr)

add_odds_ratios <- function(dt, beta="beta", se="standard_error"){
  dt$or <- exp(dt[, get(beta)])
  dt$or_lci95 <- exp(dt[, get(beta)] - 1.96 * dt[, get(se)])
  dt$or_uci95 <- exp(dt[, get(beta)] + 1.96 * dt[, get(se)])

  return(dt)
}

run.SensitivityAnalysis <- function(data, res, plt=TRUE, het=TRUE, steiger = TRUE){
  fstat_overall <- mean((data$beta.exposure)^2/(data$se.exposure)^2)
  
  res[, Fstat:=fstat_overall]
  if(het){
    # check heterogeneity with Cochran's Q-statistic
    het <- data.table::as.data.table(TwoSampleMR::mr_heterogeneity(data)) 
    res <- merge(res, het[, .(exposure, outcome, method, Q, Q_df, Q_pval)],
                 by=c("exposure", "outcome", "method"), all.x=TRUE)
  }
  if(plt){
    # Assumption 2: check pleiotropy with MR-Egger intercept
    plt <- data.table::as.data.table(TwoSampleMR::mr_pleiotropy_test(data))  #MR Egger intercept for directional pleiotropy
    res[method=="MR Egger", `:=` (egger_intercept=plt$egger_intercept, egger_intercept_se=plt$se, egger_intercept_pval=plt$pval)]
  }
  
  if(steiger){
    #Run the Steiger filtering
    data.steiger <- TwoSampleMR::steiger_filtering(data)
    if(nrow(data.steiger)==0) return(res)
    data.steiger <- subset(data.steiger, steiger_dir & steiger_pval<0.05)
    res.steiger <- data.table::as.data.table(TwoSampleMR::mr(data.steiger, method = "mr_ivw"))
    # res.steiger <- data.table::as.data.table(TwoSampleMR::generate_odds_ratios(res.steiger))
    res.steiger[, `:=`(lo_ci=NULL, up_ci=NULL, id.exposure=NULL, id.outcome=NULL)]
    res.steiger$method <- "Steiger Inverse variance weighted"
    res <- rbind(res, res.steiger, fill = TRUE)
  }
  return(res)
}

run.TwoSampleMR <- function(data){
  if(nrow(data)==0)  {
    print("No SNP in common between exposure and outcome :(")
    tmp.dt <- data.table::data.table()
    return(tmp.dt)
  }
  else if (nrow(data)==1) {
    fstat <- (data$beta.exposure)^2/(data$se.exposure)^2
    print("Only one SNP in common between exposure and outcome, running Wald ratio method")
    res <- TwoSampleMR::mr(data, method_list="mr_wald_ratio")
    res <- data.table::as.data.table(TwoSampleMR::generate_odds_ratios(res))
    res[, `:=`(lo_ci=NULL, up_ci=NULL, id.exposure=NULL, id.outcome=NULL)]
    
    tmp.dt <- data.table::as.data.table(res) %>% .[, `:=`(Fstat=fstat,Q=as.numeric(NA),Q_df=as.numeric(NA),Q_pval=as.numeric(NA),plt.egger_intercept=as.numeric(NA),plt.pval=as.numeric(NA),plt.se=as.numeric(NA))]
    return(tmp.dt)
  }
  else {
    print("Multiple SNPs in common between exposure and outcome, running MR")
    res <- TwoSampleMR::mr(data, method_list=c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
    res <- data.table::as.data.table(TwoSampleMR::generate_odds_ratios(res))
    res[, `:=`(lo_ci=NULL, up_ci=NULL, id.exposure=NULL, id.outcome=NULL)]
    
    # Run sensitivity analysis
    res <- run.SensitivityAnalysis(data, res)
    
    return(res)
  }
}

get.go.data <- function(dt){
  outcome_dat <- TwoSampleMR::format_data(type="outcome",
                                          dat = dt,
                                          chr_col = "chr",
                                          pos_col = "position",
                                          snp_col = "rsid",
                                          beta_col = "beta",
                                          se_col = "se",
                                          effect_allele_col = "ea",
                                          other_allele_col = "nea",
                                          eaf_col = "eaf",
                                          pval_col = "p",
                                          ncase_col = "ncases",
                                          samplesize_col = "n")
}

get.qtl.data <- function(exp.var){
  data <- TwoSampleMR::format_data(type = "exposure",
                                   dat = as.data.frame(exp.var),
                                   snp_col = "rsid",
                                   beta_col = "slope",
                                   se_col = "slope_se",
                                   effect_allele_col = "ALT",
                                   other_allele_col = "REF",
                                   eaf_col = "af",
                                   pval_col = "pval_nominal",
                                   chr_col = "chr",
                                   pos_col = "pos",
                                   samplesize_col = "ma_count",
                                   gene_col = "phenotype_id")
  # info_col = "eQTL")
  return(data)
}

wrap.MR <- function(tissue.lst, trait.lst, coloc.results, project.path, in.file, output.file){
  dt <- data.table::data.table()
  for (cur.tissue in tissue.lst){
    #----------- Read exposure -----------
    exposure.all <- data.table::fread(paste0(project.path, "eQTL_", cur.tissue, in.file))
    exposure.all[, alphaid:=paste(paste(sub("chr", "", chr), pos, sep=":"), pmin(ALT, REF), pmax(ALT, REF), sep="_")]
    
    for (cur.trait in trait.lst){
      #----------- Get outcome -----------
      gwas.trait <- unique(data.table::fread(paste0(project.path, cur.trait, in.file)))
      gwas.trait <- merge(gwas.trait, unique(exposure.all[, .(alphaid, rsid)]), by="alphaid")
      outcome <- get.go.data(as.data.frame(gwas.trait))
      outcome$outcome <- cur.trait
      
      cur.coloc.res <- coloc.results[GWAS_ID==cur.trait & tissue==cur.tissue]
      for (cur.gene in unique(cur.coloc.res$gene)){
        exposure <- get.qtl.data(as.data.frame(exposure.all[phenotype_id==cur.gene]))
        exposure$exposure <- paste("eQTL", cur.tissue, cur.gene, sep="_")
        
        #----------- Harmonize data -----------
        data <- tryCatch(
          TwoSampleMR::harmonise_data(exposure_dat=exposure, outcome_dat=outcome),
          error=function(e) e
        )
        
        if(inherits(data, "error")) {
          print("Harmonization failed.")
          next 
        }
        
        if(nrow(data)==0) {
          print("No SNP left after data harmonization, not able to run MR :(")
          next
        }
        
        if(nrow(subset(data, mr_keep==TRUE))==0) {
          print("No SNP left after data harmonization, not able to run MR :(")
          next
        }
        
        data$id.exposure <- data$exposure
        data$id.outcome <- cur.trait    
        data <- subset(data, mr_keep==TRUE)
        
        #----------- Run TwoSampleMR with sensitivity analyses -----------
        res.mr <- run.TwoSampleMR(data)
        dt <- rbind(dt, res.mr, fill=TRUE)
      }
    }
  }
  data.table::fwrite(dt, output.file)
  dt[method %in% c("Inverse variance weighted", "Wald ratio"), p.ivw.adj.fdr:=p.adjust(pval, method = "fdr")] %>% data.table::setnames(., old=c("b", "se", "pval", "or", "or_lci95", "or_uci95", "egger_intercept_se", "egger_intercept_pval", "p.ivw.adj.fdr"),
                                                                                                                           new=c("beta", "standard_error", "p_value", "odds_ratio", "ci_lower", "ci_upper", "egger_intercept_standard_error", "egger_intercept_pvalue", "p_value_ivw_fdr"))
  dt <- dt[,which(unlist(lapply(dt, function(x)!all(is.na(x))))),with=F]
  dt[, gene:=sub(".*_", "", exposure)]
  data.table::fwrite(dt, output.file)
}

function.processMR <- function(res.mr, weighted.mode.sensitivity = T){
  #Agreed not to look at weighted mode
  if(!weighted.mode.sensitivity){
    res.mr <- subset(res.mr, method != "Weighted mode")
  }
  
  #Subset the results to keep the IVW p-value and effect size if at least 2 snps
  res.mr.IVW <- subset(res.mr, method == "Inverse variance weighted" & nsnp >= 2)
  #Only keep Wald Ratio if only one or two SNPs
  res.mr.wald <- subset(res.mr, method == "Wald ratio")
  
  res.mr.all <- data.table::data.table()
  
  ###########################Sensitivity analyses
  #Flag if the MR Egger intercept is significant
  res.mr.MREgger <- subset(res.mr, method == "MR Egger")
  res.mr.MREgger$pval_fdr <- p.adjust(res.mr.MREgger$egger_intercept_pvalue, method = "fdr")
  MREgger.sig <- subset(res.mr.MREgger, pval_fdr < 0.05)$exposure
  # MREgger.sig <- subset(res.mr.MREgger, egger_intercept_pvalue < 0.05)$exposure
  #Look at whether heterogeneity is above threshold
  res.mr.IVW$Q_pval_fdr <- p.adjust(res.mr.IVW$Q_pval, method = "fdr")
  Het.sig <- subset(res.mr.IVW, Q_pval_fdr < 0.05)$exposure
  # Het.sig <- subset(res.mr.IVW, Q_pval < 0.05)$exposure
  #Take the beta of each sensitivity method, with MR PRESSO depending on distortion test
  #Check if MR PRESSO has been run at least one
  
  if(nrow(res.mr)==0) next
  if(any(grepl(res.mr$method, pattern = "MR_PRESSO"))){
    #Modify the distortion p-value field to have it in numeric
    subset(res.mr, method == "MR_PRESSO_IVW_Outliers")[, Distortion_pval:=ifelse(Distortion_pval == "<0.001", 1e-4, as.numeric(Distortion_pval))]
    res.mr[, Distortion_pval:=as.numeric(Distortion_pval)]
    beta.sensitivity <- data.frame(MRPRESSO = ifelse(subset(res.mr, method == "MR_PRESSO_IVW_Outliers")$Distortion_pval<0.05 & !is.na(subset(res.mr, method == "MR_PRESSO_IVW_Outliers")$Distortion_pval), subset(tmp.res, method == "MR_PRESSO_IVW_Outliers")$beta, subset(tmp.res, method == "MR_PRESSO_IVW_Raw")$beta),
                                   MREgger = ifelse(length(subset(res.mr, method == "MR Egger")$beta)!=0, subset(res.mr, method == "MR Egger")$beta, NA),
                                   WeightedMedian = ifelse(length(subset(res.mr, method == "Weighted median")$beta)!=0, subset(res.mr, method == "Weighted median")$beta, NA),
                                   Steiger = ifelse(length(subset(res.mr, method == "Steiger Inverse variance weighted")$beta)!=0, subset(res.mr, method == "Steiger Inverse variance weighted")$beta, NA))
  }else{
    beta.sensitivity <- data.frame(MREgger = ifelse(length(subset(res.mr, method == "MR Egger")$beta)!=0, subset(res.mr, method == "MR Egger")$beta, NA),
                                   WeightedMedian = ifelse(length(subset(res.mr, method == "Weighted median")$beta)!=0, subset(res.mr, method == "Weighted median")$beta, NA),
                                   Steiger = ifelse(length(subset(res.mr, method == "Steiger Inverse variance weighted")$beta)!=0, subset(res.mr, method == "Steiger Inverse variance weighted")$beta, NA))
  }
  #Flag the mQTL that don't have the all tests with same direction of effect
  Prop.SameDir <- sapply(unique(res.mr.IVW$exposure), function(z) mean(sign(beta.sensitivity[z,]) == sign(subset(res.mr.IVW, exposure == z)$beta), na.rm = T))
  Not.SameDir <- unique(res.mr.IVW$exposure)[which(Prop.SameDir<1)]
  # res.mr.IVW[outcome==out, DiffDirection:=ifelse(Prop.SameDir!=0, TRUE, FALSE)]

  #Combine the two files
  res.mr.out <- rbind(res.mr.IVW, res.mr.wald, fill=TRUE)
  
  #Flag the results that have significant het, pleio
  #Keep heteorgeneity from Radial filtering if considered, and don't consider MREgger
  res.mr.out$FlagSensitivity <- ifelse(res.mr.out$exposure %in% unique(c(MREgger.sig, Het.sig)), T, F)
  res.mr.out$FlagPleiotropy <- ifelse(res.mr.out$exposure %in% unique(c(MREgger.sig)), T, F)
  res.mr.out$FlagHeterogeneity <- ifelse(res.mr.out$exposure %in% unique(c(Het.sig)), T, F)
  
  #Flag the results that don't have same direction 
  res.mr.out$DiffDirection <- ifelse(res.mr.out$exposure %in% Not.SameDir, T, F)

  #Compute the final ajusted p-value
  res.mr.out$p_value_fdr <- p.adjust(res.mr.out$p_value, method = "fdr")
  res.mr.out[, p_value_plot:=ifelse(p_value_fdr<0.05 & DiffDirection==FALSE & FlagSensitivity==FALSE, 0, 2)]
  
  return(res.mr.out)                            
}

#### Testing
# project.path <- "/lustre/groups/itg/shared/notebooks/georgia.katsoula/RNA_seq_300/Final_RNA_seq_300_paper/Molecular_environment_OA/results/MR/"
# # project.path <- "C:/Users/ana.arruda/OneDrive - Helmholtz Zentrum München/Projects/MR_for_Georgia/"
# setwd(project.path)
# trait.lst=c("KneeOA", "TKR", "AllOA", "TJR")
# tissue.lst = c("HighGradeCartilage", "LowGradeCartilage", "Synovium")
# tissue=tissue.lst[1]
# output.file=paste0(project.path,"results/results.csv")
# library(data.table)
