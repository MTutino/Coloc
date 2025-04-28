perform_coloc = function(overlap_df, in_m_qtl, out_path, in_gwas, gwas_cc_ratio, tissue, gwas_n = NULL, GTEX_APP = TRUE, GWAS_type="quant", GWAS_ID=GWAS_ID, QTL_n=QTL_n, is.eQTL=FALSE, cores=3){
        ##############################
        number_tests = nrow(overlap_df)
        counter = 1
        ##############################
        coloc_df = do.call(rbind, mclapply(1:nrow(overlap_df), FUN = function(coloc_nr){
            ##############################
            print(paste0("Running coloc test: ", coloc_nr, " of ", number_tests))    
            ##############################
            cpg_coloc = overlap_df$cpg.id[coloc_nr]
            cpg_pos = overlap_df$cpg.pos[coloc_nr]

            gwas_signal_id_coloc = overlap_df$signal_gwas.rsid[coloc_nr]
            gwas_signal_chr = overlap_df$signal_gwas.chr[coloc_nr]
            gwas_signal_pos = overlap_df$signal_gwas.position[coloc_nr]
            
            gwas_signal_NR = as.character(overlap_df$signal_gwas.rsid[coloc_nr])
            ##############################
            # prepare mqtl_stats
            mqtl_stats = as.data.frame(in_m_qtl[[cpg_coloc]])
            
            mqtl_stats = mqtl_stats[,c("variant_id", "maf", "slope", "slope_se", "pval_nominal", "chr", "pos", "REF", "ALT" )]

            # NEA refers to REF allele
            mqtl_stats$mqtl_nea = mqtl_stats$REF
            
            # EA refers to ALT allele
            mqtl_stats$mqtl_ea = mqtl_stats$ALT

            # keep only SNPs
            nchar_EA = mqtl_stats$mqtl_ea %in% c("A", "T", "C", "G")
            nchar_NEA = mqtl_stats$mqtl_nea %in% c("A", "T", "C", "G")
            mqtl_stats = mqtl_stats[nchar_EA & nchar_NEA,]
                
         
               
            in_gwas$variant_id_chrpos <- paste(in_gwas$chr, in_gwas$position, toupper(in_gwas$ea), sep="_")
            mqtl_stats$variant_id_chrpos = paste(mqtl_stats$chr, mqtl_stats$pos, mqtl_stats$ALT, sep="_")
            
            mqtl_stats$varid_for_coloc <- mqtl_stats$variant_id_chrpos
            mqtl_stats = mqtl_stats[!is.na(mqtl_stats$varid_for_coloc),]
            mqtl_stats <- mqtl_stats %>% 
                                dplyr::filter(variant_id_chrpos %in% in_gwas$variant_id_chrpos)


	    if(dim(mqtl_stats)[1] < 100){
                return(NULL)
            }
            mqtl_stats$nvar_coloc = length(mqtl_stats$varid_for_coloc)


            print(paste0("All selected variants in gwassumstats: ", all(mqtl_stats$varid_for_coloc %in% 
                                                                        in_gwas$variant_id_chrpos)))
            
            #######################
            ## SAVE PWCOCO INPUT FILE
            
        		pwcoco_qtl_input <- mqtl_stats %>%
                                          dplyr::select("SNP"=variant_id, "A1"=ALT, "A2"=REF, "A1_freq"=maf, "beta"=slope, "se"=slope_se, "p"=pval_nominal) %>%
                                          dplyr::mutate("n"=QTL_n)
                          
            fwrite(pwcoco_qtl_input, file=file.path(out_path, paste0("PWCOCO/PWCOCO_mQTL_", tissue, "_", GWAS_ID), paste0(cpg_coloc,"_",gwas_signal_NR,"_mQTL_", tissue,"_GWAS_",GWAS_ID,".txt")), sep="\t")
            rm(pwcoco_qtl_input)
            
            
            pwcoco_gwas_input <- in_gwas %>%
                                    dplyr::filter(variant_id_chrpos %in%  mqtl_stats$varid_for_coloc) %>%
                                    dplyr::select("SNP"=rsid, "A1"=ea, "A2"=nea, "A1_freq"=eaf, "beta"=beta, "se"=se, "p"=p)

            if(GWAS_type == "cc"){
                pwcoco_gwas_input <-  pwcoco_gwas_input %>%
                                        dplyr::mutate("n"=sum(gwas_n),
                                              "s"=gwas_n[1]/sum(gwas_n))
            } else{
              pwcoco_gwas_input <-  pwcoco_gwas_input %>%
                                        dplyr::mutate("n"=gwas_n)
            }

            if(is.na(pwcoco_gwas_input$A1_freq[1])){
                pwcoco_gwas_input$A1_freq <- mqtl_stats$maf[match(mqtl_stats$variant_id, pwcoco_gwas_input$SNP)]
            }
            
            fwrite(pwcoco_gwas_input, file=file.path(out_path, paste0("PWCOCO/PWCOCO_", GWAS_ID), paste0(cpg_coloc,"_",gwas_signal_NR,"_mQTL_", tissue, "_GWAS_",GWAS_ID,".txt")), sep="\t")
            rm(pwcoco_gwas_input)
            
            
            
            #######################
            # From here
            mqtl_stats = mqtl_stats[, !(colnames(mqtl_stats) %in% 
                                        c("variant_id", "variant_id_rev", "mqtl_nea", "mqtl_ea"))]
            
            colnames(mqtl_stats) = paste0("MQTL_", colnames(mqtl_stats))
            ##############################
            gwas_stats = as.data.frame(in_gwas[in_gwas$variant_id_chrpos %in% mqtl_stats$MQTL_varid_for_coloc,])
            
            gwas_stats <- gwas_stats[,colSums(is.na(gwas_stats)) != nrow(gwas_stats)]
            
            colnames(gwas_stats) = paste0("GWAS_", colnames(gwas_stats))
            if(!is.null(gwas_n)){                                                                                                                                                 
                if(length(GWAS_n) != 1){                                                                                                                                          
                    gwas_stats$GWAS_N = sum(gwas_n)                                                                                                                               
                } else{                                                                                                                                                           
                    gwas_stats$GWAS_N = gwas_n                                                                                                                                    
                }                                                                                                                                                                 
                                                                                                                                                                                  
            }                                                                                                                                                                     
            
            mqtl_gwas_merge = merge(mqtl_stats, gwas_stats, by.x = "MQTL_varid_for_coloc", by.y = "GWAS_variant_id_chrpos")
            mqtl_gwas_merge <- mqtl_gwas_merge[!(duplicated(mqtl_gwas_merge$MQTL_varid_for_coloc)),]
            
            # Remove SNPs where the NEA is the same but the EA is different between QTL and GWAS
            mqtl_gwas_merge <- mqtl_gwas_merge[mqtl_gwas_merge$MQTL_ALT == toupper(mqtl_gwas_merge$GWAS_ea),]
            mqtl_gwas_merge <- mqtl_gwas_merge[mqtl_gwas_merge$MQTL_REF == toupper(mqtl_gwas_merge$GWAS_nea),]
            
            mqtl_gwas_merge <- mqtl_gwas_merge[complete.cases(mqtl_gwas_merge),]
            
            test_nea = identical(mqtl_gwas_merge$MQTL_REF, toupper(mqtl_gwas_merge$GWAS_nea))
            test_ea = identical(mqtl_gwas_merge$MQTL_ALT, toupper(mqtl_gwas_merge$GWAS_ea))
            ##############################
            # generate output
            n_gwas_variants_sign_5emins8 = nrow(subset(mqtl_gwas_merge, mqtl_gwas_merge$'GWAS_p' < 5e-8))
            n_gwas_variants_sign_1emins7 = nrow(subset(mqtl_gwas_merge, mqtl_gwas_merge$'GWAS_p' < 1e-7))
            n_gwas_variants_sign_1emins6 = nrow(subset(mqtl_gwas_merge, mqtl_gwas_merge$'GWAS_p' < 1e-6))
            n_gwas_variants_sign_1emins5 = nrow(subset(mqtl_gwas_merge, mqtl_gwas_merge$'GWAS_p' < 1e-5))
            n_gwas_variants_sign_1emins4 = nrow(subset(mqtl_gwas_merge, mqtl_gwas_merge$'GWAS_p' < 1e-4))
            ##############################
            if(GTEX_APP){
                ##############################
                # perform coloc.fast (this uses gtex implementation coloc.fast)
                coloc_model <- coloc.fast(beta1 = mqtl_gwas_merge$MQTL_slope_for_coloc, 
                                          se1 = mqtl_gwas_merge$MQTL_slope_se, 
                                          beta2 = mqtl_gwas_merge$GWAS_beta, se2 = mqtl_gwas_merge$GWAS_se)
                
                newline = data.frame(cpg = cpg_coloc, cpg_pos = cpg_pos, cpg_gene = as.character(m_anno[cpg_coloc,"gene"]),
                                     gwas_lead_snp = paste0(gwas_signal_chr, ":", gwas_signal_pos), 
                                     gwas_lead_snp_rs = gwas_signal_id_coloc, 
                                     nvariants = coloc_model$nvariants, 
                                     coloc_region_chr = unique(mqtl_gwas_merge$GWAS_chr), 
                                     coloc_region_start = ifelse(sign(gwas_signal_pos - 100000) == -1, 1, gwas_signal_pos - 100000), 
                                     coloc_region_end = gwas_signal_pos + 100000, 
                                     n_gwas_variants_sign_5emins8 = n_gwas_variants_sign_5emins8,
                                     n_gwas_variants_sign_1emins7 = n_gwas_variants_sign_1emins7, 
                                     n_gwas_variants_sign_1emins6 = n_gwas_variants_sign_1emins6, 
                                     n_gwas_variants_sign_1emins5 = n_gwas_variants_sign_1emins5, 
                                     n_gwas_variants_sign_1emins4 = n_gwas_variants_sign_1emins4, 
                                     PP0 = coloc_model$results$posterior[1], PP1 = coloc_model$results$posterior[2], 
                                     PP2 = coloc_model$results$posterior[3], PP3 = coloc_model$results$posterior[4],
                                     PP4 = coloc_model$results$posterior[5], test_nea = test_nea, test_ea = test_ea,    
                                     gwas_signal = gwas_signal_NR, stringsAsFactors = F)
                    ##############################
                    # save 
                    if(coloc_model$results$posterior[5] > 0.8 | (coloc_model$results$posterior[5] > 0.6 & coloc_model$results$posterior[5]/coloc_model$results$posterior[4] > 2.5)){
                      file_path = paste0(out_path, "sumstats_", cpg_coloc, "_", 
                                         as.character(paste0(gwas_signal_chr, ":", gwas_signal_pos,"_",gwas_signal_id_coloc)), ".rda")

                      save(mqtl_gwas_merge, file = file_path)
                    }
                   return(newline)
           }else{
                    ##############################
                    # another coloc model that enables 
                    # https://chr1swallace.github.io/coloc/articles/a03_enumeration.html

                    # mqtl
                    D1 = list(N = QTL_n, 
                              sdY=1,
                              MAF = mqtl_gwas_merge$MQTL_maf, 
                                    beta = mqtl_gwas_merge$MQTL_slope, 
                                    varbeta = (mqtl_gwas_merge$MQTL_slope_se)^2, 
                                    type = "quant", 
                                    snp = mqtl_gwas_merge$MQTL_variant_id_chrpos, 
                                    stringsAsFactors = F)
                    D1$MAF = ifelse(D1$MAF <= 0.5, D1$MAF, 1 - D1$MAF)    
                    # GWAS
                     D2 = list(
                              MAF = mqtl_gwas_merge$GWAS_eaf,                                                                                                                
                              sdY=1,                                                                                                                                             
                              type = GWAS_type,                                                                                                                                   
                                    beta = mqtl_gwas_merge$GWAS_beta,                                                                                                             
                                    varbeta = (mqtl_gwas_merge$GWAS_se)^2,                                                                                                        
                                    type = type,                                                                                                             
                                    snp = mqtl_gwas_merge$MQTL_varid_for_coloc, 
                                    stringsAsFactors = F)                                                                                                                         
                    if(GWAS_type == "cc"){                                                                                                                                        
                        D2$N = gwas_n[1] + gwas_n[2]                                                                                                                              
                        D2$s <- gwas_n[1]/(gwas_n[1] + gwas_n[2])                                                                                                                 
                    }else{                                                                                                                                                        
                        D2$N = gwas_n                                                                                                                                             
                                                                                                                                                                                  
                    }                                                                                                                                                             
                    D2$MAF = ifelse(D2$MAF <= 0.5, D2$MAF, 1 - D2$MAF)                                                                                                            
                              
                    # perform coloc
                    coloc_model <- coloc::coloc.abf(dataset1=D1, dataset2=D2, p12=1e-05)
                    
                    ##############################
                    if(coloc_model$summary["PP.H4.abf"] > 0.8 | (coloc_model$summary["PP.H4.abf"] > 0.6 & (coloc_model$summary["PP.H4.abf"]/coloc_model$summary["PP.H3.abf"]) > 2)){
                        
                        # save coloc frame
                        file_name = paste0("sumstats_", tissue, "_", cpg_coloc, "_", 
                                         as.character(paste0(gwas_signal_chr, ":", gwas_signal_pos,"_",gwas_signal_id_coloc)), ".rda")
                          save(mqtl_gwas_merge, file = file.path(out_path, paste0(GWAS_ID, "_coloc_rda_files"), file_name))
                        
                        
                     }
                        if(is.eQTL==TRUE){
                            # output credible set as data.table for GWAS eQTL colocalization
                            o <- order(coloc_model$results$SNP.PP.H4,decreasing=TRUE)
                            cs <- cumsum(coloc_model$results$SNP.PP.H4[o])
                            w <- which(cs > 0.95)[1]
                            credset <- coloc_model$results[o,][1:w,]$snp
                            lead.snp <- credset[1]
                            newline = data.frame(cpg = cpg_coloc, cpg_pos = cpg_pos,
                                                 coloc_lead_snp_chr = sub("_.*", "", lead.snp),
                                                 coloc_lead_snp_pos = strsplit(lead.snp, "_")[[1]][2],
                                                 coloc_lead_snp_rs = mqtl_gwas_merge[mqtl_gwas_merge$MQTL_variant_id_chrpos==lead.snp,]$GWAS_rsid,
                                                 credible_set = paste(mqtl_gwas_merge[mqtl_gwas_merge$MQTL_variant_id_chrpos %in% credset,]$GWAS_rsid, collapse=","),
                                                 gwas_lead_snp = paste0(gwas_signal_chr, ":", gwas_signal_pos), 
                                                 gwas_lead_snp_rs = gwas_signal_id_coloc, 
                                                 nvariants = coloc_model$summary["nsnps"], 
                                                 coloc_region_chr = unique(mqtl_gwas_merge$GWAS_chr), 
                                                 coloc_region_start = ifelse(sign(gwas_signal_pos - 1e6) == -1, 1, gwas_signal_pos - 1e6), 
                                                 coloc_region_end = gwas_signal_pos + 1e6, 
                                                 n_gwas_variants_sign_5emins8 = n_gwas_variants_sign_5emins8,
                                                 n_gwas_variants_sign_1emins7 = n_gwas_variants_sign_1emins7, 
                                                 n_gwas_variants_sign_1emins6 = n_gwas_variants_sign_1emins6, 
                                                 n_gwas_variants_sign_1emins5 = n_gwas_variants_sign_1emins5, 
                                                 n_gwas_variants_sign_1emins4 = n_gwas_variants_sign_1emins4, 
                                                 PP0 = coloc_model$summary["PP.H0.abf"], PP1 = coloc_model$summary["PP.H1.abf"], 
                                                 PP2 = coloc_model$summary["PP.H2.abf"], PP3 = coloc_model$summary["PP.H3.abf"],
                                                 PP4 = coloc_model$summary["PP.H4.abf"], test_nea = test_nea, test_ea = test_ea,   
                                                 gwas_signal = gwas_signal_NR, stringsAsFactors = F)
                        
                        return(newline)
             
                        } else{
                            o <- order(coloc_model$results$SNP.PP.H4,decreasing=TRUE)
                            cs <- cumsum(coloc_model$results$SNP.PP.H4[o])
                            w <- which(cs > 0.95)[1]
                            credset <- coloc_model$results[o,][1:w,]$snp
                            lead.snp <- credset[1]
                            newline = data.frame(cpg = cpg_coloc, cpg_pos = cpg_pos, cpg_gene = as.character(m_anno[cpg_coloc,"gene"]),
                                                 coloc_lead_snp_chr = sub("_.*", "", lead.snp),
                                                 coloc_lead_snp_pos = strsplit(lead.snp, "_")[[1]][2],
                                                 coloc_lead_snp_rs = mqtl_gwas_merge[mqtl_gwas_merge$MQTL_variant_id_chrpos==lead.snp,]$GWAS_rsid,
                                                 credible_set = as.list(credset),
                                                 gwas_lead_snp = paste0(gwas_signal_chr, ":", gwas_signal_pos), 
                                                 gwas_lead_snp_rs = gwas_signal_id_coloc, 
                                                 nvariants = coloc_model$summary["nsnps"], 
                                                 coloc_region_chr = unique(mqtl_gwas_merge$GWAS_chr), 
                                                 coloc_region_start = ifelse(sign(gwas_signal_pos - 100000) == -1, 1, gwas_signal_pos - 100000), 
                                                 coloc_region_end = gwas_signal_pos + 100000, 
                                                 n_gwas_variants_sign_5emins8 = n_gwas_variants_sign_5emins8,
                                                 n_gwas_variants_sign_1emins7 = n_gwas_variants_sign_1emins7, 
                                                 n_gwas_variants_sign_1emins6 = n_gwas_variants_sign_1emins6, 
                                                 n_gwas_variants_sign_1emins5 = n_gwas_variants_sign_1emins5, 
                                                 n_gwas_variants_sign_1emins4 = n_gwas_variants_sign_1emins4, 
                                                 PP0 = coloc_model$summary["PP.H0.abf"], PP1 = coloc_model$summary["PP.H1.abf"], 
                                                 PP2 = coloc_model$summary["PP.H2.abf"], PP3 = coloc_model$summary["PP.H3.abf"],
                                                 PP4 = coloc_model$summary["PP.H4.abf"], test_nea = test_nea, test_ea = test_ea,   
                                                 gwas_signal = gwas_signal_NR, stringsAsFactors = F)
                        
                        return(newline)
             
                        }
                              }
           
    }, mc.cores = cores))
    
    ##############################
    return(coloc_df)
}

coloc.fast <- function(beta1, se1, beta2, se2, 
                       priorsd1 = 1, priorsd2 = 1, priorc1 = 1e-4, priorc2 = 1e-4, priorc12 = 1e-5, 
                       rounded = 6) {
  abf1 <- abf.Wakefield(beta1, se1, priorsd1, log = TRUE)
  abf2 <- abf.Wakefield(beta2, se2, priorsd2, log = TRUE)
  w <- which(!is.na(abf1) & !is.na(abf2))
  ## Fill in zeros not filter out, thanks Karl and Karsten!!!!
  nv <- length(w)
  abf1 <- norm1(c(0, abf1[w]), log = TRUE)
  abf2 <- norm1(c(0, abf2[w]), log = TRUE)
  res <- data.frame(hypothesis = paste0("H", 0:4),
                    label = c("No association",
                      "One variant associated with phenotype 1 only",
                      "One variant associated with phenotype 2 only",
                      "Two variants separately associated with phenotypes 1 and 2",
                      "One variant associated with phenotypes 1 and 2"),
                    prior = norm1(c(1, priorc1*nv, priorc2*nv, priorc1*priorc2*nv*(nv - 1), priorc12*nv)),
                    bf = if (nv > 0) c(abf1[1]*abf2[1], 
                      sum(abf1[-1])*abf2[1]/nv, 
                      abf1[1]*sum(abf2[-1])/nv, 
                      (sum(abf1[-1])*sum(abf2[-1]) - sum(abf1[-1]*abf2[-1]))/(nv*(nv - 1)), 
                      sum(abf1[-1]*abf2[-1])/nv) else rep(NA, 5))
  res$bf <- res$bf/max(res$bf) # was: res$bf/res$bf[1] # caused division by zero for strongly colocalized signals
  res$posterior <- norm1(res$prior*res$bf)
  ## compute model averaged effect size ratios
  mw <- abf1[-1]*abf2[-1] # model weights
  alpha12 <- sum(beta1[w]/beta2[w]*mw)/sum(mw)
  alpha21 <- sum(beta2[w]/beta1[w]*mw)/sum(mw)
  if (is.finite(rounded)) {
    res$posterior = round(res$posterior, rounded)
  }
  return(list(results = res, nvariants = length(w), alpha12 = alpha12, alpha21 = alpha21))
}

# copied this function which is needed for coloc.fast package above
# https://rdrr.io/github/tobyjohnson/gtx/src/R/abf.R                    
abf.Wakefield <- function(beta, se, priorsd, log = FALSE) {
  if (log) {
    return(log(sqrt(se^2/(se^2 + priorsd^2))) + 
           (beta/se)^2/2 * priorsd^2/(se^2 + priorsd^2))
  } else {
    return(sqrt(se^2/(se^2 + priorsd^2)) * 
           exp((beta/se)^2/2 * priorsd^2/(se^2 + priorsd^2)))
  }
}           
                    
norm1 <- function(x, log = FALSE) {
  if (all(is.na(x))) return(x)
  if (log) {
    x <- x - max(x, na.rm = TRUE)
    x <- exp(x)    
  } else {
    ## This does not work if x contains NaNs or +Infs
    stopifnot(all(x >= 0, na.rm = TRUE))
    x <- x / max(x, na.rm = TRUE)
  }
  return(x / sum(x, na.rm = TRUE))
}                    

make_vcf <- function(GWASfile, chrom=NULL, pos=NULL, nea=Allele2, ea=Allele1, snp=SNP, ea_af=Freq1, effect=Effect, se=StdErr, pval=p, hg="hg38", GWAS_n=GWAS_n, variant_ann=variant_ann, lift_down=TRUE, WantToLiftOver=TRUE, output=NULL){
    print("Reading GWAS sumstats")
    vcf <- fread(GWASfile, data.table=FALSE)
    
    colnames(vcf)[colnames(vcf) == nea] <- "Allele2"
    colnames(vcf)[colnames(vcf) == ea] <- "Allele1"
   
    colnames(vcf)[colnames(vcf) == effect] <- "Effect"
    colnames(vcf)[colnames(vcf) == se] <- "StdErr"
    colnames(vcf)[colnames(vcf) == pval] <- "p"
    
       
    if(!is.null(pos) & !is.null(snp)){
        colnames(vcf)[colnames(vcf) == pos] <- "pos"
        colnames(vcf)[colnames(vcf) == chrom] <- "chrom"
        colnames(vcf)[colnames(vcf) == snp] <- "SNP"
        
        if(WantToLiftOver){
           print("Lifting over GWAS")
            res_ggr <- lifOverFunction_vcf(vcf)

            vcf$pos <- NA
            vcf$pos <- res_ggr$start[match(vcf$SNP, res_ggr$SNP)]
 
        }
        
    } else if(is.null(pos) & !is.null(snp)) {
        print("SNP position not provided, so using reference to map rsid to coordinates")
        
        colnames(vcf)[colnames(vcf) == snp] <- "SNP"
        
        vcf$chrom = variant_ann$chr[match(gsub(":.*","",vcf$SNP), variant_ann$rsid)]
        vcf$pos = variant_ann$pos[match(gsub(":.*","",vcf$SNP), variant_ann$rsid)]   
    } else if(!is.null(pos) & is.null(snp) & hg=="hg19"){
        print("SNP ID not provided, so using UKBB reference to map coordinates to rsid")
        
        colnames(vcf)[colnames(vcf) == pos] <- "pos"
        colnames(vcf)[colnames(vcf) == chrom] <- "chrom"
        
        vcf$SNP <- NA
        # Create an index for the matches
        indx <- match(paste0(vcf$chrom, vcf$pos, toupper(vcf$Allele2), toupper(vcf$Allele1)), paste0(variant_ann$chr, variant_ann$pos, variant_ann$ref, variant_ann$alt), nomatch = 0)
        vcf$SNP[indx != 0] <- variant_ann$rsid[indx]
        
        indx <- match(paste0(vcf$chrom, vcf$pos, toupper(vcf$Allele1), toupper(vcf$Allele2)), paste0(variant_ann$chr, variant_ann$pos, variant_ann$ref, variant_ann$alt), nomatch = 0)
        vcf$SNP[indx != 0] <- variant_ann$rsid[indx]
        
           
    }
        
    if(!is.null(ea_af)){
         colnames(vcf)[colnames(vcf) == ea_af] <- "Freq1"
    } else {
        print("Effect allele frequency not provided, assuming the same frequency as UKBB")
        vcf$Freq1 <- NA
        # Create an index for the matches
        indx <- match(paste0(vcf$chrom, vcf$pos, toupper(vcf$Allele2), toupper(vcf$Allele1)), paste0(variant_ann$chr, variant_ann$pos, variant_ann$ref, variant_ann$alt), nomatch = 0)
        vcf$Freq1[indx != 0] <- variant_ann$AF[indx]
        
        indx <- match(paste0(vcf$chrom, vcf$pos, toupper(vcf$Allele1), toupper(vcf$Allele2)), paste0(variant_ann$chr, variant_ann$pos, variant_ann$ref, variant_ann$alt), nomatch = 0)
        vcf$Freq1[indx != 0] <- 1 - variant_ann$AF[indx]
        
    }
    
   
    
    print("Filtering VCF for MAF 5% and known SNPs")
    vcf <- vcf %>% 
        dplyr::filter(!is.na(pos)) %>%
        dplyr::filter(!is.na(Freq1)) %>%
        dplyr::filter(!is.na(SNP)) %>%
        dplyr::filter(Freq1 > 0.05 & Freq1 < 0.95) %>%
        dplyr::filter(toupper(Allele1) %in% c("A","T","C","G")) %>%
        dplyr::mutate(SNP=gsub(":.*","",SNP),
                     p=as.numeric(p))
    
    print(paste0("Creating VCF with ", dim(vcf)[1], " SNPs"))
    
    if(length(GWAS_n) == 1){
        vcf$SS <- GWAS_n
        out <- vcf %$%
        create_vcf(chrom=chrom, pos=pos, nea=Allele2, ea=Allele1, snp=SNP, ea_af=Freq1, effect=Effect, se=StdErr, pval=p, n=SS, name="a")

        
    } else{
        vcf$SS = sum(GWAS_n)
        vcf$ncase = GWAS_n[1]
        print(colnames(vcf))
        out <- vcf %$%
        create_vcf(chrom=chrom, pos=pos, nea=Allele2, ea=Allele1, snp=SNP, ea_af=Freq1, effect=Effect, se=StdErr, pval=p, n=SS, ncase=ncase, name="a")

    }
    
   if(is.null(output)){
        dir.create("GWAS_data")
        fname <- gsub(".gz|.txt",".vcf",basename(GWASfile))
        
	output=file.path("GWAS_data", fname)
	}
    

    writeVcf(out, file=output, index=TRUE)


    
}

lifOverFunction_vcf <- function(regions, lift_down=TRUE){
	if(lift_down){
        	ch_path = "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/FunctionsAndData/hg38ToHg19.over.chain"
	  } else {
        	ch_path = "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/FunctionsAndData/hg19ToHg38.over.chain"
    	}

    ch = import.chain(ch_path)

  #print(region)
  res_ggr <- makeGRangesFromDataFrame(regions,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=TRUE,
                                    seqinfo=NULL,
                                    seqnames.field=c("chrom"),
                                    start.field="pos",
                                    end.field="pos",
                                    strand.field=NULL,
                                    starts.in.df.are.0based=FALSE)

  seqlevelsStyle(res_ggr) = "UCSC" # necessary
  res_ggr = liftOver(res_ggr, ch)
  res_ggr<-unlist(res_ggr)
  res_ggr<-as.data.table(res_ggr)
}

GWAS_sumstats_extract <- function(vcfname, top, window_length=1e6, output_path="GWAS_data"){
    GWAS_loci <- paste0(top$chr, ":", 
                     ifelse((top$position - window_length) < 0, 1, top$position - window_length) , 
                     "-", 
                     top$position + window_length)

    GWAS_associations <- c()
    x=1
    for(GWAS_locus in GWAS_loci){
        print(paste("Locus", x, "of", length(GWAS_loci)))
        param <- ScanVcfParam(which=GRanges(GWAS_locus) )
        vcf <- readVcf(vcfname, "hg19", param)

        GWAS_association <-  gwasvcf::vcf_to_granges(vcf) %>% dplyr::as_tibble() 
        
        if(dim(GWAS_association)[1] > 0){
        GWAS_association <- GWAS_association %>%
            dplyr::select(seqnames, start, ES, SE, LP, SS, id, ID, ALT, REF, AF) %>%
            dplyr::rename(
                chr = seqnames,
                position = start,
                rsid = ID,
                ea = ALT,
                nea = REF,
                eaf = AF,
                beta = ES,
                se = SE,
                p = LP,
                n = SS, 
                id = 
            ) %>%
        dplyr::mutate(p = 10^(-p) )

        GWAS_associations <- rbind(GWAS_associations, GWAS_association %>% 
                              dplyr::filter(!duplicated(rsid)))
            }
        x=x+1
        }

   
    GWAS_associations <- unique(GWAS_associations)
    
    return(GWAS_associations)

    save(GWAS_associations, file = file.path(output_path, paste0(GWAS_ID,".rda")))

} 

ld_matrix_local_GO <- function(GWAS, chr, fname, bfile=paste0("/lustre/groups/itg/shared/referenceData/ukbiobank/chip/bgen2plink/chr",chr)){
  ## write a table with SNP names
  write.table(GWAS[,c('chr.pos.b37','A2')], paste0('Coloc_susie/LD/',gsub("_mQTL_GWAS.*","", fname),'.SNP.REF.txt'), col.names=F, row.names=F, sep='\t', quote=F)
  write.table(GWAS[,c('chr.pos.b37')], paste0('Coloc_susie/LD/',gsub("_mQTL_GWAS.*","", fname),'.SNP.extract'), col.names=F, row.names=F, sep='\t', quote=F)

  system(paste0('cut -f 2 ', bfile,'.bim | sort | uniq -d > ukb_ld/', chr, '.dups'))
  system(paste0('plink --bfile ', bfile, ' --exclude ukb_ld/', chr, '.dups --const-fid --extract Coloc_susie/LD/',gsub("_mQTL_GWAS.*","", fname),'.SNP.extract --threads 1 --update-ref-allele Coloc_susie/LD/',gsub("_mQTL_GWAS.*","", fname),'.SNP.REF.txt 2 1 --r square --out Coloc_susie/LD/',gsub("_mQTL_GWAS.*","", fname),'.LD --make-just-bim'))

  ## read in the resulting plink files
  bim = read.table(paste0('Coloc_susie/LD/',gsub("_mQTL_GWAS.*","", fname),'.LD.bim'), header=F)
  #give them column names
  colnames(bim) = c("chr","rsid","cM","position","allele1","allele2")
  #read in the LD matrix (all numeric)
  ld = as.matrix(read.table(paste0('Coloc_susie/LD/',gsub("_mQTL_GWAS.*","", fname),'.LD.ld'), header=F, stringsAsFactors=F))
  #give the row and column names of this matrix the rsIDs
  colnames(ld) = bim[,2]
  rownames(ld) = bim[,2]

  ld_subset <- ld[rownames(ld) %in% GWAS$chr.pos, colnames(ld) %in% GWAS$chr.pos]
  return(ld_subset)
}

ld_matrix_local_MT <- function(QTL, chr, fname, bfile=paste0("/lustre/groups/itg/teams/zeggini/projects/child_diabesity/analysis/methyl/matrixeQTL/Results/plinkGenotypeFiles/chr",chr)){
  ## write a table with SNP names
  write.table(QTL[,c('ID','A2')], paste0('Coloc_susie/LD/',gsub("_mQTL_GWAS.*","", fname),'.SNP.REF.txt'), col.names=F, row.names=F, sep='\t', quote=F)
  write.table(QTL[,c('ID')], paste0('Coloc_susie/LD/',gsub("_mQTL_GWAS.*","", fname),'.SNP.extract'), col.names=F, row.names=F, sep='\t', quote=F)

  system(paste0('plink --bfile ', bfile, ' --const-fid --extract Coloc_susie/LD/',gsub("_mQTL_GWAS.*","", fname),'.SNP.extract --threads 1 --update-ref-allele Coloc_susie/LD/',gsub("_mQTL_GWAS.*","", fname),'.SNP.REF.txt 2 1 --r square --out Coloc_susie/LD/',gsub("_mQTL_GWAS.*","", fname),'.LD --make-just-bim'))

  ## read in the resulting plink files
  bim = read.table(paste0('Coloc_susie/LD/',gsub("_mQTL_GWAS.*","", fname),'.LD.bim'), header=F)
  #give them column names
  colnames(bim) = c("chr","rsid","cM","position","allele1","allele2")
  #read in the LD matrix (all numeric)
  ld = as.matrix(read.table(paste0('Coloc_susie/LD/',gsub("_mQTL_GWAS.*","", fname),'.LD.ld'), header=F, stringsAsFactors=F))
  #give the row and column names of this matrix the rsIDs
  colnames(ld) = bim[,2]
  rownames(ld) = bim[,2]

  ld_subset <- ld[rownames(ld) %in% QTL$ID, colnames(ld) %in% QTL$ID]
  return(ld_subset)
}



