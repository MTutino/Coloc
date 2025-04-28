project.path <- "/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/"
setwd(project.path)

library(data.table)
library(dplyr)
library(httr)
library(jsonlite)

top.genes <- fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/all_colocalizing_genes_with_names.csv")

# Add OMIM lookup
search.OMIM <- function(gene) {
  stopifnot(is.character(gene), length(gene)==1)
  gene <- tolower(gene)
  server <- "https://api.omim.org"
  r <- httr::GET(paste(server, "/api/geneMap/search?search=", gene, "&include=geneMap&include=clinicalSynopsis&apiKey=tZtgP9bNTGCcLOULBHLl7A&format=json&start=0&limit=100", sep = ""), httr::content_type("application/json"))
  httr::stop_for_status(r)
  
  content_json <- fromJSON(toJSON(content(r)))
  
  # Extract phenotype map list
  gene_maps <- content_json$omim$searchResponse$geneMapList
  if (length(gene_maps) == 0) {
    return(NA)
  }
  
  mim.numbers <- unlist(content_json$omim$searchResponse$geneMap[[1]]$phenotypeMapList[[1]]$phenotypeMap$phenotypeMimNumber)
  if(is.null(mim.numbers)) return(NA)
  r_clinical <- httr::GET(paste0("https://api.omim.org/api/clinicalSynopsis?mimNumber=", paste(mim.numbers, collapse=","), "&apiKey=tZtgP9bNTGCcLOULBHLl7A&format=json&include=all"))
  a <- tryCatch({
    httr::stop_for_status(r_clinical)},
    error=function(e){
      return("error")}  
  )
  if(class(a)!="response") return(NA)
  content_json_clinical <- fromJSON(toJSON(content(r_clinical)))
  
  # Extract phenotype map list
  clinical_features <- as.data.table(content_json_clinical$omim$clinicalSynopsisList[[1]])
  if (nrow(clinical_features) == 0) {
    return(NA)
  }
  
  possible.columns <- c("growth", "growthHeight", "skeletal", "skeletalSpine", "skeletalLimbs", "skeletalHands", "muscleSoftTissue")
  if(sum(colnames(clinical_features) %in% possible.columns)>1){
    valid.columns <- c("preferredTitle", possible.columns[possible.columns %in%  colnames(clinical_features)])
    res.dt <- clinical_features[, ..valid.columns]
    res.dt <- res.dt[, preferredTitle:=unlist(preferredTitle)]
    res.dt$gene <- toupper(gene)
    return(res.dt)
  }
}

lookup.OMIM <- function(genes) {
  genes[hgnc_symbol=="", OMIM:=0]
  omim.search <- sapply(genes[hgnc_symbol!="", hgnc_symbol], search.OMIM)
  save.output <- omim.search[!is.na(omim.search)]
  save.output <- save.output[!sapply(save.output, is.null)]
  save(save.output, file="OMIM_phenotypes.rda")
  genes[, OMIM:=ifelse(hgnc_symbol %in% names(save.output), 1, 0)]

  return(genes)
}

coloc.genes <- fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/all_colocalizing_genes_with_names.csv")
res.genes <- lookup.OMIM(genes=coloc.genes)

# Add KO mice lookup
ko.mice <- data.table::fread("BatchReport_20250227_052840.txt", sep="\t", fill=T) %>% .[, .(Input, Term)]
ko.mice <- ko.mice[!is.na(Term)]
oa.mice.terms <- c("bone", "muscle", "skeleton", "osteo", "arthritis", "muscular", "joint", "body size", "growth", "skeletal", "stature", 
                   "height", "limb", "appendage", "abnormal cartilage", "articular cartilage", "chondrocyte", "body length", "brachypodia",
                   "brachydactyly", "brachyphalangia", "femur", "tibia", "ulna", "fibula", "humerus", "radial", "spine curvature", "posture",
                   "vertebral arch", "syndactyly")
res.genes[, OA.ko.mice:=ifelse(hgnc_symbol %in% ko.mice[(ko.mice[, grepl(paste(oa.mice.terms, collapse="|"), Term, ignore.case=TRUE)]), Input], 1, 0)]

# Add MR results
mr.genes.dt <- fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/GO2_coloc_mr_genes.txt")
res.genes[, MR:=ifelse(ensembl_gene_id %in% mr.genes.dt$gene.id, 1, 0)]

# Add overlap with HiC data
hic.enhancer <- fread("overlap_credset_enhancer.csv")
hic.promoter <- fread("HiCannotation/unique_result_Chondrocyte_active_promoter.txt")
res.genes[, `:=` (promoter_enhancer_HiC=ifelse(hgnc_symbol %in% hic.enhancer$Gene.symbol, 1, 0), active_promoter_HiC=ifelse(hgnc_symbol %in% hic.promoter$gene.name, 1, 0))]

# Add colocalization information
coloc.results <- fread("/lustre/groups/itg/teams/zeggini/projects/fungen-oa/analyses/Ana_coloc_mr/coloc/coloc_per_gene_tissue_trait.csv")
res.genes[, coloc.HG:=as.integer(ifelse(ensembl_gene_id %in% unique(coloc.results[tissue=="high_grade_cartilage", gene]), 1, 0))]
res.genes[, coloc.LG:=as.integer(ifelse(ensembl_gene_id %in% unique(coloc.results[tissue=="low_grade_cartilage", gene]), 1, 0))]
res.genes[, coloc.Syn:=as.integer(ifelse(ensembl_gene_id %in% unique(coloc.results[tissue=="synovium", gene]), 1, 0))]
res.genes[, coloc.FP:=as.integer(ifelse(ensembl_gene_id %in% unique(coloc.results[tissue=="fat_pad", gene]), 1, 0))]

# Add druggability results
drugs.dt <- fread("druggability_lookup.csv")
drugs.dt[grep("osteoarthritis|arthritis|pain", indications)]
res.genes[, drug:=as.integer(0)]
res.genes[ensembl_gene_id %in% drugs.dt$ENSG_id, drug:=as.integer(1)]

# Add GO2 prioritized genes
go2.genes <- as.data.table(readxl::read_xlsx("go2_effector_genes.xlsx", skip=2))
go2.genes[SCORE>=3, .N]
go2.coloc.genes <- go2.genes[!is.na(`eQTL colocalization Synovium`) | !is.na(`eQTL colocalization chondrocyte low-grade`) | !is.na(`eQTL colocalization chondrocyte high-grade`), ENSG]
res.genes[!(ensembl_gene_id %in% go2.coloc.genes), .N]
res.genes[!(ensembl_gene_id %in% go2.genes[SCORE>=3, ENSG]), .N]
res.genes[ensembl_gene_id %in% go2.genes[SCORE>=3, ENSG], GO2.gene:=1]

res.genes[is.na(drug), drug:=0]
res.genes[is.na(GO2.gene), GO2.gene:=0]
fwrite(res.genes, "GO2_prioritized_genes_final.csv")

