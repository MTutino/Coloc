library(data.table)
setwd("C:/Users/AnaLuizaArruda/OneDrive - Helmholtz Zentrum München/Projects/molQTL300/")

targets <- as.data.table(readRDS(file="OT_move_to_lustre/data_parsed/OT_targets.rds", refhook=NULL))
colnames(targets) <- c("Gene", "ID")
# query=as.data.frame(readLines("analysis/query.tsv"))
# query=as.vector(query[,1])
query <- fread("coloc/colocalizing_genes_with_names.csv")

msg=c("--- How many genes and drugs were linked before applying our (search) conditions to the results table ? Overall listings ? ---",
      paste( "All drug-gene listings: ", length(targets[,ID]) ),
      paste( length( unique(targets[,Gene])), " unique genes" ),
      paste( length( unique(targets[,ID])), " unique drugs" )
)
print(msg)
# capture.output(msg, file=outfile, append=TRUE)


print("--- Selecting approved ---")
targets=targets[targets$drug_isApproved==TRUE,]
print("--- Selecting not withdrawn ---")
targets=targets[targets$drug_hasBeenWithdrawn==FALSE,]
print("--- Selecting listing indications ---")
targets=targets[is.na(targets$indications)==FALSE,]

msg=c(  "--- Selecting approved, not withdrawn, listing indications ---",
        paste0( "Remaining drug-gene listings: ", length(targets[,"ID"]) ),
        paste0( length(unique(targets[,"Gene"])), " unique genes" ),
        paste0( length(unique(targets[,"ID"])), " unique drugs" )
)

OT <- as.data.table(readRDS(file="OT_move_to_lustre/analysis/OT_by_gene.rds", refhook=NULL))
OT[ENSG_id %in% query$gene, .N]
OT[ENSG_id %in% query$gene]

# Add MR information
drug.dt <- fread("druggability_lookup.csv")
mr <- fread("MR/GO2_sig_processed_results_cis.csv")
all.mr <- fread("MR/GO2_processed_results_cis.csv")

new.drug.dt <- merge(drug.dt, unique(mr[gene %in% drug.dt$ENSG_id, .(gene, mr.direction=sign(beta))]), by.x="ENSG_id", by.y="gene", all.x=TRUE)
fwrite(new.drug.dt, "druggability_lookup_with_mr.csv")
