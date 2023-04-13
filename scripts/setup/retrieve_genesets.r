require(msigdbr)
require(biomaRt)

# Import the following gene sets from C2:
# https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_TRANSCRIPTIONAL_REGULATION_OF_PLURIPOTENT_STEM_CELLS.html
# https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/WP_WNT_SIGNALING_PATHWAY_AND_PLURIPOTENCY.html

get_geneset_coords <- function(mart, genesets, categ) {
    hsap.geneset.list <- lapply(categ, msigdbr, species = "Homo sapiens")
    hsap.genes <- purrr::reduce(hsap.geneset.list, rbind)
    hsap.genes <- hsap.genes %>%
                    filter(gs_name %in% genesets)
    hsap.genes <- dplyr::select(hsap.genes, human_ensembl_gene, human_entrez_gene, human_gene_symbol)

    # Get GRCh37 coordinates for each Ensembl ID - https://www.biostars.org/p/136775/
    gene.coords <- getBM(attributes = c("chromosome_name", "start_position", "end_position",
                                         "external_gene_name", "strand"),
                          filters = "ensembl_gene_id",
                          values = hsap.genes$human_ensembl_gene,
                          mart = mart)

    # Convert strand to canonical +/-, remove non-canonical alt splice variants
    gene.coords <- gene.coords %>%
                    mutate(gene_strand = case_when(
                        strand == 1 ~ "+",
                        strand == -1 ~ "-"
                    )) %>%
                    dplyr::select(-strand)
    return(gene.coords)
}

# For host, use official GRCh37 Ensembl database: https://grch37.ensembl.org/index.html
h37.mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org", dataset = "hsapiens_gene_ensembl")



#------PLURIPOTENCY AND WNT/B-CATENIN GENES-----
# Genes associated with pluripotency and Wnt/beta catenin
genesets <- c("WP_WNT_SIGNALING_PATHWAY_AND_PLURIPOTENCY",
              "REACTOME_TRANSCRIPTIONAL_REGULATION_OF_PLURIPOTENT_STEM_CELLS")

# Get gene coordinates
pluri.genes <- get_geneset_coords(h37.mart, genesets, categ = "C2")
# Filter out the following known tumor suppressors from the Wnt pathway:
#   - protein phosphatases (PP2A, AXIN, APC, GSK are part of the beta-catenin destruction complex)
#   - AXIN
#   - APC
#   - GSK
#   - RAC GTPase activating proteins
#   - TP53
#   - Naked cuticle genes, which bind Dvl and prevent its recruitment to Fzd
tumor.supp <- "AXIN.*|APC|^PP.*|GSK.*|RACG.*|TP53|NKD.*"
is.tumorsupp <- grepl(x = pluri.genes$external_gene_name, pattern = tumor.supp)
pluri.genes <- pluri.genes[!is.tumorsupp, ]

# Remove non-canonical alt splice variants on scaffolds
pluri.coords <- filter(pluri.genes, chromosome_name %in% c(seq(1, 22, 1), "X", "Y"))

# Write to .bed file
write.table(pluri.coords, file = "pluri_genes.bed", quote = FALSE,
            eol = "\n", row.names = FALSE, col.names = FALSE, sep = "\t")



#------KRAS PATHWAY--------
mapk.geneset <- c("KEGG_MAPK_SIGNALING_PATHWAY")
mapk.genes <- get_geneset_coords(h37.mart, mapk.geneset, "C2") %>%
                filter(chromosome_name %in% c(seq(1, 22, 1), "X", "Y"))

# Write to .bed
write.table(mapk.genes, file = "mapk_genes.bed", quote = FALSE,
            eol = "\n", row.names = FALSE, col.names = FALSE, sep = "\t")
