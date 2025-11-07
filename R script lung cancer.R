# تثبيت الحزمة (لو أول مرة)
if (!requireNamespace("GEOquery", quietly = TRUE))
  install.packages("BiocManager")

# تحميل GEOquery من Bioconductor
BiocManager::install("GEOquery")
# تحميل البيانات
gset <- getGEO("GSE19804", GSEMatrix = TRUE)

# تأكد كام عنصر جوه القائمة
length(gset)
gset <- gset[[1]]
expr_matrix <- exprs(gset)
metadata <- pData(gset)
head(metadata)
colnames(metadata)
"tissue:ch1"
table(metadata[,"tissue:ch1"])
group_list <- ifelse(metadata[,"tissue:ch1"] == "lung cancer", "Tumor", "Normal")
group_list <- factor(group_list)
table(group_list)
library(limma)
design <- model.matrix(~ 0 + group_list)
colnames(design) <- levels(group_list)
design
contrast.matrix <- makeContrasts(Tumor - Normal, levels=design)
fit <- lmFit(expr_matrix, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg_results <- topTable(fit2, adjust="fdr", number=Inf)
head(deg_results)
library(EnhancedVolcano)
EnhancedVolcano(deg_results,
                lab = rownames(deg_results),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = 'Volcano Plot: Lung Cancer vs Normal',
                subtitle = 'GSE19804',
                caption = 'Source: GEO')
BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
gene_symbols <- mapIds(hgu133plus2.db,
                       keys=rownames(deg_results),
                       column="SYMBOL",
                       keytype="PROBEID",
                       multiVals="first")
deg_results$GeneSymbol <- gene_symbols
EnhancedVolcano(deg_results,
                lab = deg_results$GeneSymbol,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = 'Volcano Plot: Lung Cancer vs Normal',
                subtitle = 'GSE19804',
                caption = 'Source: GEO')
library(ggrepel)

# اختيار أهم 20 Upregulated و 20 Downregulated
top_genes <- deg_gene %>%
  filter(color != "NotSignificant") %>%
  arrange(desc(-log10(adj.P.Val))) %>%
  group_by(color) %>%
  slice_head(n = 20) %>%
  ungroup()

ggplot(deg_gene, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "NotSignificant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, aes(label = GeneSymbol),
                  size = 3, box.padding = 0.3, max.overlaps = 20) +
  theme_minimal() +
  labs(title = "Volcano Plot: Lung Cancer vs Normal",
       x = "log2 Fold Change",
       y = "-log10(adj.P.Val)") +
  theme(legend.title = element_blank())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("STRINGdb")
library(STRINGdb)
library(dplyr)

# إنشاء كائن STRINGdb للبشر
string_db <- STRINGdb$new(version="11", species=9606, score_threshold=400, input_directory="")

# دمج مع DEGs
deg_top <- deg_gene %>% 
  arrange(adj.P.Val) %>% 
  slice_head(n = 50)   # أهم 50 DEG

deg_mapped <- string_db$map(deg_top, "GeneSymbol", removeUnmappedRows = TRUE)

# استخراج الشبكة
ppi_network <- string_db$get_subnetwork(deg_mapped$STRING_id)

# عرض الشبكة
string_db$plot_network(ppi_network)













