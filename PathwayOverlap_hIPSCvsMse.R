#FGSEA
library(fgsea)
library(data.table)
library(ggplot2)
library(tidyr)
library(KEGGREST)
library(readr)
library(tidyverse)

Astrocytes_Mono_VS_control_GSE191248_ <- read_csv("Astrocytes  (Mono) VS control(GSE191248).csv")
SCZVSControl_CorticolNeu_GSE182875_ <- read_csv("SCZVSControl(CorticolNeu)(GSE182875).csv")
With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_ <- read_csv("With SCZ VS Without SCZ (Astrocytes)(GSE191248).csv")
Without_SCZ_and_Control_Astrocytes_GSE191248_ <- read_csv("Without SCZ and Control(Astrocytes)(GSE191248).csv")
Mouse_Astrocyte_E2_knockdown <- read_csv("Mouse_Astrocyte_E2_knockdown.csv")
Mouse_neuron_E2_knockdown <- read_csv("Mouse_neuron_E2_knockdown.csv")

#HUMAN
#astrocyte
upreg_Astrocytes_Mono_VS_control_GSE191248_ <- 
  select(Astrocytes_Mono_VS_control_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue < 0.05, log2FoldChange > 0) %>% 
  arrange(log2FoldChange)
upreg_Astrocytes_Mono_VS_control_GSE191248_
downreg_Astrocytes_Mono_VS_control_GSE191248_ <- 
  select(Astrocytes_Mono_VS_control_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue < 0.05, log2FoldChange < 0) %>% 
  arrange(log2FoldChange)
downreg_Astrocytes_Mono_VS_control_GSE191248_
upreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_ <- 
  select(With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue < 0.05, log2FoldChange > 0) %>% 
  arrange(log2FoldChange)
upreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_
downreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_ <- 
  select(With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue < 0.05, log2FoldChange < 0) %>% 
  arrange(log2FoldChange)
downreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_
upreg_Without_SCZ_and_Control_Astrocytes_GSE191248_ <- 
  select(Without_SCZ_and_Control_Astrocytes_GSE191248_ , Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue < 0.05, log2FoldChange > 0) %>% 
  arrange(log2FoldChange)
upreg_Without_SCZ_and_Control_Astrocytes_GSE191248_ 
downreg_Without_SCZ_and_Control_Astrocytes_GSE191248_  <- 
  select(Without_SCZ_and_Control_Astrocytes_GSE191248_ , Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue < 0.05, log2FoldChange < 0) %>% 
  arrange(log2FoldChange)
downreg_Without_SCZ_and_Control_Astrocytes_GSE191248_ 
#neuron
upreg_SCZVSControl_CorticolNeu_GSE182875_ <- 
  select(SCZVSControl_CorticolNeu_GSE182875_, Symbol, pvalue, log2FoldChange) %>%
  filter(pvalue < 0.05, log2FoldChange > 0) %>% 
  arrange(log2FoldChange)
upreg_SCZVSControl_CorticolNeu_GSE182875_
downreg_SCZVSControl_CorticolNeu_GSE182875_ <- 
  select(SCZVSControl_CorticolNeu_GSE182875_, Symbol, pvalue, log2FoldChange) %>%
  filter(pvalue < 0.05, log2FoldChange < 0) %>% 
  arrange(log2FoldChange)
downreg_SCZVSControl_CorticolNeu_GSE182875_
#MOUSE
#astrocyte
upreg_Mouse_Astrocyte_E2_knockdown <- 
  select(Mouse_Astrocyte_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue < 0.05, !is.na(hsapiens_homolog_associated_gene_name), log2FoldChange > 0) %>% 
  arrange(log2FoldChange)
upreg_Mouse_Astrocyte_E2_knockdown
downreg_Mouse_Astrocyte_E2_knockdown <-
  select(Mouse_Astrocyte_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue < 0.05, !is.na(hsapiens_homolog_associated_gene_name), log2FoldChange < 0) %>% 
  arrange(log2FoldChange)
downreg_Mouse_Astrocyte_E2_knockdown
#neuron
upreg_Mouse_neuron_E2_knockdown <- 
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue < 0.05, !is.na(hsapiens_homolog_associated_gene_name), log2FoldChange > 0) %>% 
  arrange(log2FoldChange)
upreg_Mouse_neuron_E2_knockdown
downreg_Mouse_neuron_E2_knockdown <-
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue < 0.05, !is.na(hsapiens_homolog_associated_gene_name), log2FoldChange < 0) %>% 
  arrange(log2FoldChange)
downreg_Mouse_neuron_E2_knockdown

SCZ_vs_HC_up_HuAst <- matrix(upreg_Astrocytes_Mono_VS_control_GSE191248_$Symbol, ncol =1)
SCZ_vs_HC_up_HuAst
SCZ_vs_HC_down_HuAst <- matrix(downreg_Astrocytes_Mono_VS_control_GSE191248_$Symbol, ncol =1)
SCZ_vs_HC_down_HuAst
SCZ_vs_nonSCZ_up_HuAst <- matrix(upreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_$Symbol, ncol =1)
SCZ_vs_nonSCZ_up_HuAst
SCZ_vs_nonSCZ_down_HuAst <- matrix(downreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_$Symbol, ncol =1)
SCZ_vs_nonSCZ_down_HuAst
nonSCZ_vs_HC_up_HuAst <- matrix(upreg_Without_SCZ_and_Control_Astrocytes_GSE191248_$Symbol, ncol =1)
nonSCZ_vs_HC_up_HuAst
nonSCZ_vs_HC_down_HuAst <- matrix(downreg_Without_SCZ_and_Control_Astrocytes_GSE191248_$Symbol, ncol =1)
nonSCZ_vs_HC_down_HuAst
SCZ_vs_HC_up_HuN <- matrix(upreg_SCZVSControl_CorticolNeu_GSE182875_$Symbol, ncol =1)
SCZ_vs_HC_up_HuN
SCZ_vs_HC_down_HuN <- matrix(downreg_SCZVSControl_CorticolNeu_GSE182875_$Symbol, ncol =1)
SCZ_vs_HC_down_HuN

E2KO_up_MseAst <- matrix(upreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name, ncol =1)
E2KO_up_MseAst
E2KO_down_MseAst <- matrix(downreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name, ncol =1)
E2KO_down_MseAst
E2KO_up_MseN <- matrix(upreg_Mouse_neuron_E2_knockdown$hsapiens_homolog_associated_gene_name, ncol =1)
E2KO_up_MseN
E2KO_down_MseN <- matrix(downreg_Mouse_neuron_E2_knockdown$hsapiens_homolog_associated_gene_name, ncol =1)
E2KO_down_MseN

E2KO_MseAst_pathway <- read_csv("E2KO_MseAst_pathway.csv")
E2KO_MseN_pathway <- read_csv("E2KO_MseN_pathway.csv")
nonSCZ_vs_HC_HuAst_pathway <- read_csv("nonSCZ_vs_HC_HuAst_pathway.csv")
SCZ_vs_HC_HuAst_pathway <- read_csv("SCZ_vs_HC_HuAst_pathway.csv")
SCZ_vs_nonSCZ_HuAst_pathway <- read_csv("SCZ_vs_nonSCZ_HuAst_pathway.csv")
SCZ_vs_HC_HuN_pathway <- read_csv("SCZ_vs_HC_HuN_pathway.csv")

#Adjusted pvalue used because Pvalue yields too many pathways
E2KO_MseAst_pathway_filtered <- E2KO_MseAst_pathway %>% 
  filter(`Adjusted P-value` < 0.05)
E2KO_MseAst_pathway_filtered
E2KO_MseN_pathway_filtered <- E2KO_MseN_pathway %>%
  filter(`Adjusted P-value` < 0.05)
E2KO_MseN_pathway_filtered
nonSCZ_vs_HC_HuAst_pathway_filtered <- nonSCZ_vs_HC_HuAst_pathway %>%
  filter(`Adjusted P-value` < 0.05)
nonSCZ_vs_HC_HuAst_pathway_filtered
SCZ_vs_HC_HuAst_pathway_filtered <- SCZ_vs_HC_HuAst_pathway %>%
  filter(`Adjusted P-value` < 0.05)
SCZ_vs_HC_HuAst_pathway_filtered
SCZ_vs_nonSCZ_HuAst_pathway_filtered <- SCZ_vs_nonSCZ_HuAst_pathway %>%
  filter(`Adjusted P-value` < 0.05)
SCZ_vs_nonSCZ_HuAst_pathway_filtered
SCZ_vs_HC_HuN_pathway_filtered <- SCZ_vs_HC_HuN_pathway %>%
  filter(`Adjusted P-value` < 0.05) 
SCZ_vs_HC_HuN_pathway_filtered

pvalue_columns <- c("P-value.x", "P-value.y","P-value.x.x","P-value.y.y","P-value.x.x.x","P-value.y.y.y")

#overlapping of similar pathways
Pathways <- E2KO_MseAst_pathway_filtered %>% 
  full_join(E2KO_MseN_pathway_filtered, by = join_by(Term)) %>% 
  full_join(nonSCZ_vs_HC_HuAst_pathway_filtered, by = join_by(Term)) %>% 
  full_join(SCZ_vs_HC_HuAst_pathway_filtered, by = join_by(Term)) %>% 
  full_join(SCZ_vs_nonSCZ_HuAst_pathway_filtered, by = join_by(Term)) %>% 
  full_join(SCZ_vs_HC_HuN_pathway_filtered, by = join_by(Term))

#conversion to -log10
logHeatmap <- Pathways %>%  
  mutate(across(all_of(pvalue_columns), ~ -log10(.), .names = "log10_{col}" )) %>%
  select(`log10_P-value.x`, `log10_P-value.y`, `log10_P-value.x.x`, `log10_P-value.y.y`, `log10_P-value.x.x.x`, `log10_P-value.y.y.y`) %>%
  mutate_all(~ ifelse(is.na(.), 0, .))
  
#heatmap
library(ComplexHeatmap)
library(circlize)
colorscale <- colorRamp2(c(min(logHeatmap), 6, max(logHeatmap)),
                         c("white", "red", "maroon"))
Heatmap(fulljoin, col = colorscale)  

csv_file = "PathwaysMerged.csv" 
write.csv(Pathways, file = csv_file, row.names = FALSE)

Pathway_edit <- read_csv("Book8.csv")
Pathway_bundled <- Pathway_edit %>%
  mutate(across(all_of(pvalue_columns), ~ -log10(.), .names = "log10_{col}" )) %>%
  select(`log10_P-value.x`, `log10_P-value.y`, `log10_P-value.x.x`, `log10_P-value.y.y`, `log10_P-value.x.x.x`, `log10_P-value.y.y.y`) %>%
  mutate_all(~ ifelse(is.na(.), 0, .))

Heatmap(Pathway_bundled, col = colorscale, 
        row_order = order(as.numeric(gsub("row", "", rownames(Pathway_bundled)))), 
        column_order = order(as.numeric(gsub("column", "", colnames(Pathway_bundled)))))
