#Jae overlap file
install.packages("Geneoverlap")
install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install("Geneoverlap")
library("GeneOverlap")
install.packages("dplyr")
library(dplyr)
install.packages("ggplot2")
library(ggplot2)

#import csv's

#data wrangling the csv data sets and selecting the 200 highest values for upregulation and downregulation
  #selected -> gene name(symbol or hsapien_homoloh...), pvalue, and log2FoldChange, 
  #filter -> NAs and onyl significant values (pvalue <= 0.05)
  #data is arranged creating table in increasing order and taking the top 200 for the downregulation and bottom 200 for the upregulation
#HUMAN
    #astrocyte
upreg_Astrocytes_Mono_VS_control_GSE191248_ <- 
  select(Astrocytes_Mono_VS_control_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05) %>% 
  arrange(log2FoldChange) %>%
  tail(200) 
upreg_Astrocytes_Mono_VS_control_GSE191248_
downreg_Astrocytes_Mono_VS_control_GSE191248_ <- 
  select(Astrocytes_Mono_VS_control_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05) %>% 
  arrange(log2FoldChange) %>%
  head(200) 
downreg_Astrocytes_Mono_VS_control_GSE191248_
upreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_ <- 
  select(With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05) %>% 
  arrange(log2FoldChange) %>%
  tail(200) 
upreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_
downreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_ <- 
  select(With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05) %>% 
  arrange(log2FoldChange) %>%
  head(200) 
downreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_
upreg_Without_SCZ_and_Control_Astrocytes_GSE191248_ <- 
  select(Without_SCZ_and_Control_Astrocytes_GSE191248_ , Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05) %>% 
  arrange(log2FoldChange) %>%
  tail(200) 
upreg_Without_SCZ_and_Control_Astrocytes_GSE191248_ 
downreg_Without_SCZ_and_Control_Astrocytes_GSE191248_  <- 
  select(Without_SCZ_and_Control_Astrocytes_GSE191248_ , Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05) %>% 
  arrange(log2FoldChange) %>%
  head(200) 
downreg_Without_SCZ_and_Control_Astrocytes_GSE191248_ 
    #neuron
upreg_SCZVSControl_CorticolNeu_GSE182875_ <- 
  select(SCZVSControl_CorticolNeu_GSE182875_, Symbol, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05) %>% 
  arrange(log2FoldChange) %>%
  tail(200) 
upreg_SCZVSControl_CorticolNeu_GSE182875_
downreg_SCZVSControl_CorticolNeu_GSE182875_ <- 
  select(SCZVSControl_CorticolNeu_GSE182875_, Symbol, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05) %>% 
  arrange(log2FoldChange) %>%
  head(200) 
downreg_SCZVSControl_CorticolNeu_GSE182875_
#MOUSE
    #astrocyte
upreg_Mouse_Astrocyte_E2_knockdown <- 
  select(Mouse_Astrocyte_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, !is.na(hsapiens_homolog_associated_gene_name)) %>% 
  arrange(log2FoldChange) %>%
  tail(200) 
upreg_Mouse_Astrocyte_E2_knockdown
downreg_Mouse_Astrocyte_E2_knockdown <-
  select(Mouse_Astrocyte_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, !is.na(hsapiens_homolog_associated_gene_name)) %>% 
  arrange(log2FoldChange) %>%
  head(200) 
downreg_Mouse_Astrocyte_E2_knockdown
    #neuron
upreg_Mouse_neuron_E2_knockdown <- 
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, !is.na(hsapiens_homolog_associated_gene_name)) %>% 
  arrange(log2FoldChange) %>%
  tail(200) 
upreg_Mouse_neuron_E2_knockdown
downreg_Mouse_neuron_E2_knockdown <-
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, !is.na(hsapiens_homolog_associated_gene_name)) %>% 
  arrange(log2FoldChange) %>%
  head(200) 
downreg_Mouse_neuron_E2_knockdown

#making the axis of the heatmap
  #two lists one for all Human dataset and one for Mouse datasets
  #removed SCZvsnonSCZ because the gene smaple in 93-98 which is significantly less than the 200 of the others
Human <- list(mse_N_SCZ_1_up = upreg_SCZVSControl_CorticolNeu_GSE182875_$Symbol,
              mse_N_SCZ_1_down = downreg_SCZVSControl_CorticolNeu_GSE182875_$Symbol,
              hIPSC_Ast_SCZ_1_up = upreg_Astrocytes_Mono_VS_control_GSE191248_$Symbol,
              hIPSC_Ast_SCZ_1_down = downreg_Astrocytes_Mono_VS_control_GSE191248_$Symbol, 
              hIPSC_Ast_SCZ_2_up = upreg_Without_SCZ_and_Control_Astrocytes_GSE191248_$Symbol, 
              hIPSC_Ast_SCZ_2_down = downreg_Without_SCZ_and_Control_Astrocytes_GSE191248_$Symbol)
Mouse <- list(mse_N_E2KO_up = upreg_Mouse_neuron_E2_knockdown$hsapiens_homolog_associated_gene_name,
              mse_N_E2KO_down = downreg_Mouse_neuron_E2_knockdown$hsapiens_homolog_associated_gene_name, 
              mse_AST_E2KO_up = upreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name,
              mse_AST_E2KO_down = downreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name)
#creating gene overlap
gom.obj <- newGOM(Human, Mouse)
gom.obj
#creating heatmap
dev.off()
par(mar = c(1, 1, 1, 1))
drawHeatmap(gom.obj, 
            what=c("odds.ratio", "Jaccard"), log.scale=F, adj.p=F, cutoff=.05, 
            ncolused=9, grid.col=c("Greens", "Blues", "Greys", 
                                   "Oranges", "Purples", "Reds"), note.col="red")

