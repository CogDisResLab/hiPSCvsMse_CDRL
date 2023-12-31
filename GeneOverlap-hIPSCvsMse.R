#Jae overlap file
# install.packages("Geneoverlap")
# install.packages("BiocManager", repos = "https://cloud.r-project.org")
# BiocManager::install("Geneoverlap")
library(GeneOverlap)
# install.packages("dplyr")
library(dplyr)
# install.packages("ggplot2")
library(ggplot2)
library(readr)
#import csv's
Astrocytes_Mono_VS_control_GSE191248_ <- read_csv("Astrocytes  (Mono) VS control(GSE191248).csv")
SCZVSControl_CorticolNeu_GSE182875_ <- read_csv("SCZVSControl(CorticolNeu)(GSE182875).csv")
With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_ <- read_csv("With SCZ VS Without SCZ (Astrocytes)(GSE191248).csv")
Without_SCZ_and_Control_Astrocytes_GSE191248_ <- read_csv("Without SCZ and Control(Astrocytes)(GSE191248).csv")
Mouse_Astrocyte_E2_knockdown <- read_csv("Mouse_Astrocyte_E2_knockdown.csv")
Mouse_neuron_E2_knockdown <- read_csv("Mouse_neuron_E2_knockdown.csv")

#data wrangling the csv data sets and selecting the 200 highest values for upregulation and downregulation
  #selected -> gene name(symbol or hsapien_homoloh...), pvalue, and log2FoldChange, 
  #filter -> NAs and onyl significant values (pvalue <= 0.05)
  #data is arranged creating table in increasing order and taking the top 200 for the downregulation and bottom 200 for the upregulation
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

#making the axis of the heatmap
  #two lists one for all Human dataset and one for Mouse datasets
  #removed SCZvsnonSCZ because the gene smaple in 93-98 which is significantly less than the 200 of the others
Human <- list(hIPSC_N_SCZ_1_up = upreg_SCZVSControl_CorticolNeu_GSE182875_$Symbol,
              hIPSC_N_SCZ_1_down = downreg_SCZVSControl_CorticolNeu_GSE182875_$Symbol,
              hIPSC_Ast_SCZ_1_up = upreg_Astrocytes_Mono_VS_control_GSE191248_$Symbol,
              hIPSC_Ast_SCZ_1_down = downreg_Astrocytes_Mono_VS_control_GSE191248_$Symbol,
              hIPSC_Ast_SCZ_2_up = upreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_$Symbol,
              hIPSC_Ast_SCZ_2_down = downreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_$Symbol,
              hIPSC_Ast_SCZ_3_up = upreg_Without_SCZ_and_Control_Astrocytes_GSE191248_$Symbol, 
              hIPSC_Ast_SCZ_3_down = downreg_Without_SCZ_and_Control_Astrocytes_GSE191248_$Symbol)
Mouse <- list(mse_N_E2KO_up = upreg_Mouse_neuron_E2_knockdown$hsapiens_homolog_associated_gene_name,
              mse_N_E2KO_down = downreg_Mouse_neuron_E2_knockdown$hsapiens_homolog_associated_gene_name, 
              mse_Ast_E2KO_up = upreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name,
              mse_Ast_E2KO_down = downreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name)

background_size <- length(Reduce(union, c(Mouse, Human)))

#creating gene overlap
gom.obj <- newGOM(Human, Mouse, genome.size = background_size)
gom.obj
#creating heatmap
# dev.off()
par(mar = c(1, 1, 1, 1))
drawHeatmap(gom.obj, 
            what="odds.ratio",
            ncolused=9, grid.col="Oranges", note.col="black")

#All Mouse and Human gene names
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
