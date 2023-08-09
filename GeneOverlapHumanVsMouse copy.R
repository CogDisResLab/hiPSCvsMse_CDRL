#Jae overlap file
install.packages("Geneoverlap")
install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install("Geneoverlap")
library("GeneOverlap")
library(dplyr)

#Basic Gene Overlap
#Hu Data sets (symbol, pvalue, log2foldchange) pvalue <= 0.05
GO_Astrocytes_Mono_VS_control_GSE191248_ <- 
  select(Astrocytes_Mono_VS_control_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05)
GO_Astrocytes_Mono_VS_control_GSE191248_
GO_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_ <-
  select(With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05)
GO_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_
GO_Without_SCZ_and_Control_Astrocytes_GSE191248_ <- 
  select(Without_SCZ_and_Control_Astrocytes_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)) ,pvalue <= 0.05)
GO_Without_SCZ_and_Control_Astrocytes_GSE191248_
GO_SCZVSControl_CorticolNeu_GSE182875_ <- 
  select(SCZVSControl_CorticolNeu_GSE182875_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05)
GO_Without_SCZ_and_Control_Astrocytes_GSE191248_
#Mu Data sets (homolog name, pvalue, log2foldchange) pvalue <= 0.05
GO_Mouse_Astrocyte_E2_knockdown <- 
  select(Mouse_Astrocyte_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter((!is.na(hsapiens_homolog_associated_gene_name)), (pvalue <= 0.05))
GO_Mouse_Astrocyte_E2_knockdown
GO_Mouse_neuron_E2_knockdown <- 
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter((!is.na(hsapiens_homolog_associated_gene_name)), (pvalue <= 0.05))
GO_Mouse_neuron_E2_knockdown
#list for axis of heatmap
Hu<- list(HumanMonovsControlAstrocyte = GO_Astrocytes_Mono_VS_control_GSE191248_$Symbol,
          HumanSCZvsWithoutSCZ = GO_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_$Symbol,
          HumanWithoutSCZcsControlAstrocyte = GO_Without_SCZ_and_Control_Astrocytes_GSE191248_$Symbol,
          HumanSCZvsControlNeuron = GO_SCZVSControl_CorticolNeu_GSE182875_$Symbol)
Mu <- list(MouseAstrocyte = GO_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name,
           MouseNeuron = GO_Mouse_neuron_E2_knockdown$hsapiens_homolog_associated_gene_name)

D1 <- list(HumanMonovsControlAstrocyte = GO_Astrocytes_Mono_VS_control_GSE191248_$Symbol,
           HumanSCZvsWithoutSCZ = GO_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_$Symbol,
           HumanWithoutSCZcsControlAstrocyte = GO_Without_SCZ_and_Control_Astrocytes_GSE191248_$Symbol,
           HumanSCZvsControlNeuron = GO_SCZVSControl_CorticolNeu_GSE182875_$Symbol, 
           MouseAstrocyte = GO_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name,
           MouseNeuron = GO_Mouse_neuron_E2_knockdown$hsapiens_homolog_associated_gene_name)

D2 <- list(MouseAstrocyte = GO_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name,
               MouseNeuron = GO_Mouse_neuron_E2_knockdown$hsapiens_homolog_associated_gene_name, 
               HumanMonovsControlAstrocyte = GO_Astrocytes_Mono_VS_control_GSE191248_$Symbol,
               HumanSCZvsWithoutSCZ = GO_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_$Symbol,
               HumanWithoutSCZcsControlAstrocyte = GO_Without_SCZ_and_Control_Astrocytes_GSE191248_$Symbol,
               HumanSCZvsControlNeuron = GO_SCZVSControl_CorticolNeu_GSE182875_$Symbol)

#Geneoverlap obj
library(GeneOverlap)
install.packages("ggplot2")
library(ggplot2)
data(GeneOverlap)
gom.obj <- newGOM(Hu, Mu)
gom.obj
go.obj1 <- testGeneOverlap(gom.obj)
drawHeatmap(gom.obj)
par(mar = c(1, 1, 1, 1))
plot(1:10)
drawHeatmap(gom.obj, 
            what=c("odds.ratio", "Jaccard"), log.scale=F, adj.p=F, cutoff=.05, 
            ncolused=9, grid.col=c("Greens", "Blues", "Greys", 
                                   "Oranges", "Purples", "Reds"), note.col="red")
dev.off()

#GeneOverlap all data sets on each axis DOES NOT WORK
go.obj.all <- newGOM(D1, D2)
go.obj.all
go.obj2 <- testGeneOverlap(go.obj.all)
drawHeatmap(go.obj.all)
par(mar = c(1, 1, 1, 1))
plot(1:100)
drawHeatmap(go.obj.all, 
            what=c("odds.ratio", "Jaccard"), log.scale=F, adj.p=F, cutoff=.05, 
            ncolused=9, grid.col=c("Greens", "Blues", "Greys", 
                                   "Oranges", "Purples", "Reds"), note.col="red")


#Example
Exgom.obj <- newGOM(hESC.ChIPSeq.list, hESC.RNASeq.list,
                  + gs.RNASeq)
Exgom.obj
drawHeatmap(Exgom.obj)
par(mar = c(1, 1, 1, 1))
plot(1:30)
drawHeatmap(Exgom.obj, 
            what=c("odds.ratio", "Jaccard"), log.scale=F, adj.p=F, cutoff=.05, 
            ncolused=9, grid.col=c("Greens", "Blues", "Greys", 
                                   "Oranges", "Purples", "Reds"), note.col="red")


#Overlap with up and downregulation using log2FoldChange
#regulated Hu data sets NEURON (Symbol, pvalue, log2FoldChange)
#PVALUE <= 0.05 and log2foldchange >< 1
upreg_SCZVSControl_CorticolNeu_GSE182875_ <- 
  select(SCZVSControl_CorticolNeu_GSE182875_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05, log2FoldChange > 0)
upreg_SCZVSControl_CorticolNeu_GSE182875_
downreg_SCZVSControl_CorticolNeu_GSE182875_ <- 
  select(SCZVSControl_CorticolNeu_GSE182875_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05, log2FoldChange < 0)
downreg_SCZVSControl_CorticolNeu_GSE182875_
#regulated Mu data sets NEURON (homolog name, pvalue, log2FoldChange)
#PVALUE <= 0.05 and log2foldchange +- 0.26
upreg_Mouse_neuron_E2_knockdown <- 
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, log2FoldChange > 0)
upreg_Mouse_neuron_E2_knockdown
downreg_Mouse_neuron_E2_knockdown <-
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, log2FoldChange < 0)
downreg_Mouse_neuron_E2_knockdown

#list for axis of heat map
regulatedHu_neuron <- list(upHu = upreg_SCZVSControl_CorticolNeu_GSE182875_$Symbol,
                           downHu = downreg_SCZVSControl_CorticolNeu_GSE182875_$Symbol)
regulatedMu_neuron <- list(upMu = upreg_Mouse_neuron_E2_knockdown$hsapiens_homolog_associated_gene_name,
                           downMu = downreg_Mouse_neuron_E2_knockdown$hsapiens_homolog_associated_gene_name)
#NEURON HEATMAP
regNeurongom.obj <- newGOM(regulatedHu_neuron, regulatedMu_neuron)
regNeurongom.obj
par(mar = c(1, 1, 1, 1))
plot(1:30)
drawHeatmap(regNeurongom.obj, 
            what=c("odds.ratio", "Jaccard"), log.scale=F, adj.p=F, cutoff=.05, 
            ncolused=9, grid.col=c("Greens", "Blues", "Greys", 
                                   "Oranges", "Purples", "Reds"), note.col="red")
dev.off()

#regulated Hu data sets ASTROMONOvsHC (Symbol, pvalue, log2FoldChange)
#PVALUE <= 0.05 and log2foldchange >< 0
upreg_Astrocytes_Mono_VS_control_GSE191248_ <- 
  select(Astrocytes_Mono_VS_control_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05, log2FoldChange > 0)
upreg_Astrocytes_Mono_VS_control_GSE191248_
downreg_Astrocytes_Mono_VS_control_GSE191248_ <- 
  select(Astrocytes_Mono_VS_control_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05, log2FoldChange < 0)
downreg_Astrocytes_Mono_VS_control_GSE191248_
#regulated Mu data sets ASTROMONOvsHC (homolog name, pvalue, log2FoldChange)
#PVALUE <= 0.05 and log2foldchange >< 0
Aupreg_Mouse_Astrocyte_E2_knockdown <- 
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, log2FoldChange > 0)
Aupreg_Mouse_Astrocyte_E2_knockdown
Adownreg_Mouse_Astrocyte_E2_knockdown <-
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, log2FoldChange < 0)
Adownreg_Mouse_Astrocyte_E2_knockdown

#list for axis of heat map
regulatedHu_AstroMonovsHC <- list(upHu = upreg_Astrocytes_Mono_VS_control_GSE191248_$Symbol,
                           downHu = downreg_Astrocytes_Mono_VS_control_GSE191248_$Symbol)
AregulatedMu_Astro <- list(upMu = Aupreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name,
                           downMu = Adownreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name)
#ASTROMONOcsHC HEATMAP
regAstroMonovsHCgom.obj <- newGOM(regulatedHu_AstroMonovsHC, AregulatedMu_Astro)
regAstroMonovsHCgom.obj
par(mar = c(1, 1, 1, 1))
plot(1:30)
drawHeatmap(regAstroMonovsHCgom.obj, 
            what=c("odds.ratio", "Jaccard"), log.scale=F, adj.p=F, cutoff=.05, 
            ncolused=9, grid.col=c("Greens", "Blues", "Greys", 
                                   "Oranges", "Purples", "Reds"), note.col="red")
dev.off()

#regulated Hu data sets wSCZvswoutSCZ (Symbol, pvalue, log2FoldChange)
#PVALUE <= 0.05 and log2foldchange >< 0
upreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_ <- 
  select(With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05, log2FoldChange > 0)
upreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_
downreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_ <- 
  select(With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05, log2FoldChange < 0)
downreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_
#regulated Mu data sets wSCZvswoutSCZ (homolog name, pvalue, log2FoldChange)
#PVALUE <= 0.05 and log2foldchange >< 0
Bupreg_Mouse_Astrocyte_E2_knockdown <- 
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, log2FoldChange > 0)
Bupreg_Mouse_Astrocyte_E2_knockdown
Bdownreg_Mouse_Astrocyte_E2_knockdown <-
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, log2FoldChange < 0)
Bdownreg_Mouse_Astrocyte_E2_knockdown

#list for axis of heat map
regulatedHu_wSCZ_woutSCZ <- list(upHu = upreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_$Symbol,
                                  downHu = downreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_$Symbol)
BregulatedMu_Astro <- list(upMu = Bupreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name,
                          downMu = Bdownreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name)
#wSCZvswoutSCZ HEATMAP
regAstrowSCZ_woutSCZgom.obj <- newGOM(regulatedHu_wSCZ_woutSCZ, BregulatedMu_Astro)
regAstrowSCZ_woutSCZgom.obj
par(mar = c(1, 1, 1, 1))
plot(1:30)
drawHeatmap(regAstrowSCZ_woutSCZgom.obj, 
            what=c("odds.ratio", "Jaccard"), log.scale=F, adj.p=F, cutoff=.05, 
            ncolused=9, grid.col=c("Greens", "Blues", "Greys", 
                                   "Oranges", "Purples", "Reds"), note.col="red")
dev.off()

#regulated Hu data sets wSCZvsHC (Symbol, pvalue, log2FoldChange)
#PVALUE <= 0.05 and log2foldchange >< +- 0.26
upreg_Without_SCZ_and_Control_Astrocytes_GSE191248_ <- 
  select(Without_SCZ_and_Control_Astrocytes_GSE191248_ , Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05, log2FoldChange > 0.26)
upreg_Without_SCZ_and_Control_Astrocytes_GSE191248_ 
downreg_Without_SCZ_and_Control_Astrocytes_GSE191248_  <- 
  select(Without_SCZ_and_Control_Astrocytes_GSE191248_ , Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05, log2FoldChange < -0.26)
downreg_Without_SCZ_and_Control_Astrocytes_GSE191248_ 
#regulated Mu data sets wSCZvsHC (homolog name, pvalue, log2FoldChange)
#PVALUE <= 0.05 and log2foldchange >< +-0.26
Cupreg_Mouse_Astrocyte_E2_knockdown <- 
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, log2FoldChange > 0.26)
Cupreg_Mouse_Astrocyte_E2_knockdown
Cdownreg_Mouse_Astrocyte_E2_knockdown <-
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, log2FoldChange < -0.26)
Cdownreg_Mouse_Astrocyte_E2_knockdown

#list for axis of heat map
regulatedHu_wSCZ_HC <- list(upHu = upreg_SCZVSControl_CorticolNeu_GSE182875_$Symbol,
                                 downHu = downreg_SCZVSControl_CorticolNeu_GSE182875_$Symbol)
CregulatedMu_Astro <- list(upMu = Cupreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name,
                           downMu = Cdownreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name)
#wSCZvsHC HEATMAP
regAstrowSCZ_HCgom.obj <- newGOM(regulatedHu_wSCZ_HC, CregulatedMu_Astro)
regAstrowSCZ_HCgom.obj
par(mar = c(1, 1, 1, 1))
plot(1:30)
drawHeatmap(regAstrowSCZ_woutSCZgom.obj, 
            what=c("odds.ratio", "Jaccard"), log.scale=F, adj.p=F, cutoff=.05, 
            ncolused=9, grid.col=c("Greens", "Blues", "Greys", 
                                   "Oranges", "Purples", "Reds"), note.col="red")
dev.off()

#mixed Astrocyte Heatmaps
#MOUSE
upreg_Mouse_Astrocyte_E2_knockdown <- 
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, log2FoldChange > 0)
upreg_Mouse_Astrocyte_E2_knockdown
downreg_Mouse_Astrocyte_E2_knockdown <-
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, log2FoldChange < 0)
downreg_Mouse_Astrocyte_E2_knockdown
#HUMAN
allupreg_Astrocytes_Mono_VS_control_GSE191248_ <- 
  select(Astrocytes_Mono_VS_control_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05, log2FoldChange > 1.5)
allupreg_Astrocytes_Mono_VS_control_GSE191248_
alldownreg_Astrocytes_Mono_VS_control_GSE191248_ <- 
  select(Astrocytes_Mono_VS_control_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05, log2FoldChange < -1.5)
alldownreg_Astrocytes_Mono_VS_control_GSE191248_
allupreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_ <- 
  select(With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05, log2FoldChange > 1.5)
allupreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_
alldownreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_ <- 
  select(With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05, log2FoldChange < -1.5)
alldownreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_
allupreg_Without_SCZ_and_Control_Astrocytes_GSE191248_ <- 
  select(Without_SCZ_and_Control_Astrocytes_GSE191248_ , Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05, log2FoldChange > 1.5)
allupreg_Without_SCZ_and_Control_Astrocytes_GSE191248_ 
alldownreg_Without_SCZ_and_Control_Astrocytes_GSE191248_  <- 
  select(Without_SCZ_and_Control_Astrocytes_GSE191248_ , Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05, log2FoldChange < -1.5)
alldownreg_Without_SCZ_and_Control_Astrocytes_GSE191248_ 
#axis of heatmap
regulatedHu_Astro <- list(upMonovsHC = allupreg_Astrocytes_Mono_VS_control_GSE191248_$Symbol,
                          downMonovsHC = alldownreg_Astrocytes_Mono_VS_control_GSE191248_$Symbol, 
                          upSCZvsnoSCZ = allupreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_$Symbol,
                          downSCZvsnoSCZ = alldownreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_$Symbol,
                          upnoSCZvsHC= allupreg_Without_SCZ_and_Control_Astrocytes_GSE191248_$Symbol, 
                          downnoSCZvsHC= alldownreg_Without_SCZ_and_Control_Astrocytes_GSE191248_$Symbol)
regulatedMu_Astro <- list(upMu = upreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name,
                           downMu = downreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name)

#Heatmap all human astro data sets against mouse astro
allastrogom.obj <- newGOM(regulatedHu_Astro, regulatedMu_Astro)
allastrogom.obj
par(mar = c(1, 1, 1, 1))
plot(1:30)
drawHeatmap(allastrogom.obj, 
            what=c("odds.ratio", "Jaccard"), log.scale=F, adj.p=F, cutoff=.05, 
            ncolused=9, grid.col=c("Greens", "Blues", "Greys", 
                                   "Oranges", "Purples", "Reds"), note.col="red")

#Top 500
#ASTROCYTES
    #Hu
TOPupreg_Astrocytes_Mono_VS_control_GSE191248_ <- 
  select(Astrocytes_Mono_VS_control_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05) %>% 
  arrange(desc(log2FoldChange)) %>%
  head(200) 
TOPupreg_Astrocytes_Mono_VS_control_GSE191248_
TOPdownreg_Astrocytes_Mono_VS_control_GSE191248_ <- 
  select(Astrocytes_Mono_VS_control_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05) %>% 
  arrange(log2FoldChange) %>%
  head(200) 
TOPdownreg_Astrocytes_Mono_VS_control_GSE191248_
TOPupreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_ <- 
  select(With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05) %>% 
  arrange(desc(log2FoldChange)) %>%
  head(200) 
TOPupreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_
TOPdownreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_ <- 
  select(With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05) %>% 
  arrange(log2FoldChange) %>%
  head(200) 
TOPdownreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_
TOPupreg_Without_SCZ_and_Control_Astrocytes_GSE191248_ <- 
  select(Without_SCZ_and_Control_Astrocytes_GSE191248_ , Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05) %>% 
  arrange(desc(log2FoldChange)) %>%
  head(200) 
TOPupreg_Without_SCZ_and_Control_Astrocytes_GSE191248_ 
TOPdownreg_Without_SCZ_and_Control_Astrocytes_GSE191248_  <- 
  select(Without_SCZ_and_Control_Astrocytes_GSE191248_ , Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05) %>% 
  arrange(log2FoldChange) %>%
  head(200) 
TOPdownreg_Without_SCZ_and_Control_Astrocytes_GSE191248_ 
    #Mu
TOPupreg_Mouse_Astrocyte_E2_knockdown <- 
  select(Mouse_Astrocyte_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, !is.na(hsapiens_homolog_associated_gene_name)) %>% 
  arrange(desc(log2FoldChange)) %>%
  head(200) 
TOPupreg_Mouse_Astrocyte_E2_knockdown
TOPdownreg_Mouse_Astrocyte_E2_knockdown <-
  select(Mouse_Astrocyte_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, !is.na(hsapiens_homolog_associated_gene_name)) %>% 
  arrange(log2FoldChange) %>%
  head(200) 
TOPdownreg_Mouse_Astrocyte_E2_knockdown

TOPregulatedHu_Astro <- list(upMonovsHC = TOPupreg_Astrocytes_Mono_VS_control_GSE191248_$Symbol,
                          downMonovsHC = TOPdownreg_Astrocytes_Mono_VS_control_GSE191248_$Symbol, 
                          upSCZvsnoSCZ = TOPupreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_$Symbol,
                          downSCZvsnoSCZ = TOPdownreg_With_SCZ_VS_Without_SCZ_Astrocytes_GSE191248_$Symbol,
                          upnoSCZvsHC= TOPupreg_Without_SCZ_and_Control_Astrocytes_GSE191248_$Symbol, 
                          downnoSCZvsHC= TOPdownreg_Without_SCZ_and_Control_Astrocytes_GSE191248_$Symbol)
TOPregulatedMu_Astro <- list(upMu = TOPupreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name,
                          downMu = TOPdownreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name)

#Heatmap all human astro data sets against mouse astro
TOPastrogom.obj <- newGOM(TOPregulatedHu_Astro, TOPregulatedMu_Astro)
TOPastrogom.obj
par(mar = c(1, 1, 1, 1))
plot(1:30)
drawHeatmap(TOPastrogom.obj, 
            what=c("odds.ratio", "Jaccard"), log.scale=F, adj.p=F, cutoff=.05, 
            ncolused=9, grid.col=c("Greens", "Blues", "Greys", 
                                   "Oranges", "Purples", "Reds"), note.col="red")
#NEURON
    #Hu
TOPupreg_SCZVSControl_CorticolNeu_GSE182875_ <- 
  select(SCZVSControl_CorticolNeu_GSE182875_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05) %>% 
  arrange(desc(log2FoldChange)) %>%
  head(200) 
TOPupreg_SCZVSControl_CorticolNeu_GSE182875_
TOPdownreg_SCZVSControl_CorticolNeu_GSE182875_ <- 
  select(SCZVSControl_CorticolNeu_GSE182875_, Symbol, pvalue, log2FoldChange) %>%
  filter((!is.na(Symbol)), pvalue <= 0.05) %>% 
  arrange(log2FoldChange) %>%
  head(200) 
TOPdownreg_SCZVSControl_CorticolNeu_GSE182875_
    #Mu
TOPupreg_Mouse_neuron_E2_knockdown <- 
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, !is.na(hsapiens_homolog_associated_gene_name)) %>% 
  arrange(desc(log2FoldChange)) %>%
  head(200) 
TOPupreg_Mouse_neuron_E2_knockdown
TOPdownreg_Mouse_neuron_E2_knockdown <-
  select(Mouse_neuron_E2_knockdown, hsapiens_homolog_associated_gene_name, pvalue, log2FoldChange) %>%
  filter(pvalue <= 0.05, !is.na(hsapiens_homolog_associated_gene_name)) %>% 
  arrange(log2FoldChange) %>%
  head(200) 
TOPdownreg_Mouse_neuron_E2_knockdown


#list for axis of heat map
TOPregulatedHu_neuron <- list(upHu = TOPupreg_SCZVSControl_CorticolNeu_GSE182875_$Symbol,
                           downHu = TOPdownreg_SCZVSControl_CorticolNeu_GSE182875_$Symbol)
TOPregulatedMu_neuron <- list(upMu = TOPupreg_Mouse_neuron_E2_knockdown$hsapiens_homolog_associated_gene_name,
                           downMu = TOPdownreg_Mouse_neuron_E2_knockdown$hsapiens_homolog_associated_gene_name)

#NEURON HEATMAP
TOPregNeurongom.obj <- newGOM(TOPregulatedHu_neuron, TOPregulatedMu_neuron)
TOPregNeurongom.obj
par(mar = c(1, 1, 1, 1))
plot(1:30)
drawHeatmap(TOPregNeurongom.obj, 
            what=c("odds.ratio", "Jaccard"), log.scale=F, adj.p=F, cutoff=.05, 
            ncolused=9, grid.col=c("Greens", "Blues", "Greys", 
                                   "Oranges", "Purples", "Reds"), note.col="red")
dev.off()

Top200upregulatedMseAstrocyteE2KO <- column_matrix <- matrix(TOPupreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name, nrow = nrow(TOPupreg_Mouse_Astrocyte_E2_knockdown), ncol = 1)
Top200upregulatedMseAstrocyteE2KO
write.csv(Top200upregulatedMseAstrocyteE2KO, file = "Top200upregulatedMseAstrocyteE2KO.csv", row.names = FALSE)

#Final Heatmap
#
Human <- list(mse_N_SCZ_1_up = TOPupreg_SCZVSControl_CorticolNeu_GSE182875_$Symbol,
              mse_N_SCZ_1_down = TOPdownreg_SCZVSControl_CorticolNeu_GSE182875_$Symbol,
              hIPSC_Ast_SCZ_1_up = TOPupreg_Astrocytes_Mono_VS_control_GSE191248_$Symbol,
              hIPSC_Ast_SCZ_1_down = TOPdownreg_Astrocytes_Mono_VS_control_GSE191248_$Symbol, 
              hIPSC_Ast_SCZ_2_up = TOPupreg_Without_SCZ_and_Control_Astrocytes_GSE191248_$Symbol, 
              hIPSC_Ast_SCZ_2_down = TOPdownreg_Without_SCZ_and_Control_Astrocytes_GSE191248_$Symbol)
Mouse <- list(mse_N_E2KO_up = TOPupreg_Mouse_neuron_E2_knockdown$hsapiens_homolog_associated_gene_name,
              mse_N_E2KO_down = TOPdownreg_Mouse_neuron_E2_knockdown$hsapiens_homolog_associated_gene_name, 
              mse_AST_E2KO_up = TOPupreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name,
              mse_AST_E2KO_down = TOPdownreg_Mouse_Astrocyte_E2_knockdown$hsapiens_homolog_associated_gene_name)
Finalgom.obj <- newGOM(Human, Mouse)
Finalgom.obj
par(mar = c(1, 1, 1, 1))
drawHeatmap(Finalgom.obj, 
            what=c("odds.ratio", "Jaccard"), log.scale=F, adj.p=F, cutoff=.05, 
            ncolused=9, grid.col=c("Greens", "Blues", "Greys", 
                                   "Oranges", "Purples", "Reds"), note.col="red")
dev.off()