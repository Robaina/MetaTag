# Metadata from Arctic analysis
# Rafa, 21.06.2022
##############################

rm(list = ls()) #remove all the elements from your working enviroment
getwd() #gives you where you are
dir <- 'D:/###CNB/Traits/Artico/R_analisis/Data' #set the directory as an object 'dir'
setwd(file.path(dir)) #set working directory using the object and the subfolder

# Load libraries:
library(ggplot2)
library(dplyr)
library("RColorBrewer")
library(data.table)
library(reshape2)
library(tidyverse)
library(zoo)
library(SQMtools)

##############################

#Import environmental data
#Use a reformated table where I remove two points of salinity: 0 m 20/06 and 2,5 m 

metadata <- read.csv("../data/CB_dataset_Jan18_reformatted.txt", header = TRUE, sep="\t")
metadata.m <- melt(metadata, id=1:3)
metadata.m <- na.omit(metadata.m)

#Plot

Snow_PAR_Si <- ggplot(subset(metadata.m, variable %in% c("snow_depth_.cm.", "Salinity_.psu_per_mil.", "std_Si_.uM.")), 
                   aes(x=Day_of_year, y=value, color=variable)) +
  geom_line (size=2) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(panel.grid= element_blank()) +
  theme(panel.background = element_rect (fill= 'white', color= 'white')) +
  theme(legend.text = element_text(face="bold.italic")) +
  facet_grid(depth_.m. ~.  , scales="free")
Snow_PAR_Si

ggsave("../Results/SnowPARSi.pdf", plot=Snow_PAR_Si,units=c("mm"), dpi=300, scale=2)

Nutrients <- ggplot(subset(metadata.m, variable %in% c("Nox.PO4", "std_PO4_.uM.", "std_Nox_.uM.", "Chl_a_.ug_l.1.")), 
                    aes(x=Day_of_year, y=value, color=variable)) +
  geom_line (size=2) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(panel.grid= element_blank()) +
  theme(panel.background = element_rect (fill= 'white', color= 'white')) +
  theme(legend.text = element_text(face="bold.italic")) +
  facet_grid(depth_.m. ~.)
Nutrients

ggsave("../Results/Nutrients.pdf", plot=Nutrients,units=c("mm"), dpi=300, scale=2)

PAR_POC_CN <- ggplot(subset(metadata.m, variable %in% c("X_PAR_.uE_m.2_s.1.", "POC_.ug_l.1.", "C.N")), 
                    aes(x=Day_of_year, y=value, color=variable)) +
  geom_line (size=2) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(panel.grid= element_blank()) +
  theme(panel.background = element_rect (fill= 'white', color= 'white')) +
  theme(legend.text = element_text(face="bold.italic")) +
  facet_grid(depth_.m. ~.)
PAR_POC_CN

ggsave("../Results/PAR_POC_CN.pdf", plot=PAR_POC_CN,units=c("mm"), dpi=300, scale=2)

POC_ChlA <- ggplot(subset(metadata.m, variable %in% c("POC.Chla")), 
                     aes(x=Day_of_year, y=value, color=variable)) +
  geom_line (size=2) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(panel.grid= element_blank()) +
  theme(panel.background = element_rect (fill= 'white', color= 'white')) +
  theme(legend.text = element_text(face="bold.italic")) +
  facet_grid(depth_.m. ~.)
POC_ChlA

ggsave("../Results/POC_Chla_ratio.pdf", plot=POC_ChlA,units=c("mm"), dpi=300, scale=2)

#Import genetic data

GTDB_placement <- read.csv("./all_arctic_genes_assignments_scoreover0_funcolumn.tsv", header = TRUE, sep="\t")
SQM_taxonomy <- read.csv("./13.artico022_DNA.orftable_subset", header= TRUE, sep="\t")
#drop columns not needed
SQM_taxonomy_TPM<- SQM_taxonomy[,c(1,9,21:44)]
SQM_taxonomy_TPM <- SQM_taxonomy_TPM %>% select(-(contains ("9_Mar_RNA")| contains ("23_Apr_RNA")|contains ("1_May_RNA") | contains ("19_May_RNA")))

#Merge table
master_tab <- merge(GTDB_placement, SQM_taxonomy_TPM, by.x="query_name", by.y="Query_ID")
master_tab.m <- melt(master_tab, id.var= c("query_name","query_id","LWR","cluster_id","cluster_score","cluster_taxopath","taxopath","gene_fun","Tax"))
#Separate taxonomies
master_tab.m <- separate(master_tab.m, Tax, into = c("Domain_SQM","Phylum_SQM","Class_SQM","Order_SQM","Family_SQM","Genus_SQM","Species_SQM"), sep=";")
master_tab.m <- separate(master_tab.m, taxopath, into = c("Domain_GTDB","Phylum_GTDB","Class_GTDB","Order_GTDB","Family_GTDB","Genus_GTDB","Species_GTDB"), sep=";")
#Fill taxonomy
master_tab.m <- data.frame(t(apply(master_tab.m, 1, na.locf)))
master_tab.m$value <- as.numeric(master_tab.m$value)


Class_percent <- ggplot(arrange(master_tab.m, Domain_GTDB, Phylum_GTDB, Class_GTDB), aes(x=variable, y=value, 
                fill=(factor(Class_GTDB, levels=unique(Class_GTDB))))) +
  geom_bar (color="black", stat='identity',positio="fill") +
  #scale_fill_manual("Taxa", values=figure_col_tax_up)+
  #scale_fill_brewer("Domain", palette="Set1") +
  scale_fill_manual("Taxa", values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(master_tab.m$Class_GTDB)))) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=5, face="bold", color="black")) +
  theme(axis.text.y = element_text(size=12, face="bold", color="black")) +
  scale_x_discrete(expand = c(0, 0.6)) + scale_y_continuous(expand = c(0, 0))+
  ylab("TPM relative abundance (%) \n") +
  xlab("\n") + theme(legend.position="none") +
  facet_wrap(~ gene_fun) +
  theme(axis.title.y = element_text(size=12, face="bold", color="black")) +
  theme(legend.text = element_text(size=12), legend.title =element_text(size=12))

Class_percent
ggsave("../results/Class_percent_no_legend.pdf", plot=Class_percent,units=c("mm"), dpi=300, scale=2)

Class_total <- ggplot(arrange(master_tab.m, Domain_GTDB, Phylum_GTDB, Class_GTDB), aes(x=variable, y=value, 
                                                                                         fill=(factor(Class_GTDB, levels=unique(Class_GTDB))))) +
  geom_bar (color="black", stat='identity') +
  #scale_fill_manual("Taxa", values=figure_col_tax_up)+
  #scale_fill_brewer("Domain", palette="Set1") +
  scale_fill_manual("Taxa", values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(master_tab.m$Class_GTDB)))) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=5, face="bold", color="black")) +
  theme(axis.text.y = element_text(size=12, face="bold", color="black")) +
  scale_x_discrete(expand = c(0, 0.6)) + scale_y_continuous(expand = c(0, 0))+
  ylab("TPM total \n") +
  xlab("\n") + theme(legend.position="none") +
  facet_wrap(~ gene_fun, scales="free") +
  theme(axis.title.y = element_text(size=12, face="bold", color="black")) +
  theme(legend.text = element_text(size=12), legend.title =element_text(size=12))

Class_total
ggsave("../results/Class_total_no_legend.pdf", plot=Class_total,units=c("mm"), dpi=300, scale=2)

### Obtaining the general overview of gene taxonomy

general_overview <- master_tab
general_overview <- separate(general_overview, Tax, into = c("Domain_SQM","Phylum_SQM","Class_SQM","Order_SQM","Family_SQM","Genus_SQM","Species_SQM"), sep=";")
general_overview <- separate(general_overview, taxopath, into = c("Domain_GTDB","Phylum_GTDB","Class_GTDB","Order_GTDB","Family_GTDB","Genus_GTDB","Species_GTDB"), sep=";")
general_overview <- data.frame(t(apply(general_overview, 1, na.locf)))
general_overview$presence <- 1

general_overview_Class_GTDB <-general_overview %>%
  group_by(Domain_GTDB, Phylum_GTDB, Class_GTDB, gene_fun) %>%
  summarise(value = sum(presence))

#get the most abundant
general_overview_Class_GTDB_20 <- general_overview_Class_GTDB %>% unite("Merged", Domain_GTDB:Class_GTDB, remove=FALSE, sep="___") #Merge all the character columns since topAbund needs only a matrix
general_overview_Class_GTDB_20 <- dcast(general_overview_Class_GTDB_20, Merged ~ gene_fun ,fill = 0, value.var = "value" ) # dcast the melt dataframe por plotabund
rownames(general_overview_Class_GTDB_20) <- general_overview_Class_GTDB_20$Merged #the matrix for plotAbund can only accept numbers, so info goes to row.names
general_overview_Class_GTDB_20 <- general_overview_Class_GTDB_20[,c(2:ncol(general_overview_Class_GTDB_20))] #remove the Merged column
general_overview_Class_GTDB_20 <- mostAbundant(general_overview_Class_GTDB_20, N=20,  others=TRUE) #get the top 15 (in many cases there will not be 15)

general_overview_Class_GTDB_20$Info <- row.names(general_overview_Class_GTDB_20)
row.names(general_overview_Class_GTDB_20) <- 1:nrow(general_overview_Class_GTDB_20)
general_overview_Class_GTDB_20 <- separate(general_overview_Class_GTDB_20,col="Info", into=c("Domain","Phylum","Class"),sep = "___")
general_overview_Class_GTDB_20 <- general_overview_Class_GTDB_20[,c(10:12,1:9)]
general_overview_Class_GTDB_20 <- data.frame(t(apply(general_overview_Class_GTDB_20, 1, na.locf)))
general_overview_Class_GTDB_20 <- melt(general_overview_Class_GTDB_20, id.vars = c("Domain","Phylum","Class"))
general_overview_Class_GTDB_20$value <- as.numeric(general_overview_Class_GTDB_20$value)
general_overview_Class_GTDB_20$variable <- as.character(general_overview_Class_GTDB_20$variable)

#Import the simplified color file
color_general <- read.csv("./color_code_general_overview_20.txt", sep="\t", header = TRUE)
color_general <- color_general[,c(3,6)]
general_overview_Class_GTDB_20 <- merge(general_overview_Class_GTDB_20, color_general, by="Class")

figure_col_general <- c("Archaea"="#E01853","Bacteria"="#F99392","Actinobacteriota"="#FCD84A","Aquificota"="#F6BFCA","Bacteroidota"="#CD672F","Chloroflexota"="#8BC395","Desulfuromonadota"="#89CB6C","Firmicutes"="#40A635","Gemmatimonadota"="#919D5F","Nitrospinota"="#A6CEE3","Alphaproteobacteria"="#1F74CD","Gammaproteobacteria"="#5D478B","Proteobacteria"="#BFA5CF","SAR324"="#E7E099","Verrucomicrobiae"="#8B461D","Other"="#bdbdbd","Unspecified"="#636363")

#Plot the general overview

Class_general_overview <- ggplot(arrange(general_overview_Class_GTDB_20, Domain, Phylum, Class), aes(x=variable, y=value, 
                                fill=(factor(Taxa, levels=unique(Taxa))))) +
  geom_bar (color="black", stat='identity',position="fill") +
  scale_fill_manual("Taxa", values=figure_col_general)+
  #scale_fill_brewer("Domain", palette="Set1") +
  #scale_fill_manual("Taxa", values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(general_overview_Class_GTDB_20$Taxa)))) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=12, face="bold", color="black")) +
  theme(axis.text.y = element_text(size=12, face="bold", color="black")) +
  scale_x_discrete(expand = c(0, 0.6)) + scale_y_continuous(expand = c(0, 0))+
  ylab("Relative abundance (%) \n") + ggtitle("Phylogenetic distribution of placed queries") +
  xlab("\n") + #theme(legend.position="none") +
  theme(axis.title.y = element_text(size=12, face="bold", color="black")) +
  theme(legend.text = element_text(size=12), legend.title =element_text(size=12))

Class_general_overview
ggsave("../results/Class_general_overview_20.pdf", plot=Class_general_overview,units=c("mm"), dpi=300, scale=2)

Class_general_overview_total <- ggplot(arrange(general_overview_Class_GTDB_20, Domain, Phylum, Class), aes(x=variable, y=value, 
                                      fill=(factor(Taxa, levels=unique(Taxa))))) +
  geom_bar (color="black", stat='identity') +
  scale_fill_manual("Taxa", values=figure_col_general)+
  #scale_fill_brewer("Domain", palette="Set1") +
  #scale_fill_manual("Taxa", values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(general_overview_Class_GTDB_20$Taxa)))) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=12, face="bold", color="black")) +
  theme(axis.text.y = element_text(size=12, face="bold", color="black")) +
  scale_x_discrete(expand = c(0, 0.6)) + scale_y_continuous(expand = c(0, 0))+
  ylab("Total hits \n") + ggtitle("Phylogenetic distribution of placed queries") +
  xlab("\n") + theme(legend.position="none") +
  #facet_wrap( ~ variable, scales="free") +
  theme(axis.title.y = element_text(size=12, face="bold", color="black")) +
  theme(legend.text = element_text(size=12), legend.title =element_text(size=12))

Class_general_overview_total
ggsave("../results/Class_general_overview_total.pdf", plot=Class_general_overview_total,units=c("mm"), dpi=300, scale=2)


### Obtaining 15 top clasess

master_tab.m_Class_GTDB <-master_tab.m %>%
  group_by(Domain_GTDB, Phylum_GTDB, Class_GTDB, gene_fun, variable) %>%
  summarise(value = sum(value))

#Iterate to make small tables of genes and subset 15 top results in each table
genes <- c("amt","narG","nirK","amoA","nasA","nrtA","nosZ","hzsA","nirS")

genes_tables=list() # make an empty list

for (g in genes){ #fill the list with individual tables
  genes_tables[[g]] <- master_tab.m_Class_GTDB[ master_tab.m_Class_GTDB$gene_fun == g,] #create table with hits for that gene
  genes_tables[[g]] <- genes_tables[[g]] %>% unite("Merged", Domain_GTDB:gene_fun, remove=FALSE, sep="___") #Merge all the character columns since topAbund needs only a matrix
  genes_tables[[g]] <- dcast(genes_tables[[g]], Merged ~ variable ,value.var = "value" ) # dcast the melt dataframe por plotabund
  rownames(genes_tables[[g]]) <- getElement(genes_tables[[g]], "Merged") #the matrix for plotAbund can only accept numbers, so info goes to row.names
  genes_tables[[g]] <- subset(genes_tables[[g]], select=c(2:21)) #remove the Merged column
  genes_tables[[g]] <- mostAbundant(genes_tables[[g]], N=15,  others=TRUE) #get the top 15 (in many cases there will not be 15)
  rm(g)
}

#Combine dataframes

master_tab.m_Class_GTDB_n15.m <- do.call("rbind", genes_tables)

#restore the format
master_tab.m_Class_GTDB_n15.m$Info <- row.names(master_tab.m_Class_GTDB_n15.m)
row.names(master_tab.m_Class_GTDB_n15.m) <- 1:nrow(master_tab.m_Class_GTDB_n15.m)
master_tab.m_Class_GTDB_n15.m <- separate(master_tab.m_Class_GTDB_n15.m,col="Info", into=c("gene_fun","Merged"),sep = "\\.")
master_tab.m_Class_GTDB_n15.m <- separate(master_tab.m_Class_GTDB_n15.m,col="Merged", into=c("Domain","Phylum","Class","gene_fun2"),sep = "___")
master_tab.m_Class_GTDB_n15.m <- master_tab.m_Class_GTDB_n15.m[,c(21:24,1:20)]
master_tab.m_Class_GTDB_n15.m <- data.frame(t(apply(master_tab.m_Class_GTDB_n15.m, 1, na.locf)))
master_tab.m_Class_GTDB_n15.m <- melt(master_tab.m_Class_GTDB_n15.m, id.vars = c("Domain","Phylum","Class","gene_fun"))
master_tab.m_Class_GTDB_n15.m$value <- as.numeric(master_tab.m_Class_GTDB_n15.m$value)
master_tab.m_Class_GTDB_n15.m$variable <- as.character(master_tab.m_Class_GTDB_n15.m$variable)

#Import the simplified color file
color_sample <- read.csv("./color_code_sample_overview_15.txt", sep="\t", header = TRUE)
color_sample <- color_sample[,c(3,6)]
master_tab.m_Class_GTDB_n15.m <- merge(master_tab.m_Class_GTDB_n15.m, color_sample, by="Class")

figure_col_sample <- c("Archaea"="#E01853","Bacteria"="#F99392","Actinobacteriota"="#FCD84A","Bacteroidota"="#CD672F","Campylobacteriota"="#FBB268","Chloroflexota"="#8BC395","Desulfuromonadota"="#89CB6C","Firmicutes"="#40A635","Gemmatimonadota"="#919D5F","Marinisomatia"="#BFA5CF","Nitrospinota"="#A6CEE3","Nitrospiria"="#6CA9CF","Patescibacteria"="#f781bf","Planctomycetota"="#F6BFCA","Alphaproteobacteria"="#1F74CD","Gammaproteobacteria"="#5D478B","SAR324"="#E7E099","Verrucomicrobiae"="#8B461D","Other"="#bdbdbd","Unspecified"="#636363")


#Plot top 15

Class_percent15 <- ggplot(arrange(master_tab.m_Class_GTDB_n15.m, Domain, Phylum, Class), aes(x=variable, y=value, 
                          fill=(factor(Taxa, levels=unique(Taxa))))) +
  geom_bar (color="black", stat='identity',position="fill") +
  scale_fill_manual("Taxa", values=figure_col_sample)+
  #scale_fill_brewer("Domain", palette="Set1") +
  #scale_fill_manual("Taxa", values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(master_tab.m_Class_GTDB_n15.m$Class)))) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=5, face="bold", color="black")) +
  theme(axis.text.y = element_text(size=12, face="bold", color="black")) +
  scale_x_discrete(expand = c(0, 0.6)) + scale_y_continuous(expand = c(0, 0))+
  ylab("TPM relative abundance (%) \n") +
  xlab("\n") + #theme(legend.position="none") +
  facet_wrap(~ gene_fun) +
  theme(axis.title.y = element_text(size=12, face="bold", color="black")) +
  theme(legend.text = element_text(size=12), legend.title =element_text(size=12))

Class_percent15
ggsave("../results/Class_percent_no_legend.pdf", plot=Class_percent15,units=c("mm"), dpi=300, scale=2)

Class_total15 <- ggplot(arrange(master_tab.m_Class_GTDB_n15.m, Domain, Phylum, Class), aes(x=variable, y=value, 
                       fill=(factor(Class, levels=unique(Class))))) +
  geom_bar (color="black", stat='identity') +
  #scale_fill_manual("Taxa", values=figure_col_tax_up)+
  #scale_fill_brewer("Domain", palette="Set1") +
  scale_fill_manual("Taxa", values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(master_tab.m$Class)))) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=5, face="bold", color="black")) +
  theme(axis.text.y = element_text(size=12, face="bold", color="black")) +
  scale_x_discrete(expand = c(0, 0.6)) + scale_y_continuous(expand = c(0, 0))+
  ylab("TPM total \n") +
  xlab("\n") + theme(legend.position="none") +
  facet_wrap(~ gene_fun, scales="free") +
  theme(axis.title.y = element_text(size=12, face="bold", color="black")) +
  theme(legend.text = element_text(size=12), legend.title =element_text(size=12))

Class_total15
ggsave("../results/Class_total_no_legend.pdf", plot=Class_total15,units=c("mm"), dpi=300, scale=2)

