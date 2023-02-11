
# library(Hmisc)

# correlations
# prepare biological data
# subset gen_function (col8) and your TMPs per sample
master_tab_corr <- master_tab[,c(8,10:ncol(master_tab))]
# agrregate data by summing per gene on each sample
master_tab_corr_gen <- aggregate(x = master_tab_corr[names(master_tab_corr[,-1])],
                            by = master_tab_corr[c("gene_fun")],
                            FUN = sum)
row.names(master_tab_corr_gen) <- master_tab_corr_gen$gene_fun
master_tab_corr_gen.t <- as.data.frame(t(master_tab_corr_gen[2:ncol(master_tab_corr_gen)]))
master_tab_corr_gen.t$SAMPLE.NAME <- str_replace(row.names(master_tab_corr_gen.t), "TPM.", "")

# prepare environmental data
# subset metadata to sample name and continous variables only
# this depends on your data.frame
metadata_corr <- metadata[,c(1,6:ncol(metadata))]

# merge envir data and tpms per gene
master_corr <- merge(metadata_corr, master_tab_corr_gen.t, by.x="SAMPLE.NAME", by.y="SAMPLE.NAME")

row.names(master_corr) <- master_corr$SAMPLE.NAME
master_corr <- master_corr[-1]


# APUNTAD LAS COLUMNAS QUE CORRESPONDEN A VARIABLES AMBIENTALES Y LAS QUE SON LOS GENES
names(master_corr)


write.table(master_corr, file = "master_input_correlations.tsv", 
            quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

# APUNTAD LAS COLUMNAS QUE CORRESPONDEN A VARIABLES AMBIENTALES Y LAS QUE SON LOS GENES


# clean up

rm(list=ls())

# now run rhea_Correlations.r 