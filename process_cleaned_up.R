source("/mnt/raidtmp/Alejandro/functions.R")

library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)


library(SoupX)
# Load data and estimate soup profile
sc = load10X("/mnt/raidtmp/Alejandro/velocyto_florian/MKP3-DEMUX/MKP3/outs/")
# Estimate rho
sc = setContaminationFraction(sc, 0.1)
# Clean the data
out = adjustCounts(sc, roundToInt=TRUE)
saveRDS(out, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/soup_filtered_count_matrix_0.1.RDS")
#
out = readRDS("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/soup_filtered_count_matrix_0.1.RDS")


# this searches for all input matrices
#
files <- Sys.glob("/mnt/raidtmp/Alejandro/velocyto_florian/MKP3-DEMUX/MKP3/outs/raw_feature_bc_matrix/*.mtx.gz")
#files = grep(x=files, pattern="old", value=T, invert=T)
#files = grep(x=files, pattern="sepsis", value=T)

allfiles.raw = list()
allABs.raw = list()
for (file in files)
{
  samplename = "MKP3"
  #samplename = str_split(dirname(file), "/")[[1]][8] #estos numeros indican el nombre de la carpeta
  foldername = dirname(file)
  
  print(paste(samplename, foldername))
  
  h5file = Read10X(foldername,unique.features = TRUE)

  if (is.null(names(h5file)))
  {
      print(paste("WITHOUT AB", samplename))
    allfiles.raw[[samplename]] = h5file
  } else {
      print(paste("WITH AB", samplename))
    allfiles.raw[[samplename]] = h5file$`Gene Expression`
    allABs.raw[[samplename]] = h5file$`Antibody Capture`
  }

  print(paste(samplename, nrow(allfiles.raw[[samplename]]), "x", ncol(allfiles.raw[[samplename]]), "genes x cells"))
}
allfiles.raw$MKP3 <- out
joint.bcs <- intersect(colnames(allfiles.raw$MKP3), colnames(allABs.raw$MKP3))
obj.umis <- allfiles.raw$MKP3[, joint.bcs]
obj.htos <- as.matrix(allABs.raw$MKP3[, joint.bcs])

counts <- GetAssayData(seurat_obj, assay = "RNA")


# Confirm that the HTO have the correct names
rownames(obj.htos)


#
# here we create a list of seurat object. each entry corresponds to an input matrix from above
#

objlist = list()
for (x in names(allfiles.raw))
{

    matrix = allfiles.raw[[x]]
    
    # this creates a Seurat object from the count matrix. it sets the object's project to x and prepends the sample name to all cells
    # the patternlist.mouse contains patterns for mt and RP-genes
    filteredObj = makeSeuratObj(matrix, x, patternList.mouse)
    
    # this creates log-normalized count matrices in RNA assay
    filteredObj <- NormalizeData(filteredObj, verbose = FALSE)
    # this calculates the most (2000) variable features per data set. variable features are features which show a high variance between all cells of a sample
    filteredObj <- FindVariableFeatures(filteredObj, verbose = FALSE)
    
    objlist[[x]] = filteredObj

    print(x)
    print(filteredObj)
    
    
}

names(objlist)

rps_rownames = rownames(objlist$MKP3[grep('^Rps', rownames(objlist$MKP3))])
rpl_rownames = rownames(objlist$MKP3[grep('^Rpl', rownames(objlist$MKP3))])
objlist$MKP3 <- objlist$MKP3[!(row.names(objlist$MKP3) %in% rps_rownames),]
objlist$MKP3 <- objlist$MKP3[!(row.names(objlist$MKP3) %in% rpl_rownames),]


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

objlist.raw = objlist
objlist <- lapply(X = objlist.raw, FUN = function(obj) {
  # mt content: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6072887/
  print(paste("Seurat obj project", obj@project.name))
  print(obj)
  obj <- subset(obj, subset = nFeature_RNA > 250 & nFeature_RNA < 6000 & nCount_RNA > 400 & nCount_RNA < 20000)
  obj <- subset(obj, subset = percent.mt < 5)
  print(obj)
  
  return(obj)
})



for (name in names(objlist))
{
  p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  save_plot(p, paste("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq", paste(name, "filtered_violins_qc", sep="_"), sep="/"), fig.width=10, fig.height=6)

  p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, combine=F)
  p[[1]] = p[[1]] + scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000,100))
  p[[2]] = p[[2]] + scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000,100))
  p[[3]] = p[[3]] + scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,5))
  p = combine_plot_grid_list(plotlist=p, ncol=3)
  save_plot(p, paste("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq", paste(name, "filtered_violins_detail_qc", sep="_"), sep="/"), fig.width=18, fig.height=6)
  
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq", paste(name, "filtered_scatter_ncount_mt", sep="_"), sep="/"), fig.width=10, fig.height=6)
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.rp")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq", paste(name, "filtered_scatter_ncount_rp", sep="_"), sep="/"), fig.width=10, fig.height=6)
}

MKP3= objlist$MKP3

MKP3 <- FindVariableFeatures(MKP3, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(MKP3), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(MKP3)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(plot2, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/HVG_MKP3", fig.width=10, fig.height=6)


colnames(allABs.raw$MKP3)= paste("MKP3_", colnames(allABs.raw$MKP3), sep="")

obj.htos <- allABs.raw$MKP3[, colnames(MKP3)]
MKP3[["HTO"]] <- CreateAssayObject(counts = obj.htos)


MKP3 <- NormalizeData(MKP3, assay = "HTO", normalization.method = "CLR")
MKP3 <- HTODemux(MKP3, assay = "HTO", positive.quantile = 0.99)

table(MKP3$HTO_classification.global)


p= RidgePlot(MKP3, assay = "HTO", features = rownames(MKP3[["HTO"]][1:4]), ncol = 4)
save_plot(p, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/hto_classification", fig.width=40, fig.height=10)


#Singlet analysis
MKP3HTO=MKP3[,MKP3$HTO_classification.global!="Doublet"]

MKP3HTO=MKP3HTO[,MKP3HTO$HTO_classification.global!="Negative"]

cellList = colnames(MKP3HTO)

table(MKP3HTO$HTO_classification)

featVec <- vector(mode="character", length=length(cellList))
featVec = MKP3HTO$HTO_classification


featVec[featVec == "pDC_depletion"] = "pDC_depletion"
featVec[featVec == "pDC-Platelet_depletion"] = "pDC-Platelet_depletion"
featVec[featVec == "Platelet_depletion"] = "Platelet_depletion"
featVec[featVec == "Bl6"] = "Bl6"
#featVec[featVec == "Negative"] = "Negative"


MKP3HTO$CSclassification=featVec


###################################################### INTEGRATE SEURAT OBJECT HERE
MKP3HTO <- CellCycleScoring(MKP3HTO, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
MKP3HTO <- ScaleData(MKP3HTO, vars.to.regress = c("percent.mt", "percent.rp", "percent.rps", "percent.rpl","S.Score", "G2M.Score"))
MKP3HTO= RunPCA(MKP3HTO)
MKP3HTO= RunUMAP(MKP3HTO, dims = 1:30, reduction.key = "UMAP_")
MKP3HTO <- FindNeighbors(MKP3HTO, dims = 1:30)
MKP3HTO = readRDS("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/MKP3HTO.RDS")
MKP3HTO_new <- FindClusters(MKP3HTO_new, resolution = 0.3)


p=DimPlot(MKP3HTO, shuffle = T, seed = 1, group.by= "CSclassification", raster=FALSE)
save_plot(p, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/version9_wnn_ig_dimplot", 12, 8)
p=DimPlot(MKP3HTO, shuffle = T, seed = 1, group.by= "HTO_classification", split.by= "HTO_classification", ncol=3, raster=FALSE)
save_plot(p, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/version9_wnn_HTO_dimplot", 36, 24)
p=DimPlot(MKP3HTO, group.by= "seurat_clusters", raster=FALSE, label = TRUE)
save_plot(p, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/version9_wnn_pca_ig_dimplot", 12, 8)


####
deResTT = makeDEResults(MKP3HTO, group.by="seurat_clusters", assay="RNA", test="t")
exprdfTT = getDEXpressionDF(MKP3HTO, deResTT, assay="RNA", group.by="seurat_clusters")
write.table(exprdfTT,"/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/version9_expr_test_clusters_t.tsv", sep="\t", row.names=F, quote = F)
write_xlsx(deResTT, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/version9_expr_test_clusters_t.xlsx")


#exprdfTT<-read_tsv("singlets/expr_test_clusters_t.tsv")
DefaultAssay(MKP3HTO) <- "RNA"
markers.use.tt= subset(deResTT , avg_log2FC>0&p_val_adj<0.05&!startsWith(gene, "mt-")&!startsWith(gene, "rp"))
finalMarkers.use.tt = markers.use.tt %>% arrange(p_val_adj, desc(abs(pct.1)*abs(avg_log2FC))) %>% group_by(clusterID) %>% dplyr::slice(1:20)
finalMarkers.use.tt


data_dupli= finalMarkers.use.tt[!duplicated(finalMarkers.use.tt[ , "gene"]), ]


events= data_dupli %>% dplyr::count(clusterID)
inser=cumsum(events$n)+0.5-events$n
insert=replace(inser, inser==0.5, 0)

xmi<- insert
xmin<- xmi[c(FALSE, TRUE)]
xma<- insert+events$n
xmax<- xma[c(FALSE, TRUE)]
ymi<- 0*events$n
ymin<- ymi[c(FALSE, TRUE)]
yma<- rep(length(events$n)+0.5, each=length(events$n))
ymax<- yma[c(FALSE, TRUE)]



p_dp_genes_idents = DotPlot(MKP3HTO, features = data_dupli$gene, assay="RNA", dot.scale = 5, group.by="seurat_clusters")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
    annotate("rect", xmin=xmin, xmax=xmax, ymin=ymin , ymax=ymax, alpha=0.2, fill="blue") #rep(c("blue", "grey"), times= length(events$n)/2)
  save_plot(p_dp_genes_idents, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/version9_dotplot_cluster_genes_colored", 45, 8)





#########################Post-analysis
##
save.image("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/MKP3HTO.Rdata")
saveRDS(MKP3HTO, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/MKP3HTO.RDS")
MKP3HTO = readRDS("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/MKP3HTO.RDS")



new.cluster.ids <- c("cluster_0", "cluster_1","cluster_2","cluster_3","cluster_4","cluster_5","cluster_6","cluster_7","cluster_8","cluster_9")
names(new.cluster.ids) <- levels(MKP3HTO)
MKP3HTO <- RenameIdents(MKP3HTO, new.cluster.ids)

MKP3HTO$Idents = Idents(MKP3HTO)
MKP3HTO_new <- subset(x = MKP3HTO, idents = c("cluster_0", "cluster_1","cluster_2","cluster_4","cluster_8"))

deGroup = "CSclassification"
allConds = c("Bl6","Platelet-depletion","pDC-Platelet-depletion")
print(allConds)

all.cells = list()

for (cond in allConds)
{
  all.cells[[cond]] = cellIDForClusters(MKP3HTO_new, deGroup, c(cond)) 
}



for (i in 1:(length(allConds)-1))
{

  for (j in (i+1):length(allConds))
  {
    print(paste(i,"<->",j))

    condI = allConds[i]
    condJ = allConds[j]

    condNameI = tolower(str_replace_all(condI, " ", "_"))
    condNameJ = tolower(str_replace_all(condJ, " ", "_"))

    print(paste(condJ, condNameJ))
    print(paste(condI, condNameI))

    deMarkers= de.condi.condj = compareCellsByCluster(MKP3HTO_new, all.cells[[condJ]], all.cells[[condI]], condNameJ, condNameI,
                                                outfolder=paste("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/Pseudobulk_01248", "de", sep="/"), fcCutoff=0.1)


    deMarkersW= dewilcox.condi.condj = compareCellsByCluster(MKP3HTO_new, all.cells[[condJ]], all.cells[[condI]], condNameJ, condNameI,
                                                    outfolder=paste("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/Pseudobulk_01248", "dewilcox", sep="/"), test="wilcox", fcCutoff=0.1)

    compName = paste(condNameJ, condNameI, sep="_")

    saveRDS(de.condi.condj, file=paste("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/Pseudobulk_01248", "de", paste("de_", compName, ".rds", sep=""), sep="/"))
    saveRDS(dewilcox.condi.condj, file=paste("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/Pseudobulk_01248", "dewilcox", paste("der_", compName, ".rds", sep=""), sep="/"))
    
    makeVolcanos(de.condi.condj, paste("DE", condJ, "vs", condI), paste("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/Pseudobulk_01248", "de_volcano", compName, sep="/"),
                turnExpression=F, FCcutoff=0.1, pCutoff = 0.05)

    makeVolcanos(dewilcox.condi.condj, paste("DE", condJ, "vs", condI), paste("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/Pseudobulk_01248", "dewilcox_volcano", compName, sep="/"),
                turnExpression=F, FCcutoff=0.1, pCutoff = 0.05)
    }

}



###### trajectory analysis monocle3

library(monocle3)
MKP3HTO = readRDS("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/MKP3HTO.RDS")


cds <- order_cells(cds, root_pr_nodes='Y_40')
p = plot_cells(cds, color_cells_by = "pseudotime", label_branch_points=TRUE, label_leaves=TRUE)
save_plot(p, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/monocle3_MKP3/version_4/monocle3_MKP3_pseudotime_Y40", 15, 15)


# ...1 Convert to cell_data_set object ------------------------
new.cluster.ids <- c("0", "1","2","3","4","5","6","7","8","9")
names(new.cluster.ids) <- levels(MKP3HTO)
MKP3HTO <- RenameIdents(MKP3HTO, new.cluster.ids)

MKP3HTO$Idents = Idents(MKP3HTO)
MKP3HTO <- subset(MKP3HTO, idents = c("0","1","2","4","8"))

cds <- SeuratWrappers::as.cell_data_set(MKP3HTO)
cds

# to get cell metadata
colData(cds)
# to gene metdata
fData(cds)
rownames(fData(cds))[1:10]

# since it misses the gene_short_name column, let's add it
fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
counts(cds)



# ...2. Cluster cells (using clustering info from seurat's UMAP)---------------------------
# let's use the clustering information have

# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)


cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 

list_cluster <- MKP3HTO@active.ident
cds@clusters$UMAP$clusters <- list_cluster


# Assign UMAP coordinate - cell embeddings

cds@int_colData@listData$reducedDims$UMAP <- MKP3HTO@reductions$umap@cell.embeddings



# plot

cluster.before.trajectory <- plot_cells(cds,
           color_cells_by = 'cluster',
           label_groups_by_cluster = FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size = 5) +
  theme(legend.position = "right")
save_plot(cluster.before.trajectory, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/monocle3_MKP3/version_4/cluster", 12, 8)  


cluster.names <- plot_cells(cds,
           color_cells_by = "partition",
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  theme(legend.position = "right")
save_plot(cluster.names, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/monocle3_MKP3/version_4/cluster_partition", 12, 8)  




# ...3. Learn trajectory graph ------------------------
cds <- learn_graph(cds, use_partition = FALSE)

p = plot_cells(cds,
           color_cells_by = 'Idents',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)
save_plot(p, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/monocle3_MKP3/version_4/cluster_graph", 12, 8)  



# ...4. Order the cells in pseudotime -------------------

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == 1]))

p = plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)
save_plot(p, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/monocle3_MKP3/version_4/pseudotime", 12, 8)  



# cells ordered by monocle3 pseudotime

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

p = ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(Idents, monocle3_pseudotime, median), fill = Idents)) +
  geom_boxplot()
save_plot(p, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/monocle3_MKP3/version_4/pseudotime_monocle3", 12, 8)  




# ...5. Finding genes that change as a function of pseudotime --------------------
deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)
deg_bcells = deg_bcells[order(deg_bcells$morans_I, decreasing = TRUE), ]
genes = head(rownames(deg_bcells), 200)
deg_bcells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()

png("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/monocle3_MKP3/version_4/de_genes_v4.png", 1700, 15000)
FeaturePlot(MKP3HTO, features = genes)
dev.off()

install.packages("pheatmap")
library(pheatmap)


gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))
partitions = partitions(cds, reduction_method = "UMAP")
partitions
clusters = clusters(cds, reduction_method = "UMAP")
clusters

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=clusters(cds)[colnames(cds)])
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("clusters ", colnames(agg_mat))
png(file = "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/monocle3_MKP3/version_4/HEATMAP_CLUSTERS_v4.png", 450, 650)
pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)
dev.off()                   

p = plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% 1:17),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
save_plot(p, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/monocle3_MKP3/version_4/modules_v4", 20, 15)  

x= as.data.frame(gene_module_df)
rownames(x) = x$id
genes_module_9_10_11 = rownames(x[x$module %in% c(9,10,11), ])
png("/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/monocle3_MKP3/version_4/genes_in_modules_1_4_8.png", 1700, 13000)
FeaturePlot(MKP3HTO, features = genes_module_1_4_8)
dev.off()       


# visualizing pseudotime in seurat

MKP3HTO$pseudotime <- pseudotime(cds)
Idents(MKP3HTO) <- MKP3HTO$Idents
p =FeaturePlot(MKP3HTO, features = "pseudotime", label = T)
save_plot(p, "/mnt/raidtmp/Alejandro/velocyto_florian/soup_scRNAseq/without_ribosomal_genes/monocle3_MKP3/version_4/pseudotime_monocle3_function", 12, 8)  






