############# scRNAseqAnalysis.R #######################
#### processes aligned scRNA-seq for Sars-CoV-2 paper

################### Packages ###########################
library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(cowplot)

################### Settings ###########################
options(future.globals.maxSize = 8000 * 1024^2) #important for SCT
dataDir = "" #wherever the folder of sample folders of matrix files are 
homedir = ""
SPECIES = "HUMAN"

################### Constants ###########################
DARKCOLORS <- c(RColorBrewer::brewer.pal(8,"Set1"), RColorBrewer::brewer.pal(8,"Dark2"), "#FFFFFF")
LIGHTGRAY = "#D4D3D2"
DARKBLUE = "#000099"
COLORSCHEME2 = c("#3E80A0", "#5E3EA0", "#A03E80", "#3ea05e", "#a05e3e", "#A1903F", "#a23f3f" ,"#3fa27a", "#a23f68")
HEMOGENES.MMU <- c("Hbb-bt", "Hbb-bs", "Hba-a1", "Hba-a2", "Hbq1a", "Hbq1b")
HEMOGENES.HSA <- c("HBB", "HBA1", "HBA2")
if (SPECIES == "HUMAN") {
  HEMOGENES = HEMOGENES.HSA
} else if (SPECIES == "MOUSE") {
  HEMOGENES = HEMOGENES.MMU 
} else {
  print("only mmu and hsa supported for species")
}

############### Functions #####################

pp <- function(seurat_object, HEMOGENES=NULL, subset=TRUE, SCT = FALSE, VAR.FEATURES=3000, MAX_FEATURES_PER_CELL=2500, MIN_FEATURES_PER_CELL=100, MAX_PERCENT_MT=5, DIMS=40) {
  genes <- rownames(seurat_object@assays$RNA@counts)
  genes.mt <- genes[grep("^M[Tt]-", genes)]
  print(paste0("found ", length(genes.mt), " genes: ", paste(genes.mt, collapse = " ")))
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^M[Tt]-")
  
  genes.ribo <- genes[grep("^R[Pp][Ss]|R[Pp][Ll]", genes)]
  print(paste0("found ", length(genes.ribo), " genes: ", paste(genes.ribo, collapse = " ")))
  seurat_object[["percent.ribo"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^R[Pp][Ss]|R[Pp][Ll]")
  
  if (is.null(HEMOGENES)) {
    genes.hemo <- genes[grep("^H[Bb][Bb]|^H[Bb][Aa][12]", genes)]
  } else genes.hemo <- HEMOGENES
  print(paste0("found ", length(genes.hemo), " genes: ", paste(genes.hemo, collapse = " ")))
  seurat_object[["percent.hbb"]] <- PercentageFeatureSet(object = seurat_object, features = genes.hemo)
  featureinfo <- seurat_object@meta.data[,c("nCount_RNA", "nFeature_RNA", "orig.ident", "percent.mt", "percent.ribo", "percent.hbb")]
  #subset on filtered cells and continue
  if (subset) {
    seurat_object <- subset(x = seurat_object, subset = nFeature_RNA > MIN_FEATURES_PER_CELL & nFeature_RNA < MAX_FEATURES_PER_CELL & percent.mt < 5)
  }
  if (SCT) {
    seurat_object <- SCTransform(seurat_object, vars.to.regress = "percent.mt") #use this to replace scaledata if it's not good enough
  } else {
    seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 1e4)
    seurat_object <- FindVariableFeatures(object = seurat_object, selection.method = 'vst', nfeatures = VAR.FEATURES)
    seurat_object <- ScaleData(object = seurat_object, vars.to.regress = c('percent.mt') )
  }
  seurat_object <- RunPCA(object = seurat_object, features = VariableFeatures(object = seurat_object))
  seurat_object <- RunUMAP(object = seurat_object, dims=1:DIMS) #VariableFeatures(seurat_object) makes it error -> need dims for some reason 
  seurat_object <- FindNeighbors(object = seurat_object, dims=1:DIMS)
  seurat_object <- FindClusters(object = seurat_object, dims=1:DIMS, resolution = 0.4)
  message("saving...")
  saveRDS(seurat_object, paste0(seurat_object@project.name, "_pp.rds"))
  return(seurat_object)
}


allUMAPs <- function(seurat_object) {
  UMAPPlot(seurat_object, group.by="seurat_clusters", label=TRUE)
  ggsave("umap_byCluster.png", device = "png", units = "in", width = 7, height = 6)
  UMAPPlot(seurat_object, group.by="condition")
  ggsave("umap_byCondition.png", device = "png", units = "in", width = 7, height = 6)
  UMAPPlot(seurat_object, group.by="orig.ident")
  ggsave("umap_bySample.png", device = "png", units = "in", width = 7, height = 6)
}

clusterDE <- function(seurat_object, name, cluster_col=NULL, assay="RNA") {
  library(dplyr)
  if (is.null(cluster_col)) {
    Idents(seurat_object) <- seurat_object@meta.data$seurat_clusters
  } else {
    Idents(seurat_object) <- seurat_object@meta.data[,cluster_col]
  }
  seurat_object.markers <- FindAllMarkers(seurat_object, min.pct = 0.25, logfc.threshold = 0.25, random.seed = 777)
  write.table(seurat_object.markers, paste0("markers_", name, ".txt"), sep = "\t", quote = F, col.names = NA)
  top10 <- seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  DoHeatmap(seurat_object, features = as.vector(top10$gene)) 
  height80dpi <- (length(as.vector(top10$gene))*12)*(1/80)+1
  ggsave(paste0("heatmap_", name, ".png"), height= height80dpi, width = 8)
}

PercentAbove <- function(x, threshold=1){
  return(length(x = x[x > threshold]) / length(x = x))
}

getDotplotData <- function(seurat_object, features) {
  data.features <- FetchData(object = seurat_object, vars = features)
  data.features$id <- Idents(seurat_object)
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      print(data.use)
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  return(data.plot)
}

#Gets avg.exp and pct.exp for each combination in each identity from seurat object
# avg.exp is actually just mean over feat1's expression in that ident
# pct.exp is percent of both being expressed over that ident
#featureCombos <- c("ACE2_TMPRSS2", "ACE2_CTSB", "ACE2_CTSL")
coexpressionDotplotDat <- function(seurat_object, featureCombos) {
  coexpression_dat <- c()
  for (fc in featureCombos) {
    feat1 <- unlist(strsplit(fc, "_"))[[1]]
    feat2 <- unlist(strsplit(fc, "_"))[[2]]
    if (!(feat1 %in% rownames(seurat_object@assays$RNA@data)) || !(feat2 %in% rownames(seurat_object@assays$RNA@data))) {
      if (!(feat1 %in% rownames(seurat_object@assays$RNA@data))) {
        print(paste0(feat1, " not in object"))
      } else {
        print(paste0(feat2, " not in object"))
      }
      next
    }
    fd1 <- FetchData(seurat_object, feat1)
    fd2 <- FetchData(seurat_object, feat2)
    fd.both <- cbind.data.frame(fd1, fd2)
    fd.both$id <- Idents(seurat_object)
    fd.both$pct.exp <- NA
    fd.both$avg.exp <- NA
    for (id in unique(fd.both$id)) {
      frame <- fd.both[which(fd.both$id == id),]
      pct <- (calcPctCoexp(frame[,1], frame[,2]) / nrow(frame))
      fd.both$pct.exp[which(rownames(fd.both) %in% rownames(frame))] <- pct
      avg <- mean(frame[,feat1])
      fd.both$avg.exp[which(rownames(fd.both) %in% rownames(frame))] <- avg
    }
    fd.both.pre <- fd.both[,c((ncol(fd.both)-2), (ncol(fd.both)-1), ncol(fd.both))]
    fd.both.final <- unique(fd.both.pre[,1:3])
    fd.both.final$features.plot <- fc
    fd.both.final <- fd.both.final[,c(3,2,4,1)]
    rownames(fd.both.final) <- paste0(fd.both.final$features.plot, 1:length(fd.both.final$features.plot))
    coexpression_dat <- rbind(coexpression_dat, fd.both.final)
  }
  return(coexpression_dat)
}

calcPctCoexp <- function(x,y) {
  count <- 0
  for (i in length(x)) {
    if (min(x[i], y[i]) > 0) {
      count = count + 1
    }
  }
  return(count)
}

#combined expression dotplot
sarsCov2Coexpression <- function(seurat_object, outdir) {
  sarsCov2.features.gene <- c("TMPRSS2", "CTSB", "CTSL", "ACE2")
  featureCombos <- c("ACE2_TMPRSS2", "ACE2_CTSB", "ACE2_CTSL")
  notin <- sarsCov2.features.gene[!(sarsCov2.features.gene %in% rownames(seurat_object@assays$RNA@data))]
  data.plt <- getDotplotData(seurat_object, sarsCov2.features.gene)
  temp <- c()
  for (f in notin) {
    for (id in unique(data.plt$id)) {
      filler <- c("avg.exp"=0, "pct.exp"=0, "features.plot"=f, "id"=id)
      temp <- rbind(temp, filler)
    }
    rownames(temp) <- paste0(f,1:length(unique(data.plt$id)))
  }
  data.plt <- rbind(data.plt, temp)
  data.plt$avg.exp <- as.numeric(data.plt$avg.exp)
  data.plt$pct.exp <- as.numeric(data.plt$pct.exp)
  data.coexp.plt <- coexpressionDotplotDat(seurat_object, featureCombos)
  data.plt.all <- rbind(data.plt, data.coexp.plt)
  data.plt.all$pct.exp[which(data.plt.all$pct.exp == 0)] <- NA
  data.plt.all$avg.exp[which(data.plt.all$avg.exp == 0)] <- NA
  data.plt.all$features.plot <- gsub("_", " and ", data.plt.all$features.plot)
  data.plt.all$features.plot <- factor(data.plt.all$features.plot, 
                                       levels = c(sarsCov2.features.gene, gsub("_", " and ", featureCombos)))
  ggplot(data.plt.all, aes(x=features.plot, y=id)) + 
    geom_point(aes(size = pct.exp, color = avg.exp)) +
    # geom_point(aes(size = pct.exp), shape = 1, color = "#000000") +
    facet_wrap(~features.plot, ncol=length(unique(data.plt.all$features.plot)), scales="free_x") +
    scale_radius(breaks=c(1,5,10,40), trans = "log2") +
    scale_color_viridis( trans = "log2") +
    xlab("Feature") + 
    ylab("Cluster and Condition") +
    labs(color="Average Expression", size="Percent Expressed") +
    theme(axis.text.x = element_text(hjust=1,vjust=1,angle=45,face=2,size=9), 
          axis.text.y = element_text(face=2,size=9), 
          axis.title = element_text(face=2,size=12), 
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill="#FFFFFF",color="#000000", size=1),
          strip.background = element_blank(), strip.text = element_blank())
  ggsave(paste0(outdir, "/sars_coexpression_dotplot.pdf" ))
}

#preprocess using method from paper
filterForComplexity <- function(x) {
  top200 <- sort(x, decreasing = T)[1:200]
  if (sum(top200)/sum(x) < .85) {
    return(TRUE)
  } else return(FALSE)
}

newtsne <- function(seed) {
  inds123 <- RunTSNE(inds123, reduction = "pca", dims=1:2, seed.use = 36)
  inds123 <- FindNeighbors(inds123, reduction = "pca", dims=1:2)
  inds123 <- FindClusters(inds123, resolution = 0.1, random.seed = 0)
  TSNEPlot(inds123)
}

####################################################


# importing matrix.mtx, setting up Experiment and Creating objects
# get sample names from matrix path x/y/sampleName/matrix.mtx
## import all the data
setwd("/path/to/C_sec/and/GRP/countsdata")
seurat_project <- list()
for (f in list.files(dataDir)[grep("C_sec|GRP", list.files(dataDir))]) { #changed to .jfgc
  print(paste0(f, "/", list.files(paste0(dataDir, "/", f, "/Solo.out")))) #assuming we use STAR solo
  curfolder <- paste0(dataDir, "/", f,"/Solo.out")
  if (length(list.files(curfolder)) > 0) {
    print(paste0("importing: ", f))
    temp.data <- Read10X(curfolder)
    seurat_project[[f]] <- CreateSeuratObject(counts = temp.data, project = f, min.cells = MIN_CELLS_PER_FEATURE, min.features = MIN_FEATURES_PER_CELL)
  }
}

#### pp from the paper ####
# seurat_project.paper <- seurat_project[c(1,4)] #1 is C_sec, 4 is GRP
seurat_project.paper <- seurat_project #revised so the only objects int the project are C_sec and GRP
for (i in 1:length(seurat_project.paper)) {
  seurat_project.paper[[i]][["percent.mt"]] <- PercentageFeatureSet(object = seurat_project.paper[[i]], pattern = "^MT-")
  seurat_project.paper[[i]] <- subset(seurat_project.paper[[i]], subset = nFeature_RNA > MIN_FEATURES_PER_CELL & nFeature_RNA < MAX_FEATURES_PER_CELL &
                            percent.mt < 25) #wow this is really high
  seurat_project.paper[[i]] <- NormalizeData(seurat_project.paper[[i]], normalization.method = "LogNormalize", scale.factor = 1e4, verbose=F)
  seurat_project.paper[[i]] <- FindVariableFeatures(object = seurat_project.paper[[i]], selection.method = 'vst', nfeatures = VAR.FEATURES)
}
seurat_project.anchors <- FindIntegrationAnchors(object.list = seurat_project.paper, dims = 1:30)
ovary_unsorted.int <- IntegrateData(anchorset = seurat_project.anchors, dims = 1:30)
ovary_unsorted.int <- ScaleData(ovary_unsorted.int, verbose = F)
ovary_unsorted.int <- RunPCA(object = ovary_unsorted.int, npcs=13, verbose=F)
ovary_unsorted.int <- RunUMAP(object = ovary_unsorted.int, reduction = "pca", dims=1:13)
ovary_unsorted.int <- FindNeighbors(object = ovary_unsorted.int, dims=1:13)
ovary_unsorted.int <- FindClusters(object = ovary_unsorted.int, resolution = 0.1)

# ovary_unsorted.int <- readRDS("saves/ovary_unsorted_filtered_integrated_37982cells.rds")
ovary_unsorted.int.markers <- FindAllMarkers(ovary_unsorted.int, assay = "integrated", min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
write.table(ovary_unsorted.int.markers, paste0(homedir, "R_Analysis/ovary_unsorted_int_markers.txt"), 
            sep = "\t", quote = F, col.names = NA)
ovary_unsorted.int.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- as.vector(top10$gene)

#heatmap
DoHeatmap(ovary_unsorted.int, features = unique(top10genes)) #do this
height80dpi <- (length(top10genes)*12)*(1/80)+1 #it's a pdf so it ends up being fine (for dpi)
ggsave("ovary_unsorted_int_heatmap.pdf", height=height80dpi, width=8)

#dotplot
COLORSCHEME2 = c("#3E80A0", "#5E3EA0", "#A03E80", "#3ea05e", "#a05e3e", "#A1903F", "#a23f3f" ,"#3fa27a", "#a23f68")
DotPlot(ovary_unsorted.int, features = unique(top10genes), cols = c(LIGHTGRAY, COLORSCHEME2[3])) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1))
ggsave("ovary_unsorted_int_dotplot.pdf", width = 20, height=6)

ovary_unsorted.int@meta.data$condition <- ovary_unsorted.int@meta.data$orig.ident
UMAPPlot(ovary_unsorted.int, group.by="seurat_clusters", label=TRUE)
ggsave("umap_byCluster.png", device = "png", units = "in", width = 7, height = 6)
UMAPPlot(ovary_unsorted.int, group.by="condition")
ggsave("umap_byCondition.png", device = "png", units = "in", width = 7, height = 6)
UMAPPlot(ovary_unsorted.int, group.by="orig.ident")
ggsave("umap_bySample.png", device = "png", units = "in", width = 7, height = 6)

#### After doing find markers on the clusters 0:8, 0 and 1 are the same celltype, merging.
ovary_unsorted.int@meta.data$seurat_clusters[which(ovary_unsorted.int@meta.data$seurat_clusters == 1)] <- 0
Idents(ovary_unsorted.int) <- ovary_unsorted.int@meta.data$seurat_clusters
ovary_unsorted.int <- RenameIdents(ovary_unsorted.int, '0'="stroma", '2'="pv", '3'="endo", '4'="immune_nk_t", '5'="granulosa", '6'="immune", 
                                   '7'="endo", '8'="oocyte")
ovary_unsorted.int@meta.data$celltype <- Idents(ovary_unsorted.int)
# saveRDS(ovary_unsorted.int, "saves/ovary_unsorted_int_with_celltype_names.rds")

#### UPDATE: changing naming from 0,2,3,4... to 0,1,2,3...  2020-05-14
ovary_unsorted.int@meta.data$seurat_clusters[which(ovary_unsorted.int@meta.data$seurat_clusters == 2)] <- 1
ovary_unsorted.int@meta.data$seurat_clusters[which(ovary_unsorted.int@meta.data$seurat_clusters == 3)] <- 2
ovary_unsorted.int@meta.data$seurat_clusters[which(ovary_unsorted.int@meta.data$seurat_clusters == 4)] <- 3
ovary_unsorted.int@meta.data$seurat_clusters[which(ovary_unsorted.int@meta.data$seurat_clusters == 5)] <- 4
ovary_unsorted.int@meta.data$seurat_clusters[which(ovary_unsorted.int@meta.data$seurat_clusters == 6)] <- 5
ovary_unsorted.int@meta.data$seurat_clusters[which(ovary_unsorted.int@meta.data$seurat_clusters == 7)] <- 6
ovary_unsorted.int@meta.data$seurat_clusters[which(ovary_unsorted.int@meta.data$seurat_clusters == 8)] <- 7
UMAPPlot(ovary_unsorted.int, group.by = "seurat_clusters", label=T)
ggsave("umap_by_cluster.pdf")

#### sars-cov2 genes ####
sarsCov2.features <- c("TMPRSS2", "CTSB", "CTSL", "ACE2")
DotPlot(ovary_unsorted.int, assay = "RNA", features=sarsCov2.features, group.by = "celltype", cols = c(LIGHTGRAY, COLORSCHEME2[4])) + #facet_wrap(~ seurat_clusters, ncol=1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1))
ggsave("sars_cov2/dotplot.pdf")

#Feature plots and umap to confirm.
UMAPPlot(ovary_unsorted.int, group.by = "celltype")
ggsave("umap_by_celltype.pdf")
UMAPPlot(ovary_unsorted.int, group.by = "seurat_clusters")
ggsave("umap_by_cluster.pdf")
DefaultAssay(ovary_unsorted.int) <- "RNA"

DARKBLUE = "#000099"
FeaturePlot(ovary_unsorted.int, features = sarsCov2.features,cols = c(LIGHTGRAY, DARKBLUE))
ggsave(paste0(homedir, "umap_features_sars.pdf"), width = 8, height = 7.82)


# REFERENCE FROM GITHUB #
ovary_unsorted.int <- RenameIdents(ovary_unsorted.int, '5' = '1_oocytes', '4' = '2_immune', '3' = '3_gran', '2' = '4_endo', 
                                   '1' = '5_pv', '0' = '6_stroma')




#### end pp from paper####



#### falopian paper ####
# let's just try taking the counts 
# QC
fallopian_counts <- read.table(paste0(homedir, "fallopian/GSE139079_counts_benign.txt", sep = "\t"))
ft <- CreateSeuratObject(fallopian_counts, project = "ft", min.cells = 0, min.features = 0) 
ft[["percent.mito"]] <- PercentageFeatureSet(object = ft, pattern = "^MT") # mito genes are mixed in, no dash after "MT"
ft@meta.data$orig.ident <- "fallopian_tube"
hist(ft@meta.data$nCount_RNA, breaks=200, main = "Histogram of UMI") 
hist(ft@meta.data$nFeature_RNA, breaks=200, main = "Histogram of Genes")
#flawed because percent.mito doesn't only have mitochdondrial genes, but it's soooo low anyway (max is 1.21%)
hist(ft@meta.data$percent.mito, breaks=200, main = "Histogram of Percent Mitochondrial Reads") 


fallopian_counts <- read.table(paste0(homedir, "fallopian/GSE139079_counts_benign.txt"), sep = "\t")
fallopian_counts <- fallopian_counts[,apply(fallopian_counts, 2, filterForComplexity)] #same number
ft <- CreateSeuratObject(fallopian_counts, project = "ft", min.cells = 0, min.features = 0) 
ft[["percent.mito"]] <- PercentageFeatureSet(object = ft, pattern = "^MT") # mito genes are mixed in, no dash after "MT"
ft@meta.data$orig.ident <- gsub("sc_|_[A-Z]+.*$", "",rownames(ft@meta.data))
Idents(ft) <- ft@meta.data$orig.ident
ft <- subset(ft, subset = nFeature_RNA > 1200 & nFeature_RNA < 7500 & nCount_RNA > 100000) #1857 before, 1708 after
ft <- NormalizeData(ft, normalization.method = "LogNormalize", scale.factor = 1e4, verbose=F)
ft <- ScaleData(ft) 
ft <- FindVariableFeatures(ft, nfeatures = 2000)
ft <- RunPCA(ft, npcs = 50, nfeatures.print = 5, ndims.print = 1:5, features = ft@assays$RNA@var.features) 
ft <- RunUMAP(ft, dims = 1:20) 
UMAPPlot(ft, label=T)
ft <- FindNeighbors(ft)
ft <- FindClusters(ft, resolution = 0.2)
ft.markers <- FindAllMarkers(ft, assay = "RNA", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1)
ft.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- as.vector(top10$gene)
DoHeatmap(ft, features = unique(top10genes)) #do this
height80dpi <- (length(top10genes)*12)*(1/80)+1 #it's a pdf so it ends up being fine (for dpi)
ggsave(paste0(homedir, "fallopian/R_Analysis/benign_heatmap.pdf"))
saveRDS(ft, paste0(homedir, "fallopian/R_Analysis/fallopianTubeBenign_processed_1708total.rds"))


#### annotation after
ft@meta.data$celltype <- NA
ft@meta.data$celltype[which(ft@meta.data$seurat_clusters == 0)] <- "Secretory"
ft@meta.data$celltype[which(ft@meta.data$seurat_clusters == 1)] <- "Ciliated"
ft@meta.data$celltype[which(ft@meta.data$seurat_clusters == 2)] <- "Secretory"
ft@meta.data$celltype[which(ft@meta.data$seurat_clusters == 3)] <- "Leukocytes"
ft@meta.data$celltype[which(ft@meta.data$seurat_clusters == 4)] <- "Stromal"
setwd(paste0(homedir, "fallopian/R_Analysis/afterannotation"))

sarsCov2.features.gene <- c("TMPRSS2", "CTSB", "CTSL", "ACE2")
DotPlot(ft, assay = "RNA", features=sarsCov2.features.gene, group.by = "celltype", cols = c(LIGHTGRAY, COLORSCHEME2[4])) + #facet_wrap(~ seurat_clusters, ncol=1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1))
ggsave(paste0(homedir, "fallopian/R_Analysis/afterannotation/dotplot_byCluster.pdf"), width=10, height=10)

#Feature plots and umap to confirm.
UMAPPlot(ft, group.by = "celltype")
  # theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), 
  #       panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18), 
  #       axis.line.x = element_line(arrow  = arrow(angle = 30, length=unit(0.08, "inches"))), 
  #       axis.line.y = element_line(arrow  = arrow(angle = 30, length=unit(0.08, "inches")))) 
ggsave(paste0(homedir, "fallopian/R_Analysis/afterannotation/umap_byCluster.pdf"), width=10, height=10)
UMAPPlot(ft, group.by = "orig.ident", label=T)
ggsave(paste0(homedir, "fallopian/R_Analysis/afterannotation/umap_bySample.pdf"), width=10, height=10)

FeaturePlot(ft, features = sarsCov2.features.gene, pt.size = 0.5)
ggsave(paste0(homedir, "fallopian/R_Analysis/afterannotation/sars_features.pdf"), width=8, height=6)

#other feature plot (more customizable)
DefaultAssay(ft) <- "RNA"
plots <- list()
for (i in 1:length(sarsCov2.features.gene)) {
  if (!(sarsCov2.features.gene[i] %in% rownames(ft@assays$RNA@data))) {
    print(paste0("gene ", sarsCov2.features.gene[i], " not in RNA@data, continuing..."))
    next
  }
  cur.expressiondata <- ft@assays$RNA@data[which(rownames(ft@assays$RNA@data) == sarsCov2.features.gene[i]),]
  cur.plotdata <- cbind(ft@meta.data, "exprs"=cur.expressiondata)
  cur.plotdata.umap <- merge(cur.plotdata, ft@reductions$umap@cell.embeddings, by='row.names')
  
  #umap
  cur.plotdata.umap.order <- cur.plotdata.umap[order(cur.plotdata.umap$exprs,decreasing = F),]
  plots[[sarsCov2.features.gene[i]]] <- ggplot(data=cur.plotdata.umap, aes_string(x="UMAP_1", y="UMAP_2")) +
    geom_point(aes(color=exprs), size=0.3) + scale_colour_gradient(low = LIGHTGRAY,high = DARKBLUE, ) +
    ggtitle(sarsCov2.features.gene[i]) + 
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), 
          panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18)) 
}
if (length(plots) == 3) {
  plot_grid(plots[[1]], NULL, plots[[2]] + 
              theme(axis.title = element_text(), 
              axis.line.x = element_line(arrow  = arrow(angle = 30, length=unit(0.08, "inches"))), 
              axis.line.y = element_line(arrow  = arrow(angle = 30, length=unit(0.08, "inches")))),
          plots[[3]], ncol = 2, nrow = 2) 
} else if (length(plots) == 4){
  plot_grid(plots[[1]], plots[[2]], plots[[3]] #+ theme for arrows
            , plots[[4]], ncol = 2, nrow = 2)  
}
ggsave(paste0(homedir, "fallopian/R_Analysis/sars_features.pdf"), width=8, height=6)
cluster_col <- "celltype"
name = "fallopian_tube_benign"
if (is.null(cluster_col)) {
  Idents(ft) <- ft@meta.data$seurat_clusters
} else {
  Idents(ft) <- ft@meta.data[,cluster_col]
}
ft.markers <- FindAllMarkers(ft, min.pct = 0.25, logfc.threshold = 0.25, random.seed = 777)
write.table(ft.markers, paste0("markers_", name, ".txt"), sep = "\t", quote = F, col.names = NA)
top10 <- ft.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(ft, features = as.vector(top10$gene)) 
height80dpi <- (length(as.vector(top10$gene))*12)*(1/80)+1
ggsave(paste0("heatmap_", name, ".png"), height= height80dpi, width = 8)
Idents(ft) <- ft@meta.data$celltype
x <- DotPlot(ft, assay = "RNA", features=sarsCov2.features.gene, group.by = "celltype") + 
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1))


#### end falopian paper data ####



#### breast epithelial cell paper ####

####### individuals 1-3
setwd(paste0(homedir, "breastepi/R_Analysis"))
datadirs <- paste0(homedir, "breastepi/", 
                   c("GSE113197_RAW/", "GSE113099_RAW/", "GSE113198_RAW/", "GSE113127_RAW/"))
for (d in datadirs) {
  print(d)
  samples.files <- list.files(d)
  counts.table <- c()
  for (i in 1:length(samples.files)) {
    file.location <- paste0(d, samples.files[i])
    if (i == 1) {
      genes.names <- read.table(file.location)
      genes.names <- genes.names[-1,1]
      counts.table <- as.vector(genes.names)
    }
    sample.name <- paste(unlist(strsplit(samples.files[i], "_"))[c(2,3,4)], collapse = "_")
    fpkm <- read.table(file.location)
    fpkm <- fpkm[-1,2]
    counts.table <- cbind(counts.table, as.vector(as.matrix(fpkm)))
    colnames(counts.table)[ncol(counts.table)] <- sample.name
  }
  rownames(counts.table) <- counts.table[,1]
  counts.table <- counts.table[,-1]
  counts.table <- as.data.frame(counts.table)
  saveRDS(counts.table, paste0("saves/counts_", unlist(strsplit(d, "/"))[length(unlist(strsplit(d, "/")))],".rds"))
}

# all of these are full expression matrices
datadir <- paste0(homedir, "breastepi/")
#human_ens_GRCh38_annot.extended.txt is taken from ensemble's biomart with gene symbol, ensembl id, and gene description
annot <- read.table("~/data/human_ens_GRCh38_annot.extended.txt", sep="\t", quote="", header=T, row.names=1, stringsAsFactors=F, fill=T)
annot <- annot[,c(1,2,3)] #just for gene transfer

counts_GSE113099_RAW <- readRDS("saves/counts_GSE113099_RAW.rds") # individual 1 
counts_GSE113127_RAW <- readRDS("saves/counts_GSE113127_RAW.rds") # individual 2 
counts_GSE113198_RAW <- readRDS("saves/counts_GSE113198_RAW.rds") # individual 3 
counts_GSE113197_RAW <- readRDS("saves/counts_GSE113197_RAW.rds")

so.GSE113099_RAW <- CreateSeuratObject(counts = counts_GSE113099_RAW, project = "GSE113099", min.cells = 3, min.features = 900)
so.GSE113127_RAW <- CreateSeuratObject(counts = counts_GSE113127_RAW, project = "GSE113127", min.cells = 3, min.features = 900)
so.GSE113198_RAW <- CreateSeuratObject(counts = counts_GSE113198_RAW, project = "GSE113198", min.cells = 3, min.features = 900)
inds123 <- merge(so.GSE113099_RAW, list(so.GSE113127_RAW, so.GSE113198_RAW), 
                 add.cell.ids = c("so.GSE113099_RAW", "so.GSE113127_RAW", "so.GSE113198_RAW"),
                 project = "first3individuals")
FeatureScatter(inds123, feature1="nCount_RNA", feature2 = "nFeature_RNA")
inds123 <- NormalizeData(inds123, normalization.method = "LogNormalize", scale.factor = 10000)
inds123 <- FindVariableFeatures(inds123, nfeatures = 2500)
inds123 <- ScaleData(inds123)
inds123 <- RunPCA(inds123, npcs=10, ndims.print = 2, seed.use = 42)
newtsne(36)

#adjustments to figure -> red and green are actually correct colors
FeaturePlot(inds123, "ENSG00000111057") #KRT18
ggsave("ind123_feature_krt18.pdf")
FeaturePlot(inds123, "ENSG00000186847") #KRT14
ggsave("ind123_feature_krt14.pdf")

inds123.markers <- FindAllMarkers(inds123, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
inds123.markers.out <- merge(inds123.markers, annot, by = "row.names")
inds123.markers.out <- inds123.markers.out[,c(1,(ncol(inds123.markers.out)-1), ncol(inds123.markers.out), 2:(ncol(inds123.markers.out)-2)) ]
inds123.markers.out <- inds123.markers.out[,-c((ncol(inds123.markers.out)-1):(ncol(inds123.markers.out)))]
rownames(inds123.markers.out) <- inds123.markers.out[,1]
inds123.markers.out <- inds123.markers.out[,-1]
write.table(inds123.markers.out, paste0(homedir, "breastepi/R_Analysis/inds123_markers.txt"), 
            sep = "\t", quote = F, col.names = NA)
inds123.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- as.vector(top10$gene)

#heatmap
DoHeatmap(inds123, features = unique(top10genes)) #do this
height80dpi <- (length(top10genes)*12)*(1/80)+1 #it's a pdf so it ends up being fine (for dpi)
ggsave(paste0(homedir, "breastepi/R_Analysis/inds123_heatmap.pdf"))

#sars dotplot
sarsCov2.features.gene <- c("TMPRSS2", "CTSB", "CTSL", "ACE2")
sarsCov2.features <- c("ENSG00000184012", "ENSG00000164733", "ENSG00000135047", "ENSG00000130234")
DotPlot(inds123, assay = "RNA", features=sarsCov2.features, group.by = "seurat_clusters", cols = c(LIGHTGRAY, COLORSCHEME2[4])) + #facet_wrap(~ seurat_clusters, ncol=1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1))
ggsave(paste0(homedir, "breastepi/R_Analysis/sars_cov2/sars_dotplot.pdf"))

#Feature plots and umap to confirm.
TSNEPlot(inds123, group.by = "seurat_clusters")
ggsave(paste0(homedir, "breastepi/R_Analysis/tsne_by_cluster.pdf"))
TSNEPlot(inds123, group.by = "orig.ident")
ggsave(paste0(homedir, "breastepi/R_Analysis/tsne_by_sample.pdf"))
DefaultAssay(inds123) <- "RNA"

DARKBLUE = "#000099"
FeaturePlot(inds123, features = sarsCov2.features,cols = c(LIGHTGRAY, DARKBLUE))
plots <- list()
for (i in 1:length(sarsCov2.features)) {
  if (!(sarsCov2.features[i] %in% rownames(inds123@assays$RNA@data))) {
    print(paste0("gene ", sarsCov2.features[i], " not in RNA@data, continuing..."))
    next
  }
  cur.expressiondata <- inds123@assays$RNA@data[which(rownames(inds123@assays$RNA@data) == sarsCov2.features[i]),]
  cur.plotdata <- cbind(inds123@meta.data, "exprs"=cur.expressiondata)
  cur.plotdata.tsne <- merge(cur.plotdata, inds123@reductions$tsne@cell.embeddings, by='row.names')
  
  #tsne
  cur.plotdata.tsne.order <- cur.plotdata.tsne[order(cur.plotdata.tsne$exprs,decreasing = F),]
  plots[[sarsCov2.features.gene[i]]] <- ggplot(data=cur.plotdata.tsne.order, aes_string(x="tSNE_1", y="tSNE_2")) +
    geom_point(aes(color=exprs), size=0.3) + scale_colour_gradient(low = LIGHTGRAY,high = DARKBLUE) +
    ggtitle(sarsCov2.features.gene[i]) + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), 
          panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18)) 
}
plot_grid(plots[[sarsCov2.features.gene[1]]], plots[[sarsCov2.features.gene[2]]],
          plots[[sarsCov2.features.gene[3]]], plots[[sarsCov2.features.gene[4]]], ncol = 2, nrow = 2)
ggsave(paste0(homedir, "breastepi/R_Analysis/tsne_features_sars.pdf"), width = 8, height = 7.82)
saveRDS(inds123, paste0(homedir, "breastepi/R_Analysis/saves/inds123_processed_732cells.rds"))

####### individuals 4-7 (only GSE113196_RAW samples)
setwd(paste0(homedir, "breastepi/R_Analysis"))
d <- paste0(homedir, "breastepi/GSE113196_RAW/")
samples.files <- list.files(paste0(homedir, "breastepi/GSE113196_RAW/"))
seurat_project <- list()
for (i in 1:length(samples.files)) { #changed to .jfgc
  name = unlist(strsplit(samples.files[i], "_"))[2]
  counts <- read.table(paste0(d, samples.files[i]), header = T, row.names = 1)
  seurat_project[[name]] <- CreateSeuratObject(counts = counts, project = name, min.cells = 3, min.features = 500)
}

#trying to integrate (they didn't in the paper because they used seurat 2.0)
# saveRDS(seurat_project.inds4567, "saves/seurat_project_inds4567_unprocessed.rds")
seurat_project.inds4567 <- seurat_project
for (i in 1:length(seurat_project.inds4567)) {
  seurat_project.inds4567[[i]][["percent.mt"]] <- PercentageFeatureSet(object = seurat_project.inds4567[[i]], pattern = "^MT-")
  seurat_project.inds4567[[i]] <- subset(seurat_project.inds4567[[i]], 
                                         subset = nFeature_RNA > 500 & nFeature_RNA < 6000 &
                                        percent.mt < 10) 
  seurat_project.inds4567[[i]] <- NormalizeData(seurat_project.inds4567[[i]], normalization.method = "LogNormalize", 
                                                scale.factor = 1e4, verbose=F)
  seurat_project.inds4567[[i]] <- FindVariableFeatures(object = seurat_project.inds4567[[i]], 
                                                       selection.method = 'vst', nfeatures = 2000)
}
seurat_project.anchors <- FindIntegrationAnchors(object.list = seurat_project.inds4567, dims = 1:30)
inds4567.int <- IntegrateData(anchorset = seurat_project.anchors, dims = 1:30)
inds4567.int <- ScaleData(inds4567.int, verbose = F, vars.to.regress = c("percent.mt"))
inds4567.int <- RunPCA(object = inds4567.int, npcs=10, verbose=F, seed.use = 777)
inds4567.int <- RunUMAP(object = inds4567.int, reduction = "pca", dims=1:10, seed.use = 123)
inds4567.int <- RunTSNE(object = inds4567.int, reduction = "pca", dims=1:10)
inds4567.int <- FindNeighbors(object = inds4567.int, dims=1:10)
inds4567.int <- FindClusters(object = inds4567.int, resolution = 0.2)

saveRDS(inds4567.int, paste0(homedir, "saves/inds4567_integrated_processed.rds"))
#adjustments to figure -> red and green are actually correct colors
FeaturePlot(inds4567.int, "KRT18") 
ggsave("inds4567_feature_KRT18.pdf")
FeaturePlot(inds4567.int, "KRT14") 
ggsave("inds4567_feature_KRT14.pdf")
FeaturePlot(inds4567.int, "ANKRD30A") 
ggsave("inds4567_feature_ANKRD30A.pdf")
FeaturePlot(inds4567.int, "SLPI") 
ggsave("inds4567_feature_SLPI.pdf")

inds4567.int.markers <- FindAllMarkers(inds4567.int, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
write.table(inds4567.int.markers, paste0(homedir, "breastepi/R_Analysis/inds4567_int_markers.txt"), 
            sep = "\t", quote = F, col.names = NA)
inds4567.int.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
top10genes <- as.vector(top10$gene)

#heatmap
DefaultAssay(inds4567.int) <- "integrated"
DoHeatmap(inds4567.int, features = unique(top10genes)) #do this
height80dpi <- (length(top10genes)*12)*(1/80)+1 #it's a pdf so it ends up being fine (for dpi)
ggsave(paste0(homedir, "breastepi/R_Analysis/inds4567_int_heatmap.pdf"))

# 20200522 -> rename idents and redo all of the sars plots
inds4567.int@meta.data$celltype <- NA
inds4567.int@meta.data$celltype[which(inds4567.int@meta.data$seurat_clusters == 0)] <- "Myofibroblasts"
inds4567.int@meta.data$celltype[which(inds4567.int@meta.data$seurat_clusters == 1)] <- "Luminal Epithelium"
inds4567.int@meta.data$celltype[which(inds4567.int@meta.data$seurat_clusters == 2)] <- "Luminal Epithelium 2"
inds4567.int@meta.data$celltype[which(inds4567.int@meta.data$seurat_clusters == 3)] <- "Myofibroblasts"
inds4567.int@meta.data$celltype[which(inds4567.int@meta.data$seurat_clusters == 4)] <- "Basal"
inds4567.int@meta.data$celltype[which(inds4567.int@meta.data$seurat_clusters == 5)] <- "Luminal Epithelium"
inds4567.int@meta.data$celltype[which(inds4567.int@meta.data$seurat_clusters == 6)] <- "X"

#sars dotplot
sarsCov2.features.gene <- c("TMPRSS2", "CTSB", "CTSL", "ACE2")
# sarsCov2.features <- c("ENSG00000184012", "ENSG00000164733", "ENSG00000135047", "ENSG00000130234")
DotPlot(inds4567.int, assay = "RNA", features=sarsCov2.features.gene, group.by = "celltype", cols = c(LIGHTGRAY, COLORSCHEME2[4])) + #facet_wrap(~ seurat_clusters, ncol=1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1))
# ggsave(paste0(homedir, "breastepi/R_Analysis/sars_cov2/sars_dotplot.pdf"))
ggsave(paste0(homedir, "breastepi/R_Analysis/celltype_rename/sars_cov2/sars_dotplot.pdf"))

#Feature plots and umap to confirm.
UMAPPlot(inds4567.int, group.by = "celltype")
# ggsave(paste0(homedir, "breastepi/R_Analysis/umap_by_cluster.pdf"))
ggsave(paste0(homedir, "breastepi/R_Analysis/celltype_rename/umap_by_cluster.pdf"), width=10, height=10)
UMAPPlot(inds4567.int, group.by = "orig.ident")
# ggsave(paste0(homedir, "breastepi/R_Analysis/umap_by_sample.pdf"))
ggsave(paste0(homedir, "breastepi/R_Analysis/celltype_rename/umap_by_sample.pdf"), width=10, height=10)
DefaultAssay(inds4567.int) <- "RNA"

DARKBLUE = "#000099"
plots <- list()
for (i in 1:length(sarsCov2.features.gene)) {
  if (!(sarsCov2.features.gene[i] %in% rownames(inds4567.int@assays$RNA@data))) {
    print(paste0("gene ", sarsCov2.features.gene[i], " not in RNA@data, continuing..."))
    next
  }
  cur.expressiondata <- inds4567.int@assays$RNA@data[which(rownames(inds4567.int@assays$RNA@data) == sarsCov2.features[i]),]
  cur.plotdata <- cbind(inds4567.int@meta.data, "exprs"=cur.expressiondata)
  cur.plotdata.tsne <- merge(cur.plotdata, inds4567.int@reductions$tsne@cell.embeddings, by='row.names')
  
  #tsne
  cur.plotdata.tsne.order <- cur.plotdata.tsne[order(cur.plotdata.tsne$exprs,decreasing = F),]
  plots[[sarsCov2.features.gene[i]]] <- ggplot(data=cur.plotdata.tsne.order, aes_string(x="tSNE_1", y="tSNE_2")) +
    geom_point(aes(color=exprs), size=0.3) + scale_colour_gradient(low = LIGHTGRAY,high = DARKBLUE) +
    ggtitle(sarsCov2.features.gene[i]) + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), 
          panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18)) 
}
plot_grid(plots[[sarsCov2.features.gene[1]]], plots[[sarsCov2.features.gene[2]]],
          plots[[sarsCov2.features.gene[3]]], plots[[sarsCov2.features.gene[4]]], ncol = 2, nrow = 2)


FeaturePlot(inds4567.int, features = sarsCov2.features.gene,cols = c(LIGHTGRAY, DARKBLUE))
# ggsave(paste0(homedir, "breastepi/R_Analysis/umap_features_sars.pdf"), width = 8, height = 7.82)
ggsave(paste0(homedir, "breastepi/R_Analysis/celltype_rename/umap_features_sars.pdf"), width = 8, height = 7.82)
# saveRDS(inds123, paste0(homedir, "breastepi/R_Analysis/saves/inds123_processed_732cells.rds"))
saveRDS(inds4567.int, paste0(homedir, "breastepi/R_Analysis/saves/inds4567_processed_celltypenamed_24464cells.rds"))








#### end breast epi cell paper ####

####  uterus endo paper ####
uterus.dge <- read.table(paste0(homedir, "uterus_endo/GSM4008666_Adult-Uterus1_dge.txt"), header = T, row.names = 1)
uterus_endo.less <- uterus.dge[,colSums(uterus.dge)<500 & colSums(uterus.dge)> 100]
uterus_endo.more <- uterus.dge[,colSums(uterus.dge)>=500]
uterus_endo.anno <- data.frame(Cell_barcode= colnames(uterus_endo.more),
                                Sample      = replicate("UterusEndo",n=ncol(uterus_endo.more)),
                                Batch       = replicate("UterusEndo",n=ncol(uterus_endo.more)))
uterus_endo.anno[,"Cell_id"]  <- colnames(uterus_endo.more)
uterus_endo.anno[,"Cluster_id"] = replicate("1",n=ncol(uterus_endo.more))
uterus_endo.anno$Development_stage<-"Adult"
uterus_endo.anno$Method<-rep("Microwell-seq")
uterus_endo.anno$Gender<-"Male"
uterus_endo.anno$Source<-rep("HCL")
uterus_endo.anno$Biomaterial<-rep("UterusEndo")
uterus_endo.anno$Name<-rep("UterusEndo")
name<-"UterusEndo"
name_background <- paste(name,"background", sep="_")
name_500more  <-  paste(name,"500more", sep="_")
name_500less  <-  paste(name,"500less", sep="_")
par(mfrow=c(2,1)) #here
hist(colSums(uterus_endo.more),breaks = 200)
abline(v=300)
hist(colSums(uterus_endo.more>0),breaks = 200)
abline(v=300)
summary(colSums(uterus_endo.more>0))

library(Seurat)
seqwell <- CreateSeuratObject(counts = Matrix(as.matrix(uterus_endo.more),sparse=T)
                              ,min.cells = 3,min.features = 300,names.delim = "\\.")  ##no normarlize
seqwell[["percent.mito"]] <- PercentageFeatureSet(object = seqwell, pattern = "^MT-")
seqwell@project.name <- "UterusEndo_20"
#percent.mito in 'example' was 2%, which no cells in this dataset have, using 20% instead (this is most of cells)
seqwell <- subset(seqwell , subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mito < 20)  
seqwell <- NormalizeData(object = seqwell, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
seqwell<- ScaleData(seqwell, vars.to.regress=c("nCount_RNA", "percent.mito"))
par(mfrow=c(1,1))
seqwell<- FindVariableFeatures(object = seqwell, mean.function = ExpMean, dispersion.function = LogVMR 
                            ,x.low.cutoff = 0.01, 
                            x.high.cutoff = 6, y.cutoff = 0.5)
length(seqwell@assays$RNA@var.features)# 2000
var.gene<-seqwell@assays$RNA@var.features
var.gene<-var.gene[!grepl(pattern = "*RPS",x=var.gene)] # 6 removed
var.gene<-var.gene[!grepl(pattern = "*RPL",x=var.gene)] # 8 removed
var.gene<-var.gene[!grepl(pattern = "*MT",x=var.gene)] # 29 removed
length(var.gene) #1957
seqwell <- RunPCA(object = seqwell,  features = var.gene, npcs = 50, do.print = TRUE, 
               pcs.print = 1:5, nfeatures.print = 5)
png(paste0(homedir, "uterus_endo/pcheatmap_1_15.png"), units = "in", res = 300,
    width = 8, height = 15)
PCHeatmap(object = seqwell, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()
png(paste0(homedir, "uterus_endo/pcheatmap_16_30.png"), units = "in", res = 300,
    width = 8, height = 15)
PCHeatmap(object = seqwell, dims = 16:30, cells = 500, balanced = TRUE)
dev.off()
png(paste0(homedir, "uterus_endo/pcheatmap_31_50.png"), units = "in", res = 300,
    width = 8, height = 15)
PCHeatmap(object = seqwell, dims = 31:50, cells = 500, balanced = TRUE)
dev.off()
seqwell <- RunTSNE(object = seqwell, dims = 1:20, do.fast = TRUE)
TSNEPlot(object = seqwell,label = T, pt.size = 1,label.size = 5)
seqwell <- RunUMAP(object = seqwell, dims = 1:20)
UMAPPlot(object = seqwell,label = T, pt.size = 1,label.size = 5)
seqwell <- FindNeighbors(object = seqwell, reduction = "pca", dims = 1:15)
seqwell <- FindClusters(object = seqwell, reduction = "pca", dims = 1:15, 
                     resolution =c(0.6,0.8,1,1.4,2,2.5,4))

#cluster ANNOTATION
celltype_annot <- read.table()
seqwell.mt20@meta.data$celltype <- NA
seqwell.mt20@meta.data$celltype[which(seqwell.mt20@meta.data$RNA_snn_res.0.6 == 0)] <- "Endothelial"
seqwell.mt20@meta.data$celltype[which(seqwell.mt20@meta.data$RNA_snn_res.0.6 == 1)] <- "Endothelial"
seqwell.mt20@meta.data$celltype[which(seqwell.mt20@meta.data$RNA_snn_res.0.6 == 2)] <- "Endothelial"
seqwell.mt20@meta.data$celltype[which(seqwell.mt20@meta.data$RNA_snn_res.0.6 == 3)] <- "Smooth Muscle"
seqwell.mt20@meta.data$celltype[which(seqwell.mt20@meta.data$RNA_snn_res.0.6 == 4)] <- "Smooth Muscle"
seqwell.mt20@meta.data$celltype[which(seqwell.mt20@meta.data$RNA_snn_res.0.6 == 5)] <- "Stromal"
seqwell.mt20@meta.data$celltype[which(seqwell.mt20@meta.data$RNA_snn_res.0.6 == 6)] <- "Fibroblast"
seqwell.mt20@meta.data$celltype[which(seqwell.mt20@meta.data$RNA_snn_res.0.6 == 7)] <- "Endothelial"
seqwell.mt20@meta.data$celltype[which(seqwell.mt20@meta.data$RNA_snn_res.0.6 == 8)] <- "Luminal Epithelium"
seqwell.mt20@meta.data$celltype[which(seqwell.mt20@meta.data$RNA_snn_res.0.6 == 9)] <- "Macrophages"

sarsCov2.features.gene <- c("TMPRSS2", "CTSB", "CTSL", "ACE2")
DotPlot(seqwell.mt20, assay = "RNA", features=sarsCov2.features.gene, group.by = "celltype", cols = c(LIGHTGRAY, COLORSCHEME2[4])) + #facet_wrap(~ seurat_clusters, ncol=1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1))
ggsave(paste0(homedir, "uterus_endo/mt_lessthan20/after_annotation/dotplot_byCluster.pdf"), width=10, height=10)

#Feature plots and umap to confirm.
UMAPPlot(seqwell.mt20, group.by = "celltype")
ggsave(paste0(homedir, "uterus_endo/mt_lessthan20/after_annotation/umap_byCluster.pdf"), width=10, height=10)

DefaultAssay(seqwell.mt20) <- "RNA"
plots <- list()
for (i in 1:length(sarsCov2.features.gene)) {
  if (!(sarsCov2.features.gene[i] %in% rownames(seqwell.mt20@assays$RNA@data))) {
    print(paste0("gene ", sarsCov2.features.gene[i], " not in RNA@data, continuing..."))
    next
  }
  cur.expressiondata <- seqwell.mt20@assays$RNA@data[which(rownames(seqwell.mt20@assays$RNA@data) == sarsCov2.features.gene[i]),]
  cur.plotdata <- cbind(seqwell.mt20@meta.data, "exprs"=cur.expressiondata)
  cur.plotdata.umap <- merge(cur.plotdata, seqwell.mt20@reductions$umap@cell.embeddings, by='row.names')
  
  #umap
  cur.plotdata.umap.order <- cur.plotdata.umap[order(cur.plotdata.umap$exprs,decreasing = F),]
  plots[[sarsCov2.features.gene[i]]] <- ggplot(data=cur.plotdata.umap, aes_string(x="UMAP_1", y="UMAP_2")) +
    geom_point(aes(color=exprs), size=0.3) + scale_colour_gradient(low = LIGHTGRAY,high = DARKBLUE) +
    ggtitle(sarsCov2.features.gene[i]) + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), 
          panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18)) 
}
plot_grid(plots[[sarsCov2.features.gene[1]]], plots[[sarsCov2.features.gene[2]]],
          plots[[sarsCov2.features.gene[3]]], plots[[sarsCov2.features.gene[4]]], ncol = 2, nrow = 2)
ggsave(paste0(homedir, "uterus_endo/mt_lessthan20/after_annotation/sars_features.pdf"), width=8, height=6)
clusterDE(seqwell.mt20, "fallopian_tube_benign", cluster_col = "celltype")


#before annotation
seqwell.mt20@meta.data$condition <- "UterusEndometrium"
allUMAPs(seqwell.mt20)
sarsCov2.features.gene <- c("TMPRSS2", "CTSB", "CTSL", "ACE2")
# sarsCov2.features <- c("ENSG00000184012", "ENSG00000164733", "ENSG00000135047", "ENSG00000130234")
DotPlot(seqwell.mt20, assay = "RNA", features=sarsCov2.features.gene, group.by = "RNA_snn_res.0.6", cols = c(LIGHTGRAY, COLORSCHEME2[4])) + #facet_wrap(~ seurat_clusters, ncol=1) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1))
ggsave(paste0(homedir, "uterus_endo/mt_lessthan20/dotplot_bycelltype_mt20.pdf"), width=10, height=10)

#Feature plots and umap to confirm.
UMAPPlot(seqwell.mt20, group.by = "RNA_snn_res.0.6", label=T)
if (seqwell.mt20@project.name == "UterusEndo") {
  ggsave(paste0(homedir, "uterus_endo/mt_lessthan7/umap_bycelltype_mt7.pdf"), width=10, height=10)
} else {
  ggsave(paste0(homedir, "uterus_endo/mt_lessthan20/umap_bycelltype_mt20.pdf"), width=10, height=10)
}
DefaultAssay(seqwell.mt20) <- "RNA"
plots <- list()
for (i in 1:length(sarsCov2.features.gene)) {
  if (!(sarsCov2.features.gene[i] %in% rownames(seqwell.mt20@assays$RNA@data))) {
    print(paste0("gene ", sarsCov2.features.gene[i], " not in RNA@data, continuing..."))
    next
  }
  cur.expressiondata <- seqwell.mt20@assays$RNA@data[which(rownames(seqwell.mt20@assays$RNA@data) == sarsCov2.features.gene[i]),]
  cur.plotdata <- cbind(seqwell.mt20@meta.data, "exprs"=cur.expressiondata)
  cur.plotdata.tsne <- merge(cur.plotdata, seqwell.mt20@reductions$tsne@cell.embeddings, by='row.names')
  
  #tsne
  cur.plotdata.tsne.order <- cur.plotdata.tsne[order(cur.plotdata.tsne$exprs,decreasing = F),]
  plots[[sarsCov2.features.gene[i]]] <- ggplot(data=cur.plotdata.tsne.order, aes_string(x="tSNE_1", y="tSNE_2")) +
    geom_point(aes(color=exprs), size=0.3) + scale_colour_gradient(low = LIGHTGRAY,high = DARKBLUE) +
    ggtitle(sarsCov2.features.gene[i]) + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), 
          panel.grid = element_blank(), plot.background = element_rect(fill="white"), text = element_text(face="bold",size=18)) 
}
plot_grid(plots[[sarsCov2.features.gene[1]]], plots[[sarsCov2.features.gene[2]]],
          plots[[sarsCov2.features.gene[3]]], plots[[sarsCov2.features.gene[4]]], ncol = 2, nrow = 2)

FeaturePlot(seqwell.mt20, features = sarsCov2.features.gene,cols = c(LIGHTGRAY, DARKBLUE))
ggsave(paste0(homedir, "uterus_endo/mt_lessthan20/sars_cov2_features.mt20.pdf"), width=8, height=7)
saveRDS(seqwell.mt20, paste0(homedir, "uterus_endo/mt_lessthan20/uterus_endo.mt20.rds"))
clusterDE(seqwell.mt20, "uterus_endo_mt20", cluster_col = "RNA_snn_res.0.6")

#### end uterust endo paper ####


#### myometrium ####
myo <- readRDS(paste0(homedir, "saves/myo_all5samples_unfiltered.rds"))
#separate the object into the 5 
myo.samples <- list()
myo.samples[["myometrium061919"]] <- subset(myo, subset=orig.ident == "myometrium061919")
myo.samples[["myometrium11911"]] <- subset(myo, subset=orig.ident == "myometrium11911")
myo.samples[["myometrium55_11564"]] <- subset(myo, subset=orig.ident == "myometrium55_11564")
myo.samples[["myometrium55_12640"]] <- subset(myo, subset=orig.ident == "myometrium55_12640")
myo.samples[["myometrium55_12745"]] <- subset(myo, subset=orig.ident == "myometrium55_12745")

myo.samples <- lapply(myo.samples, pp,
       HEMOGENES=HEMOGENES.HSA, subset=TRUE, SCT = TRUE, VAR.FEATURES=2000, 
          MAX_FEATURES_PER_CELL=2500, MIN_FEATURES_PER_CELL=200, MAX_PERCENT_MT=7
       )
myo.samples.feats <- SelectIntegrationFeatures(object.list = myo.samples, nfeatures = 2000)
myo.samples.feats <- myo.samples.feats[grep("^M[Tt]-|^H[Bb][BbAa][12]*", myo.samples.feats, invert = TRUE)]
myo.samples <- PrepSCTIntegration(object.list = myo.samples, anchor.features = myo.samples.feats, verbose = FALSE)
myo.samples.anchors <- FindIntegrationAnchors(object.list = myo.samples, normalization.method = "SCT", 
                                           anchor.features = myo.samples.feats, verbose = FALSE)
myo.sct.int <- IntegrateData(anchorset = myo.samples.anchors, normalization.method = "SCT", verbose = FALSE)
myo.sct.int <- ScaleData(myo.sct.int)
myo.sct.int <- RunPCA(object = myo.sct.int, features = VariableFeatures(object = myo.sct.int))
myo.sct.int <- RunUMAP(object = myo.sct.int, dims=1:20)
myo.sct.int <- FindNeighbors(object = myo.sct.int, dims=1:20)
myo.sct.int <- FindClusters(object = myo.sct.int, dims=1:20, resolution = 0.4)
saveRDS(myo.sct.int, paste0(paste0(homedir, "",
        "myometrium_from_bigproject/R_Analysis/saves/myo_sct_int.rds")))
setwd(paste0(paste0(homedir, "",
      "myometrium_from_bigproject/R_Analysis/")))
myo.sct.int@meta.data$condition <- "myometrium"
clusterDE(myo.sct.int, "myo_sct_int", cluster_col=NULL, assay="RNA")

# remove cluster 11 -> epithelial
myo.sct.int.backup <- myo.sct.int
myo.sct.int <- subset(myo.sct.int, subset = seurat_clusters != 11)
myo.sct.int <- FindVariableFeatures(myo.sct.int, nfeatures = 2000)
myo.sct.int <- RunPCA(object = myo.sct.int, features = VariableFeatures(object = myo.sct.int))
myo.sct.int <- RunUMAP(object = myo.sct.int, dims=1:20)
myo.sct.int <- FindNeighbors(object = myo.sct.int, dims=1:20)
myo.sct.int <- FindClusters(object = myo.sct.int, dims=1:20, resolution = 0.4)
allUMAPs(myo.sct.int)
clusterDE(myo.sct.int, "myo_sct_int", cluster_col="seurat_clusters", assay="integrated")

# add celltype names and remake all figures
myo.sct.int@meta.data$celltype <- NA
myo.sct.int@meta.data$celltype[which(myo.sct.int@meta.data$seurat_clusters == 0)] <- "Smooth Muscle"
myo.sct.int@meta.data$celltype[which(myo.sct.int@meta.data$seurat_clusters == 1)] <- "Endothelial"
myo.sct.int@meta.data$celltype[which(myo.sct.int@meta.data$seurat_clusters == 2)] <- "Endothelial"
myo.sct.int@meta.data$celltype[which(myo.sct.int@meta.data$seurat_clusters == 3)] <- "Endothelial"
myo.sct.int@meta.data$celltype[which(myo.sct.int@meta.data$seurat_clusters == 4)] <- "NK"
myo.sct.int@meta.data$celltype[which(myo.sct.int@meta.data$seurat_clusters == 5)] <- "Smooth Muscle"
myo.sct.int@meta.data$celltype[which(myo.sct.int@meta.data$seurat_clusters == 6)] <- "NK"
myo.sct.int@meta.data$celltype[which(myo.sct.int@meta.data$seurat_clusters == 7)] <- "Myeloid"
myo.sct.int@meta.data$celltype[which(myo.sct.int@meta.data$seurat_clusters == 8)] <- "T"
myo.sct.int@meta.data$celltype[which(myo.sct.int@meta.data$seurat_clusters == 9)] <- "Fibroblast"
myo.sct.int@meta.data$celltype[which(myo.sct.int@meta.data$seurat_clusters == 10)] <- "Lymphatic Endothelial"
myo.sct.int@meta.data$celltype[which(myo.sct.int@meta.data$seurat_clusters == 11)] <- "Myeloid"
myo.sct.int@meta.data$celltype[which(myo.sct.int@meta.data$seurat_clusters == 12)] <- "Smooth Muscle"
myo.sct.int@meta.data$celltype[which(myo.sct.int@meta.data$seurat_clusters == 13)] <- "Myeloid"

Idents(myo.sct.int) <- myo.sct.int@meta.data$celltype
DefaultAssay(myo.sct.int) <- "RNA"
x <- UMAPPlot(myo.sct.int, pt.size = 0.3, group.by="seurat_clusters", label=TRUE) +
  theme(panel.grid = element_blank(), 
        plot.background = element_rect(fill="white"), 
        legend.text = element_text(face="bold",size=18)) +
  guides(color=guide_legend(ncol=1, override.aes = list(size=10)))
y <- ggplot_build(x)
y$data[[2]]$fontface <- 2
y$data[[2]]$size <- 4
z <- ggplot_gtable(y)
pdf(paste0(paste0(homedir, "myometrium_from_bigproject/R_Analysis/",
           "myo_sct_int_celltypeumap.pdf")))
plot(z) 
dev.off()

UMAPPlot(myo.sct.int, pt.size = 0.3, group.by="celltype") +
  theme(panel.grid = element_blank(), 
        plot.background = element_rect(fill="white"), 
        legend.text = element_text(face="bold",size=18)) +
  guides(color=guide_legend(ncol=1, override.aes = list(size=10)))
ggsave(paste0(paste0(homedir, "myometrium_from_bigproject/R_Analysis/",
           "myo_sct_int_celltypeumap.pdf")), width = 10, height=6.7)

myo.sct.int <- NormalizeData(myo.sct.int, normalization.method = "LogNormalize", scale.factor = 1e4, verbose=F)
myo.sct.int <- ScaleData(myo.sct.int, assay = "RNA", vars.to.regress = c("percent.mt", "percent.hbb"))
clusterDE(myo.sct.int, "myo_sct_int", cluster_col="celltype", assay="RNA")
myo.sct.int.markers <- read.table(paste0(homedir, "myometrium_from_bigproject/R_Analysis/markers_myo_sct_int.txt"), 
  header=T, sep = "\t", row.names = 1)
top10 <- myo.sct.int.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top2 <- myo.sct.int.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
DoHeatmap(myo.sct.int, features = as.vector(top2$gene)) 
height80dpi <- (length(as.vector(top2$gene))*12)*(1/80)+1
ggsave(paste0(homedir, "myometrium_from_bigproject/R_Analysis/heatmap_celltype_top2.pdf"), height=height80dpi, width=8)
sarsCov2Coexpression(myo.sct.int, paste0(homedir, "myometrium_from_bigproject/R_Analysis"))
FeaturePlot(myo.sct.int, features = sarsCov2.features.gene,cols = c(LIGHTGRAY, DARKBLUE))
ggsave(paste0(homedir, "myometrium_from_bigproject/R_Analysis/sars_featureplot.pdf"))

#### end myometrium ####



#### sars cov2 coexpression ####
all.myo <- readRDS(paste0(homedir, "myometrium_from_bigproject/R_Analysis/saves/all_myometrium4sarscov2Paper.rds"))
inds4567.int <- readRDS(paste0(homedir, "breastepi/R_Analysis/saves/inds4567_processed_celltypenamed_24464cells.rds"))
ovary.int <- readRDS(paste0(homedir, "ovary/R_Analysis/saves/ovary_unsorted_int_with_celltype_names.rds"))
Idents(inds4567.int) <- inds4567.int@meta.data$celltype
Idents(ft) <- ft@meta.data$celltype
Idents(seqwell.mt20) <- seqwell.mt20@meta.data$celltype
Idents(ovary.int) <- ovary.int@meta.data$celltype
Idents(all.myo) <- all.myo@meta.data$celltype
DefaultAssay(inds4567.int) <- "RNA"
DefaultAssay(ft) <- "RNA"
DefaultAssay(seqwell.mt20) <- "RNA"
DefaultAssay(ovary.int) <- "RNA"
DefaultAssay(all.myo) <- "RNA"

sarsCov2Coexpression(ft, paste0(homedir, "fallopian/R_Analysis/sarscov2"))
sarsCov2Coexpression(seqwell.mt20,paste0(paste0(homedir,"uterus_endo/mt_lessthan20/after_annotation")))
sarsCov2Coexpression(inds4567.int, paste0(paste0(homedir, "breastepi/R_Analysis/celltype_rename")))
sarsCov2Coexpression(ovary.int, paste0(homedir, "ovary/R_Analysis/sars_cov2"))
sarsCov2Coexpression(all.myo, paste0(homedir, "myometrium_from_bigproject/R_Analysis"))

#### final umaps ####
seqwell.mt20@meta.data$seurat_clusters <- seqwell.mt20@meta.data$RNA_snn_res.0.6
oblist <- list(ft, seqwell.mt20, inds4567.int, ovary.int, all.myo)
names(oblist) <- c("ft", "seqwell.mt20", "inds4567.int", "ovary.int", "all.myo")
for (i in 1:length(oblist)) {

  x <- UMAPPlot(oblist[[i]], pt.size = 0.3, group.by="seurat_clusters", label=TRUE) +
    theme(panel.grid = element_blank(), 
          plot.background = element_rect(fill="white"), 
          legend.text = element_text(face="bold",size=18)) +
    guides(color=guide_legend(ncol=1, override.aes = list(size=10)))
  y <- ggplot_build(x)
  y$data[[2]]$fontface <- 2
  y$data[[2]]$size <- 4
  z <- ggplot_gtable(y)
  pdf(paste0(paste0(homedir, "finalUmaps/", names(oblist)[i], "_betterumap.pdf")))
  plot(z) 
  dev.off()
}


# human data from https://doi.org/10.1016/j.cell.2020.04.035 ----
# download is here: https://drive.google.com/drive/folders/1bxCIqNeZ7wLuVOT16gphwj98_cc9KhrV
setwd("")
homedir <- getwd()
raw <- read.table("Human_lung_epithelial_cell_raw_counts.txt", header=T, row.names = 1, sep = "\t")
meta <- read.csv("Human_lung_epithelial_metadata.csv", header=T, row.names = 1)
lung <- CreateSeuratObject(as.matrix(raw))
length(intersect(rownames(meta), colnames(lung@assays$RNA@counts)))
meta <- meta[rownames(lung@meta.data),]
lung@meta.data <- merge(lung@meta.data, meta, by="row.names")
rownames(lung@meta.data) <- lung@meta.data[,1]
HEMOGENES=HEMOGENES.HSA
lung <- pp(lung, 
           HEMOGENES=HEMOGENES.HSA, 
           subset=TRUE, 
           SCT = FALSE, 
           VAR.FEATURES=2000, 
           MAX_FEATURES_PER_CELL=2500, 
           MIN_FEATURES_PER_CELL=200)

#umap
x <- UMAPPlot(lung, pt.size = 1, group.by="seurat_clusters", label=TRUE) +
  theme(panel.grid = element_blank(), 
        plot.background = element_rect(fill="white"), 
        legend.text = element_text(face="bold",size=18)) +
  guides(color=guide_legend(ncol=1, override.aes = list(size=10)))
y <- ggplot_build(x)
y$data[[2]]$fontface <- 2
y$data[[2]]$size <- 4
z <- ggplot_gtable(y)
pdf(paste0(homedir, "/betterumap.pdf"))
plot(z)
dev.off()

DefaultAssay(lung) #RNA

#coexpression dotplot
sarsCov2Coexpression(lung, paste0(homedir, "/"))

#features to show why ace2 + tmprss2 aren't in the same cells (it's because there aren't many cells)
FeaturePlot(lung, features=c("ACE2", "TMPRSS2"), blend=T)
ggsave("lung_ace2_tmprss2_blend.pdf", width=16,height=4)
FeaturePlot(lung, features = sarsCov2.features.gene,cols = c(LIGHTGRAY, DARKBLUE))
ggsave(paste0(homedir, "/sars_featureplot.pdf"), height=10,width=10)

#cluster DE
lung.markers <- FindAllMarkers(lung, only.pos = T)
write.table(lung.markers, "lung_cluster_markers.txt", sep = "\t", quote = F, col.names = NA)

#heatmap
lung.markers$rank <- (lung.markers$pct.1 * lung.markers$avg_logFC) / lung.markers$pct.2
top5 <- lung.markers %>% group_by(cluster) %>% top_n(5, rank)
height80dpi <- (length(as.vector(top5$gene))*12)*(1/80)+1
pdf("heatmap.pdf", width=5.5,height=height80dpi)
print(DoHeatmap(lung, features = as.vector(top5$gene)))
dev.off()




