ibrary(Giotto)
library(tidyverse)
options(stringsAsFactors = F)
library(dplyr)
library(limma)
library(NMF)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(CellChat)
#sc intercellular
expdata <- read.csv("C:/Users/sun/Desktop/tutorial/data/example_data/demo1/predata/output_1/demo1_new_sc_data.csv",row.names = 1)
anndata <- read.csv("C:/Users/sun/Desktop/tutorial/data/example_data/demo1/predata/output_1/demo1_new_sc_celltype.csv",row.names = 1)
my_working_dir = 'E:/BScSpaCall-R/demo1'
setwd(my_working_dir)
write.table(expdata,'E:/BScSpaCall-R/demo1/sc_expression.txt')
expr_path <- 'E:/BScSpaCall-R/demo1/sc_expression.txt'
python_path <- 'C:/Users/sun/AppData/Local/r-miniconda/envs/giotto_env/python.exe'
instrs <- createGiottoInstructions(save_plot = TRUE, show_plot = FALSE, save_dir = my_working_dir, python_path = python_path)
scRNA_giotto <- createGiottoObject(raw_exprs = expdata, instructions = instrs)
cellmeta <- data.frame(ID = anndata$Cell , cell_types = anndata$Cell_type)
scRNA_giotto <- addCellMetadata(scRNA_giotto, new_metadata = cellmeta, by_column = TRUE, column_cell_ID = 'ID')
scRNA_giotto <- normalizeGiotto(gobject = scRNA_giotto, scalefactor = 6000, verbose = TRUE)
scRNA_giotto <- addStatistics(gobject = scRNA_giotto)
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")  
LR_data <- data.frame(humanLigand = CellChatDB.use[["interaction"]][["ligand"]] , humanReceptor= CellChatDB.use[["interaction"]][["receptor"]])
all_scores <- c()
feat_ID <- scRNA_giotto@feat_metadata[["cell"]][["rna"]]@metaDT[["feat_ID"]]
library(data.table)
LR_data <- as.data.table(LR_data)
LR_data[, ligand_det := ifelse(humanLigand %in% feat_ID, TRUE, FALSE)]
LR_data[, receptor_det := ifelse(humanReceptor %in% feat_ID, T, F)]
LR_data_det = LR_data[ligand_det == T & receptor_det == T]
select_ligands = LR_data_det$humanLigand
select_receptors = LR_data_det$humanReceptor
expr_only_scores = exprCellCellcom(scRNA_giotto,
                                   cluster_column = 'cell_types',
                                   random_iter = 500,
                                   gene_set_1 = select_ligands,
                                   gene_set_2 = select_receptors )

#st intercellular
library(Giotto)
library(tidyverse)
options(stringsAsFactors = F)
library(dplyr)
expdata <- read.csv("C:/Users/sun/Desktop/tutorial/data/example_data/demo1/result/data_demo1_5.csv",row.names = 1)
anndata <- read.csv("C:/Users/sun/Desktop/tutorial/data/example_data/demo1/result/meta_demo1_5.csv",row.names = 1)
my_working_dir = 'E:/BScSpaCall-R/demo1'
setwd(my_working_dir)
write.table(expdata,'E:/BScSpaCall-R/demo1/sc_expression.txt')
expr_path <- 'E:/BScSpaCall-R/demo1/sc_expression.txt'
python_path <- 'C:/Users/sun/AppData/Local/r-miniconda/envs/giotto_env/python.exe'

anndata <- anndata %>% mutate(View = 0) %>% rename(c("ID" = "Cell","X" = "Cell_xcoord", "Y" = "Cell_ycoord", "cell_types" = "Cell_type", "FOV" = "View"))
starmap_meta <- anndata
instrs = createGiottoInstructions(save_plot = TRUE,
                                  show_plot = FALSE,
                                  save_dir = my_working_dir,
                                  python_path = python_path)
my_offset_file = data.table::data.table(field = c(0),
                                        x_offset = c(0),
                                        y_offset = c(0))
starmap_meta_location = data.frame(starmap_meta[,c(1,6,7,8)])
location_file = starmap_meta_location
stitch_file = stitchFieldCoordinates(location_file = starmap_meta_location,
                                     offset_file = my_offset_file,
                                     cumulate_offset_x = T,
                                     cumulate_offset_y = F,
                                     field_col = 'FOV',
                                     reverse_final_x = F,
                                     reverse_final_y = T)
stitch_file    = stitch_file[,.(ID, X_final, Y_final)]
my_offset_file = my_offset_file[,.(field, x_offset_final, y_offset_final)]
cellmeta <- starmap_meta %>% select(ID,FOV,cell_types)
starmap <- createGiottoObject(raw_exprs = expr_path,
                              spatial_locs = stitch_file,
                              offset_file = my_offset_file,
                              instructions = instrs)
starmap <- addCellMetadata(starmap,
                           new_metadata = cellmeta,
                           by_column = T,
                           column_cell_ID = 'ID')
starmap <- normalizeGiotto(gobject = starmap, scalefactor = 6000, verbose = T)
starmap <- addStatistics(gobject = starmap)
total_expr <- starmap@cell_metadata[["cell"]][["rna"]]@metaDT[["total_expr"]]
nr_feats <- starmap@cell_metadata[["cell"]][["rna"]]@metaDT[["nr_feats"]]
starmap <- adjustGiottoMatrix(gobject = starmap, expression_values = c('normalized'),
                              batch_columns = NULL, covariate_columns = c('nr_feats', 'total_expr'),
                              return_gobject = TRUE,
                              update_slot = c('custom'))

starmap <- createSpatialNetwork(gobject = starmap, method = 'kNN', k = 5, name = 'spatial_network')
p1 <- spatPlot(gobject = starmap, show_network = T,
               network_color = "#7CFC00" ,spatial_network_name = 'spatial_network' ,
               point_size = 8, cell_color = 'cell_types') + ggtitle("spatial_network" )

CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB 
LR_data <- data.frame(humanLigand = CellChatDB.use[["interaction"]][["ligand"]] , humanReceptor= CellChatDB.use[["interaction"]][["receptor"]])
spatial_all_scores <- c()
feat_ID <- starmap@feat_metadata[["cell"]][["rna"]]@metaDT[["feat_ID"]]
library(data.table)
LR_data <- as.data.table(LR_data)
LR_data[, ligand_det := ifelse(humanLigand %in% feat_ID, TRUE, FALSE)]
LR_data[, receptor_det := ifelse(humanReceptor %in% feat_ID, T, F)]
LR_data_det = LR_data[ligand_det == T & receptor_det == T]
select_ligands = LR_data_det$humanLigand
select_receptors = LR_data_det$humanReceptor

spatial_all_scores = spatCellCellcom(starmap,
                                     spatial_network_name = 'spatial_network',
                                     cluster_column = 'cell_types',
                                     random_iter = 300,
                                     feat_set_1 = select_ligands,
                                     feat_set_2 = select_receptors,
                                     adjust_method = 'fdr',
                                     do_parallel = T,
                                     cores = 4,
                                     verbose = 'a little')
#intracellular
mt <- CreateNichConObject(data, min.feature = 3,
                          source = "TPM",
                          scale.factor = 10^6,
                          Org = "Homo sapiens",
                          project = "Microenvironment")
mt@meta.data[["celltype"]]<-meta[,2]
mt <- TransCommuProfile(object = mt,
                        pValueCor = 0.05,
                        CorValue = 0.1,
                        topTargetCor=1,
                        p.adjust = 0.05,
                        use.type="median", 
                        probs = 0.9,
                        method="weighted",
                        IS_core = TRUE,
                        Org = 'Homo sapiens')
ligand <- expr_only_scores$ligand
lig_cell_type <- expr_only_scores$lig_cell_type
lig_expr <- expr_only_scores$lig_expr
unique_combinations1 <- interaction(ligand, lig_cell_type, sep = " <-> ")
mean_expr1 <- aggregate(lig_expr ~ unique_combinations1, data = data.frame(unique_combinations1, lig_expr), FUN = mean)
mean_expr1$lig_cell_type <- sub(" <-> .*", "", mean_expr1$unique_combinations1)
mean_expr1$ligand <- sub(".* <-> ", "", mean_expr1$unique_combinations1)
mean_expr1 <- mean_expr1[, -1]
library(tidyr)
new_df1 <- spread(mean_expr1, key = ligand, value = lig_expr)
rownames(new_df1) <- new_df1$lig_cell_type
new_df1 <- new_df1[, -1]
new_df_transposed1 <- new_df1
new_df_transposed1<-as.data.frame(new_df_transposed1)
receptor <- expr_only_scores$receptor
rec_cell_type <- expr_only_scores$rec_cell_type
rec_expr <- expr_only_scores$rec_expr
unique_combinations2 <- interaction(receptor, rec_cell_type, sep = " <-> ")
mean_expr2 <- aggregate(rec_expr ~ unique_combinations2, data = data.frame(unique_combinations2, rec_expr), FUN = mean)
mean_expr2$receptor <- sub(" <-> .*", "", mean_expr2$unique_combinations2)
mean_expr2$rec_cell_type <- sub(".* <-> ", "", mean_expr2$unique_combinations2)
mean_expr2 <- mean_expr2[, -1]
library(tidyr)
new_df2 <- spread(mean_expr2, key = receptor, value = rec_expr)
rownames(new_df2) <- new_df2[,1]
new_df2 <- new_df2[, -1]
new_df_transposed2 <- t(new_df2)
new_df_transposed2<-as.data.frame(new_df_transposed2)
LR_comb <- expr_only_scores$LR_comb
LR_expr <- expr_only_scores$LR_expr
LR_cell_comb <- expr_only_scores$LR_cell_comb
LR_cell_comb <- gsub("--", "-", LR_cell_comb)
unique_combinations <- interaction(LR_comb, LR_cell_comb, sep = " <-> ")
mean_expr <- aggregate(LR_expr ~ unique_combinations, data = data.frame(unique_combinations, LR_expr), FUN = mean)
mean_expr$LR_comb <- sub(" <-> .*", "", mean_expr$unique_combinations)
mean_expr$LR_cell_comb <- sub(".* <-> ", "", mean_expr$unique_combinations)
mean_expr <- mean_expr[, -1]
library(tidyr)
new_df <- spread(mean_expr, key = LR_comb, value = LR_expr)
rownames(new_df) <- new_df$LR_cell_comb
new_df <- new_df[, -1]
new_df_transposed <- t(new_df)
new_df_transposed<-as.data.frame(new_df_transposed)
mt@data[["expr_l_r_log2_scale"]] <- new_df_transposed
mt@data[["expr_l_r"]] <- new_df_transposed
mt@data[["expr_l_r_log2"]] <- new_df_transposed
mt@data[["softmax_ligand"]]<-new_df_transposed1
mt@data[["softmax_receptor"]]<-new_df_transposed2
n <- mt@data$expr_l_r_log2_scale
cell_color <- data.frame(color=c('#62b58f','#bc95df', '#67cdf2', '#ffc533'), stringsAsFactors = FALSE)
rownames(cell_color) <- c("T cell", "B cell", "Monocyte", "Dendritic cell")

mt <- LR2TF(object = mt1, sender_cell="T cell", recevier_cell="B cell",
            slot="expr_l_r_log2_scale", org="Homo sapiens")

head(mt@reductions$sankey)
if(!require(networkD3)){
  BiocManager::install("networkD3")
}
library(networkD3)
cell_color <- data.frame(color = c('#62b58f', '#bc95df', '#67cdf2', '#ffc533'), 
                         stringsAsFactors = FALSE)
rownames(cell_color) <- c("T cell", "B cell", "Monocyte", "Dendritic cell")
sank <- LRT.Dimplot(mt, fontSize = 18, nodeWidth = 50, height = 800, width = 1200, 		 
                    sinksRight=FALSE, DIY.color = FALSE)
sank
networkD3::saveNetwork(sank, "~/Dendritic-Monocyte_full1.html")


