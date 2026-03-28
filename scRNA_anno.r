library(Seurat)
library(dplyr)
library(ggplot2)

sc_input <- readRDS('/data/work/MP/scRNA/scRNA_total.rds')

Idents(sc_input) <- sc_input$`seurat_clusters`
main_type_anno <- c("0"="Erythrocytes",
"1"="Endothelial cells",
"2"="Fibroblasts",
"3"="Endothelial cells",
"4"="Fibroblasts",
"5"="T cells",
"6"="Pericytes",
"7"="Epithelial cells",
"8"="B cells",
"9"="Monocytes",
"10"="Neutrophils",
"11"="T cells",
"12"="Monocytes",
"13"="Endothelial cells",
"14"="Platelets",
"15"="Macrophages",
"16"="Epithelial cells",
"17"="T cells",
"18"="Epithelial cells",
"19"="Epithelial cells",
"20"="Endothelial cells",
"21"="Epithelial cells",
"22"="Epithelial cells",
"23"="Endothelial cells")
                   

#names(main_type_anno) 
sc_input <- RenameIdents(sc_input, main_type_anno)
sc_input$total_cell_type_res_0.5 <- Idents(sc_input)

cell_colors <- c(
'Erythrocytes' = '#d2e0ac',    
'Endothelial cells' = '#74a893',
'Fibroblasts' = '#7587b1',   
'T cells' = '#edeaa4',  
'Pericytes' = '#e3d1db',
'Epithelial cells' = '#F19294',
'B cells' = '#C0937E',
'Monocytes' = '#684797',
'Neutrophils' = '#BDA7CB',
 'Platelets' = "#5EC0AD",   
'Macrophages' = '#e0bc58'
)
DimPlot(sc_input, reduction = "umap", label = T,group.by = "total_cell_type_res_0.5",cols =cell_colors)
ggsave("/data/work/MP/scRNA/total_cell_type_res_0.5_umap.pdf",bg = "transparent", width =9.5, height =8)

mainmarkers <- c("PECAM1", "VWF", "DCN", "MGP", "IGFBP5", "CREG1", "SLC4A1", "IL7R", "BCL11B",
"EPCAM", "KRT18", "PDGFRB", "KCNJ8", "RGS5", "C1QA", "C1QB", "CRIP1", "CSF1R",
"EBF1", "CD79B", "CD74", "CSF3R", "G0S2", "MMP9", "CD9", "F13A1", "GP1BA", "GP1BB",
"ITGB3", "MPL", "STMN1", "MKI67", "TOP2A", "CENPF")

sc_input$total_cell_type_res_0.5 <- factor(sc_input$total_cell_type_res_0.5, 
                                                   levels = c(
'Endothelial cells',
'Fibroblasts',                                                       
'Erythrocytes',
'T cells',
"Epithelial cells",
'Pericytes',
'Monocytes',
'B cells',
'Neutrophils',
'Platelets',
'Macrophages'
))

p <- DotPlot(sc_input, features = mainmarkers, group.by = "total_cell_type_res_0.5", scale = TRUE) +
  theme(axis.text = element_text(size = 20))  +
  scale_color_gradientn(colors = c('#0571b0', '#92c5de', '#f7f7f7', '#f4a582', '#ca0020')) +
  theme(text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.5, "inches"),
        axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1)) 

dat <- p$data

ggplot(dat, aes(features.plot, id,size=pct.exp, fill=avg.exp.scaled)) + 
  geom_point(shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill=NA))) + 
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.major.y = element_line(color = "grey80"),
    axis.title = element_blank(),
    axis.text.y = element_text(color='black',size=12),
    axis.text.x = element_text(color='black',size=12, angle = 90, hjust = 1, vjust = 0.5))+
  scale_fill_gradientn(colours = c('#7CACCE', '#FFFFFF', '#EA85A8'))
ggsave("/data/work/MP/scRNA/total_cell_type_res_0.5.pdf", width = 11, height = 8)