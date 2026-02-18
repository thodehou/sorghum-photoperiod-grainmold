source("functions/QTL_colocalization.R")

qtl_data <- read.csv("data/QTL-colocalization.csv")
qtl_data$Trait <- as.factor(qtl_data$Trait)
qtl_data$Chr_string <- factor(qtl_data$Chr_string,
                              levels = paste0("Chr", sprintf("%02d", 1:10)))
qtl_data$Start_QTL_v3 <- as.numeric(qtl_data$Start_QTL_v3)
qtl_data$End_QTL_v3   <- as.numeric(qtl_data$End_QTL_v3)

mycol <- c("FLA"="royalblue","Mold"="firebrick","YLD"="gold")

ma_genes <- data.frame(
  Gene = c("Ma1","Ma2","Ma3","Ma5","Ma6"),
  Chr  = c("Chr06","Chr02","Chr01","Chr01","Chr06"),
  Pos  = c(40304883, 67882606, 68034103, 6748036, 697459)
)

chr_levels <- levels(qtl_data$Chr_string)

jpeg("outputs/QTL_colocalization.jpg", width=1600, height=1200, res=150)

segplot(
  reorder(factor(QTL_Id), Chr_string) ~ Start_QTL_v3 + End_QTL_v3 | Chr_string,
  data = qtl_data,
  groups = Trait,
  level = qtl_data$Trait,
  layout = c(10,1),
  draw.bands = FALSE,
  lty = 1, lwd = 6,
  col.regions = mycol,
  ylab = "",
  xlab = "",
  colorkey = NULL,
  
  scales = list(x=list(draw=FALSE)),
  panel = add_ma_genes(ma_genes, chr_levels)
)

dev.off()
