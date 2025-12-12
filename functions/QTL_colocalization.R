# === Fichier functions_plot_genes.R ===
library(lattice)
library(latticeExtra)

# Fonction pour ajouter des lignes verticales et noms des gènes
add_ma_genes <- function(ma_genes, chr_levels, col_line="gray", cex_text=1, font_text=2) {
  panel_func <- function(x, y, groups, ...) {
    panel.segplot(x, y, groups=groups, ...)
    
    chr_here <- panel.number()
    chr_name <- chr_levels[chr_here]
    
    g <- subset(ma_genes, Chr == chr_name)
    if(nrow(g) > 0){
      y_min <- min(y)
      y_offset <- 0.05 * diff(range(y))  # 5% sous le min
      y_text <- y_min - y_offset
      
      for(i in 1:nrow(g)){
        # ligne verticale
        panel.abline(v=g$Pos[i], col=col_line, lty=2, lwd=2)
        # texte sous le plot
        panel.text(
          x = g$Pos[i],
          y = y_text,
          labels = g$Gene[i],
          cex = cex_text,
          col = col_line,
          font = font_text
        )
      }
    }
  }
  return(panel_func)
}
