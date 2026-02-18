library(dplyr)
qtl <- read.csv("data/sheet5.csv", header = TRUE)

chr_length <- data.frame(
  Chromosome = 1:10,
  Length = c(80884392, 77742459, 74386277, 68658214, 71854669,
             61277060, 65505356, 62686529, 59416394, 61233695)
)

chr_length$Length_plot <- max(chr_length$Length)

buffer <- 12000000 / 2

# Function: genome position (GRAPHICAL)----

get_genome_pos_plot <- function(chr, pos, chr_length, buffer = 0) {
  offsets <- c(0, cumsum(chr_length$Length_plot[-nrow(chr_length)] + buffer))
  offsets[chr] + (pos / chr_length$Length[chr]) * chr_length$Length_plot[chr]
}


# QTL positions----
qtl$genome_start <- get_genome_pos_plot(qtl$Chromosome, qtl$QTL_start, chr_length, buffer)
qtl$genome_end <- get_genome_pos_plot(qtl$Chromosome, qtl$QTL_end, chr_length, buffer)

chr_length$genome_start <- get_genome_pos_plot(chr_length$Chromosome, 0, chr_length, buffer)
chr_length$genome_end <- get_genome_pos_plot(chr_length$Chromosome, chr_length$Length, chr_length, buffer)

x_max <- max(chr_length$genome_end)


trait_col <- c(FLA = "red", GMR = "darkgoldenrod", YLD = "darkcyan")

genes <- data.frame(
  Gene = c("Ma1", "Ma3", "Ma6", "Ma2", "Ma5", "Hd1", "CO", 
           "Ehd1", "GI", "CN12", "LHY", "CRY2"),
  Chr = c(6, 1, 6, 2, 1, 4, 10, 1, 3, 3, 7, 10),
  Start = c(40304883, 68034103, 697459, 67882606, 
            6748036, 56625893, 12353901, 21860030, 
            3821973, 62747946, 4745530, 22156857)
)

genes$genome_pos <- get_genome_pos_plot(genes$Chr, genes$Start, chr_length, buffer)

y_offsets <- rep(1, nrow(genes))
y_offsets[genes$Chr == 1] <- c(1,1,0.3)
y_offsets[genes$Chr == 2] <- 1
y_offsets[genes$Chr == 3] <- c(0.3, 1)
y_offsets[genes$Chr == 4] <- 1
y_offsets[genes$Chr == 10] <- c(0.3, 1)

pdf("outputs/QTL_colocalization.pdf", width = 6, height = 3)
#png("outputs/QTL_colocalization.png", width = 700, height = 400, res = 300)

par(
  mai = c(0.8,0.8, 0.3, 0.3),  # Marges réduites pour rapprocher les titres
  mgp = c(0.8, 0.3, 0),            # Ajuste la position des titres des axes
  cex.main = 1.2,              # Taille du titre principal
  cex.lab = 1.0,                # Taille des titres des axes
  cex.axis = 0.8                # Taille des étiquettes des axes
)
plot(
  qtl$QTL_number ~ qtl$genome_start,
  type = "n",
  xaxt = "n",
  yaxt = "n",
  xlab = "Chromosome",
  ylab = "Trait",
  ylim = c(0, max(qtl$QTL_number) +1),
  xlim = c(0, x_max)) #,
  #main = "QTL colocalization with known loci")

############################
# QTL rectangles
############################
rect(
  qtl$genome_start,
  qtl$QTL_number,
  qtl$genome_end,
  qtl$QTL_number + 0.35,
  col = trait_col[qtl$Trait],
  border = trait_col[qtl$Trait]
)

############################
# Chromosome bars
############################
for (i in 1:10) {
  rect(
    chr_length$genome_start[i],
    -0.5,
    chr_length$genome_end[i],
    0.15,
    col = c("grey30", "grey10")[i %% 2 + 1],
    border = NA
  )
}

############################
# Colocalizing genes (avec décalages)
############################
# Dessiner les segments (pointillés) jusqu'à la hauteur des étiquettes
for (i in 1:nrow(genes)) {
  segments(
    x0 = genes$genome_pos[i],
    y0 = 0,
    x1 = genes$genome_pos[i],
    y1 = max(qtl$QTL_number) + y_offsets[i],
    col = "gray",
    lty = 3,
    lwd = 1.2
  )
}

for (i in 1:nrow(genes)) {
  if (genes$Chr[i] == 10) {
    text(
      genes$genome_pos[i],
      max(qtl$QTL_number) + y_offsets[i],
      labels = genes$Gene[i],
      cex = 0.7,
      adj = 0,  
      pos = 2
    )
  } else {
    text(
      genes$genome_pos[i],
      max(qtl$QTL_number) + y_offsets[i],
      labels = genes$Gene[i],
      cex = 0.7,
      adj = 0
    )
  }
}

############################
# Legend
############################
legend(
  "topright",
  legend = names(trait_col),
  fill = trait_col,
  bty = "n",
  cex = 0.5
)
dev.off()

