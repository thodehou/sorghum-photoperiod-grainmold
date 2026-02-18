#Load library----
library(qtl)
# library(googlesheets4)
# library(googledrive)
# library(readxl)
# 
# #Load data from Google sheet
# #sheet_url<-"https://docs.google.com/spreadsheets/d/1scmIOMcw54usmJNhgj-z68Z1_rDVytwQtccfqrLHD5M/edit?usp=sharing"

GrGn_qtl <-read.cross(format = "csv", dir = "/Users/thierry/Documents/ThierryLab/these/data/",
                   file = "Rqtl_mapping.csv", genotypes = c("A", "B", "H"),
                   alleles = c("A", "B"), crosstype = "riself")

GrGn<- jittermap(GrGn_qtl, amount = 1e-6)

names(GrGn_qtl$pheno)

phenogeno<- calc.genoprob(GrGn_qtl, 
  step=0, off.end=0.0, error.prob=1.0e-4,stepwidth = "fixed", map.function="kosambi")

phenogeno_cross <- sim.geno(phenogeno,
  n.draws=32, step=0, off.end=0.0, error.prob=1.0e-4, stepwidth = "fixed", map.function="kosambi")


#FLA Bameby 22 ----
scan.cim_FLA_BY22 <- cim(phenogeno_cross, pheno.col=1, map.function="kosambi")
perm_FLA_B22 <- cim(phenogeno_cross,pheno.col= 1, map.function="kosambi", n.perm=1000)
summary(scan.cim_FLA_BY22)
summary(scan.cim.perm_FLA_B22)
plot(scan.cim_FLA_BY22, ylab = "LOD",  main = "QTL profile in BY22 for FLA")
summary(scan.cim_FLA_BY22, perms =  scan.cim.perm_FLA_B22, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =1 , infer = F, marker = "Sbv3.1_06_40312464", ylab = "FLA (BY22)")

qtl <- makeqtl(phenogeno_cross, chr=6, pos= 119.3 ,what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=1, 
                 formula = y~Q1 +GrGn_qtl$pheno$FLA_BY22 , qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

# Exemple corrigé
fitqtl2 <- fitqtl(
  cross = GrGn_qtl, 
  pheno.col = 14,  # ta variable de rendement ou autre
  qtl = qtl, 
  formula = y ~ Q1 + FLA_BY22,   # simple nom ici
  covar = data.frame(FLA_BY22 = GrGn_qtl$pheno$FLA_BY22), 
  method = "hk", 
  get.ests = TRUE
)

summary(fitqtl2)





lodint(results =  scan.cim_FLA_BY22, chr =6 , expandtomarkers = T)

# > summary(scan.cim_FLA_BY22)
# chr   pos    lod
# S1_68317108          1  82.7  5.297
# S2_3042598           2  83.7  1.586
# S3_66095714          3 329.5  2.233
# S4_68260182          4 269.0  1.887
# S5_62754800          5  23.6  6.162
# Sbv3.1_06_40312464   6 119.3 56.923
# S7_54789016          7  75.8  0.880
# S8_55843231          8  96.1  1.710
# S9_49803067          9 192.3  1.142
# snpSB00146          10 211.6  0.964
# > summary(scan.cim.perm_FLA_B22)
# LOD thresholds (1000 permutations)
# [,1]
# 5%  3.60
# 10% 3.21
# > summary(scan.cim_FLA_BY22, perms =  scan.cim.perm_FLA_B22, alpha=0.1, pvalues=TRUE)
# chr   pos   lod  pval
# S1_68317108          1  82.7  5.30 0.002
# S5_62754800          5  23.6  6.16 0.000
# Sbv3.1_06_40312464   6 119.3 56.92 0.000
# > 

# > summary(fitqtl)
# 
# fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 286 
# 
# Full model result
#   Model formula: y ~ Q1 
# 
# df         SS       MS      LOD     %var Pvalue(Chi2)   Pvalue(F)
# Model   1   843.3031 843.3031 1.622055 2.578019  0.006274108 0.006505783
# Error 284 31867.9758 112.2112                                           
# Total 285 32711.2789                                                    
# 
# 
# Estimated effects:
#   est      SE       t
# Intercept 76.3547  0.6312 120.972
# 1@82.7     1.7846  0.6510   2.741
# 
# > lodint(results =  scan.cim_FLA_BY22, chr =1 , expandtomarkers = T)
# chr      pos       lod
# S1_67862324   1 81.36691 3.5877476
# S1_68317108   1 82.73083 5.2967020
# S1_64013631   1 91.18710 0.1469714

# > summary(fitqtl)
# 
# fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 286 
# 
# Full model result
#   Model formula: y ~ Q1 
# 
# df        SS        MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  4268.491 4268.4913 8.683727 13.04899 2.552609e-10 3.051526e-10
# Error 284 28442.788  100.1507                                            
# Total 285 32711.279                                                      
# 
# 
# Estimated effects:
#   -----------------
#   est      SE       t
# Intercept 77.0049  0.5955 129.304
# 5@23.6    -3.9648  0.6073  -6.528
# 
# > lodint(results =  scan.cim_FLA_BY22, chr =5 , expandtomarkers = T)
# chr      pos         lod
# S5_4754454    5 11.64515 0.124026762
# S5_62754800   5 23.58505 6.161820272
# S5_64829267   5 29.58272 0.001886974

# > qtl <- makeqtl(phenogeno_cross, chr=6, pos= 119.3 ,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=1, 
#                    +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# > summary(fitqtl)
# 
# fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 286 
# 
# Full model result
# ----------------------------------  
#   Model formula: y ~ Q1 
# 
# df       SS          MS      LOD     %var Pvalue(Chi2) Pvalue(F)
# Model   1 18367.90 18367.90463 51.20061 56.15159            0         0
# Error 284 14343.37    50.50484                                         
# Total 285 32711.28                                                     
# 
# 
# Estimated effects:
#   -----------------
#   est      SE      t
# Intercept 75.2242  0.4261 176.54
# 6@119.3   -8.1993  0.4299 -19.07
# 
# > lodint(results =  scan.cim_FLA_BY22, chr =6 , expandtomarkers = T)
# chr      pos      lod
# S6_42183773          6 118.4694 52.58639
# Sbv3.1_06_40312464   6 119.3410 56.92301
# S6_39929198          6 119.7650 54.20940
# >


# Interactions QTL-----
mar <- find.marker(phenogeno, chr=6, pos=119.3)
mar2 <- find.marker(phenogeno, chr=1, pos=82.7)

plotPXG(phenogeno, marker=mar)
effectplot(phenogeno, mname1=mar)


effectplot(phenogeno, mname2="6@119.3", mname1="6@132")
effectplot(phenogeno, mname2="6@119.3", mname1="6@118")
effectplot(phenogeno, mname2="6@132", mname1="6@118")

##Multi) analyses----
phenogeno<- calc.genoprob(GrGn_qtl, 
  step=0, off.end=0.0, error.prob=1.0e-4,stepwidth = "fixed", map.function="kosambi")# to re-run again 

qtl <- makeqtl(phenogeno, chr=c(6,6,6), pos=c(118, 119, 132), what="prob")#create QTL object containing the loci on chr 1 and 6.

qtl

out.fq <- fitqtl(phenogeno, qtl=qtl, method="hk")#fit the two locus additive model as follows
summary(out.fq)

summary(fitqtl(phenogeno, qtl=qtl, method="hk", get.ests=TRUE, dropone=FALSE))#the estimated effects of the QTL

out.fqi1 <- fitqtl(phenogeno, qtl=qtl, method="hk", formula=y~Q1*Q2)#Check if interaction
out.fqi2 <- fitqtl(phenogeno, qtl=qtl, method="hk", formula=y~Q1+Q2+Q1:Q2)#another way
summary(out.fqi1)
summary(out.fqi2)

addint(phenogeno, qtl=qtl, method="hk")#addint function also to check interactions (usefull for more than 2 QTL)

rqtl <- refineqtl(phenogeno, qtl=qtl, method="hk")#refine the location of QTL (Multi-QTL analyses)
summary(rqtl)

summary(out.fqr <- fitqtl(phenogeno, qtl=rqtl, method="hk"))#re-run fitqtl to get the revised drop-one-term table.
#out.hk <- scanone(phenogeno, method="hk")

plotLodProfile(rqtl)#plot LOD 

plot(phenogeno, chr=c(6,6,6), col="red", add=TRUE)

qtl <- makeqtl(phenogeno_cross, chr=6, pos= 119.3 ,what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=1, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results =  scan.cim_by22, chr =6 , expandtomarkers = T)


out.em <- scanone(phenogeno)
summary(out.em)
phenogeno<- calc.genoprob(GrGn_qtl, 
  step=0, off.end=0.0, error.prob=1.0e-4,stepwidth = "fixed", map.function="kosambi")

qtl <- makeqtl(phenogeno, chr=c(6, 6), pos=c(119.3, 132), what="prob")
qtl

out.fq <- fitqtl(phenogeno, qtl=qtl, method="hk")
summary(out.fq)

rqtl <- refineqtl(phenogeno, qtl=qtl, method="hk", keeplodprofile = TRUE)


# Plot the refined QTL profile as the base plot
plotLodProfile(rqtl, col = "darksalmon", labelsep = TRUE,
               lwd = 2, lty = 1, main = "LOD Profiles for Two QTLs", 
               ylab = "LOD Score", xlab = "Position (cM)", xlim = c(0, 150))

abline(h = lod_threshold, col = "darkgreen", lty = 2, lwd = 2)




# Define your threshold level (e.g., 3 for a significant LOD score)
lod_threshold <- 3

# Define a base plot with the first QTL profile
plotLodProfile(rqtl, col = "darksalmon", labelsep = TRUE, lwd = 2, lty = 1, 
               main = "LOD Profiles for QTLs", ylab = "LOD Score", xlab = "Position (cM)", 
               xlim = c(0, 150))

# Add a threshold line across the plot
abline(h = lod_threshold, col = "darkgreen", lty = 2, lwd = 2)

# Add any additional QTL profiles if available, using `lines()` for distinct colors
# (Assuming rqtl2 and rqtl3 are other QTL profiles)
# lines(rqtl2$pos, rqtl2$lod, col = "steelblue", lwd = 2, lty = 2)
# lines(rqtl3$pos, rqtl3$lod, col = "purple", lwd = 2, lty = 3)

# Customize the legend
legend("topright", legend = c("QTL 1 (Refined)", "QTL 2", "QTL 3", paste("LOD Threshold =", lod_threshold)),
       col = c("darksalmon", "steelblue", "purple", "darkgreen"), lty = c(1, 2, 3, 2), lwd = c(2, 2, 2, 2))


## -------------------------------------------------------------------------------
names(GrGn_qtl$pheno)
scan.cim_FLA_b23 <- cim(phenogeno_cross, pheno.col=4, map.function="kosambi")
scan.cim.perm_FLA_b23 <- cim(phenogeno_cross,pheno.col= 4, map.function="kosambi", n.perm=1000)

summary(scan.cim_FLA_b23)
summary(scan.cim.perm_FLA_b23)
plot(scan.cim_FLA_b23, ylab = "LOD", main = "QTL profile in BY23 for FLA")

summary(scan.cim_FLA_b23, perms =  scan.cim.perm_FLA_b23, alpha=0.1, pvalues=TRUE)


plotPXG(phenogeno_cross, pheno.col =4 , marker = "Sbv3.1_06_40312464", ylab = "FLA (BY23)")

qtl <- makeqtl(phenogeno_cross, chr=6, pos= 119.3 ,what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=4, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results =  scan.cim_by23, chr =6 , expandtomarkers = T)


## -------------------------------------------------------------------------------
plotPXG(phenogeno_cross, pheno.col =4 , marker = "S10_13924090", ylab = "FLA (BY23)")

qtl <- makeqtl(phenogeno_cross, chr=10, pos= 207 ,what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=4, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results =  scan.cim_by23, chr =10 , expandtomarkers = T)

# FLA SM22-----
names(GrGn_qtl$pheno)
scan.cim_FLA_SM22 <- cim(phenogeno_cross, pheno.col=2, map.function="kosambi")
scan.cim.perm_FLA_SM22 <- cim(phenogeno_cross,pheno.col= 2, map.function="kosambi", n.perm=1000)

summary(scan.cim.perm_FLA_SM22)
summary(scan.cim_FLA_SM22)
summary(scan.cim_FLA_SM22, perms =  scan.cim.perm_FLA_SM22, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =2 , marker = "S7_12049713", ylab = "FLA (SM22)")

qtl <- makeqtl(phenogeno_cross, chr=7, pos=56,what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=2, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results =  scan.cim_FLA_SM22, chr =7 , expandtomarkers = T)



# > summary(scan.cim_FLA_SM22, perms =  scan.cim.perm_FLA_SM22, alpha=0.1, pvalues=TRUE)
#             chr   pos   lod  pval
# snpSB00330    2  75.3  3.58 0.072
# S3_52787725   3  33.0  6.57 0.000
# S5_59574508   5  19.2  5.57 0.002
# S6_1839684    6 131.8 16.30 0.000
# S7_12049713   7  56.0  3.64 0.065
# > qtl <- makeqtl(phenogeno_cross, chr=3, pos=33,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=2, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 1 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 285 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df       SS         MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1 1763.920 1763.91980 12.12985 17.79879 7.782663e-14 9.969803e-14
# Error 283 8146.416   28.78592                                            
# Total 284 9910.336                                                       
# 
# 
# Estimated effects:
# -----------------
#               est      SE       t
# Intercept 67.8839  0.3296 205.961
# 3@33.0     2.6595  0.3397   7.828
# 
# > qtl <- makeqtl(phenogeno_cross, chr=5, pos=19.2,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=2, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 1 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 285 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df       SS         MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1 1213.722 1213.72208 8.085192 12.24703  1.04788e-09 1.238468e-09
# Error 283 8696.614   30.73009                                            
# Total 284 9910.336                                                       
# 
# 
# Estimated effects:
# -----------------
#               est      SE       t
# Intercept 68.7845  0.3302 208.328
# 5@19.2    -2.1496  0.3420  -6.285
# 
# > lodint(results =  scan.cim_FLA_SM22, chr =3 , expandtomarkers = T)
#             chr      pos        lod
# S3_8607931    3 20.11218 0.09621547
# S3_52787725   3 33.01659 6.56777875
# S3_4555653    3 38.08470 2.24591698
# > lodint(results =  scan.cim_FLA_SM22, chr =5 , expandtomarkers = T)
#             chr      pos       lod
# S5_4754454    5 11.64515 0.4226116
# S5_59574508   5 19.18057 5.5713445
# S5_64396901   5 26.94659 2.1721124
# > lodint(results =  scan.cim_FLA_SM22, chr =3 , expandtomarkers = T)
#             chr      pos        lod
# S3_8607931    3 20.11218 0.09621547
# S3_52787725   3 33.01659 6.56777875
# S3_4555653    3 38.08470 2.24591698
# > lodint(results =  scan.cim_FLA_SM22, chr =é , expandtomarkers = T)
# Error: object 'é' not found
# 
# > lodint(results =  scan.cim_FLA_SM22, chr =2 , expandtomarkers = T)
#             chr      pos       lod
# S2_61537157   2 71.42158 0.4125169
# snpSB00330    2 75.31200 3.5759893
# S2_14172689   2 79.47247 2.0147507
# > qtl <- makeqtl(phenogeno_cross, chr=7, pos=56,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=2, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 1 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 285 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS        MS      LOD     %var Pvalue(Chi2)  Pvalue(F)
# Model   1  190.1692 190.16919 1.199089 1.918897   0.01877867 0.01930587
# Error 283 9720.1670  34.34688                                          
# Total 284 9910.3361                                                    
# 
# 
# Estimated effects:
# -----------------
#               est      SE       t
# Intercept 68.6939  0.3513 195.555
# 7@56.0     0.8338  0.3544   2.353
# 
# > lodint(results =  scan.cim_FLA_SM22, chr =7 , expandtomarkers = T)
#             chr      pos       lod
# S7_10020265   7 48.87481 0.6017292
# S7_12049713   7 56.00662 3.6363422
# S7_540003     7 86.04898 0.1062768
# >

#Multi) analyses----
phenogeno<- calc.genoprob(GrGn_qtl, 
  step=0, off.end=0.0, error.prob=1.0e-4,stepwidth = "fixed", map.function="kosambi")

qtl <- makeqtl(phenogeno, chr=6, pos=119, what="prob")
qtl

qtl1 <- makeqtl(phenogeno, chr=6, pos=131.81, what="prob")
qtl1

rqtl <- refineqtl(phenogeno, qtl=qtl, method="hk")

rqtl1 <- refineqtl(phenogeno, qtl=qtl1, method="hk")

plotLodProfile(rqtl1)
plot(out.em, out.hk, out.imp, col=c("blue", "red", "green"), chr=c(7,15))

plotLodProfile(rqtl, rqtl1, col=c("blue", "red"))
# Plot LOD profile for the first QTL
plotLodProfile(rqtl, col = "blue", lty = 1, main = "LOD Profiles for Two QTLs", ylab = "LOD Score", xlab = "Position (cM)")

# Add LOD profile for the second QTL
plotLodProfile(rqtl1, col = "red", lty = 2)

# Add a legend
legend("topright", legend = c("QTL 1", "QTL 2"), col = c("blue", "red"), lty = c(1, 2))


plotLodProfile(rqtl)
plot(out.hk, chr=c(7,15), col="red", add=TRUE)

#plot(out.em, out.hk, out.imp, col=c("blue", "red", "green"), chr=c(7,15))


#FLA SD22------
names(GrGn_qtl$pheno)
scan.cim_FLA_SD22 <- cim(phenogeno_cross, pheno.col=3, map.function="kosambi")
scan.cim.perm_FLA_SD22 <- cim(phenogeno_cross,pheno.col= 3, map.function="kosambi", n.perm=1000)

summary(scan.cim.perm_FLA_SD22)
summary(scan.cim_FLA_SD22)
plot(scan.cim_FLA_SD22, ylab = "LOD", main = "QTL profile in SD22 for FLA")
summary(scan.cim_FLA_SD22, perms =  scan.cim.perm_FLA_SD22, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =6 ,infer = F, marker = "S1_71516550", ylab = "FLA (SD22)")

qtl <- makeqtl(phenogeno_cross, chr=1, pos= 70,what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=3, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results =  scan.cim_FLA_SD22, chr =1 , expandtomarkers = T)


# > scan.cim_FLA_SD23 <- cim(phenogeno_cross, pheno.col=6, map.function="kosambi")
# > scan.cim_FLA_SD23 <- cim(phenogeno_cross, pheno.col=3, map.function="kosambi")
# > scan.cim.perm_FLA_SD23 <- cim(phenogeno_cross,pheno.col= 3, map.function="kosambi", n.perm=1000)
# > summary(scan.cim.perm_FLA_SD23)
# LOD thresholds (1000 permutations)
#     [,1]
# 5%  3.83
# 10% 3.42
# > summary(scan.cim_FLA_SD23)
#              chr    pos   lod
# S1_71516550    1  70.03  2.23
# S2_75160149    2  39.50  1.52
# S3_1239653     3 556.39  3.27
# S4_466340      4   0.00  0.93
# S5_11309581    5   2.93  1.24
# S6_1839684     6 131.81 11.60
# S7_54708717    7  71.72  1.06
# S8_58972170    8  55.36  1.83
# S9_2886889     9 141.64  2.03
# S10_58521233  10 642.53  1.06
# > summary(scan.cim_FLA_SD23, perms =  scan.cim.perm_FLA_SD23, alpha=0.1, pvalues=TRUE)
#            chr pos  lod pval
# S6_1839684   6 132 11.6    0
# > plotPXG(phenogeno_cross, pheno.col =6 , marker = "S6_1839684", ylab = "FLA (SD23)")
# > qtl <- makeqtl(phenogeno_cross, chr=6, pos= 132,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=3, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 3 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 283 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df       SS         MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1 1907.264 1907.26370 13.53843 19.77273  2.88658e-15 3.774758e-15
# Error 281 7738.668   27.53974                                            
# Total 282 9645.932                                                       
# 
# 
# Estimated effects:
# -----------------
#               est      SE       t
# Intercept 81.5900  0.3327 245.256
# 6@131.8   -2.8652  0.3443  -8.322
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 283 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df       SS         MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1 1907.264 1907.26370 13.53843 19.77273  2.88658e-15 3.774758e-15
# Error 281 7738.668   27.53974                                            
# Total 282 9645.932                                                       
# 
# 
# Estimated effects:
# -----------------
#               est      SE       t
# Intercept 81.5900  0.3327 245.256
# 6@131.8   -2.8652  0.3443  -8.322
# 
# > lodint(results =  scan.cim_SD23, chr =6 , expandtomarkers = T)
# Error: object 'scan.cim_SD23' not found
# 
# > lodint(results =  scan.cim_FLA_SD23, chr =6 , expandtomarkers = T)
#             chr      pos       lod
# S6_2454378    6 128.0386  4.821417
# S6_1839684    6 131.8088 11.595097
# S6_1839684    6 131.8088 11.595097
# 
# 
# > summary(fitqtl)
# 
# fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 283 
# 
# Full model result
# ----------------------------------  
#   Model formula: y ~ Q1 
# 
# df        SS        MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  610.6182 610.61822 4.018735 6.330319  1.69285e-05 1.844299e-05
# Error 281 9035.3138  32.15414                                            
# Total 282 9645.9320                                                      
# 
# 
# Estimated effects:
#   -----------------
#   est     SE       t
# Intercept 82.778  0.341 242.721
# 1@70.0     1.530  0.351   4.358
# 
# > lodint(results =  scan.cim_FLA_SD22, chr =1 , expandtomarkers = T)
# chr      pos        lod
# S1_74099911   1 59.74481 0.81277241
# S1_71516550   1 70.02649 3.64569337
# S1_74580234   1 77.85136 0.01499609

#FLA BY23-----
names(GrGn_qtl$pheno)
scan.cim_FLA_BY23 <- cim(phenogeno_cross, pheno.col=4, map.function="kosambi")
scan.cim.perm_FLA_BY23 <- cim(phenogeno_cross, pheno.col= 4, map.function="kosambi", n.perm=1000)

summary(scan.cim.perm_FLA_BY23)
summary(scan.cim_FLA_BY23)
plot(scan.cim_FLA_BY23, ylab = "LOD", main = "QTL profile in SM22 for FLA")
summary(scan.cim_FLA_BY23, perms =  scan.cim.perm_FLA_BY23, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =4 ,infer = F, marker = "S10_7557386", ylab = "FLA (BY23)")

qtl <- makeqtl(phenogeno_cross, chr=10, pos= 190,what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=4, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results =  scan.cim_FLA_BY23, chr = 10 , expandtomarkers = T)

# > summary(scan.cim_FLA_BY23)
#                    chr   pos    lod
# S1_66026375          1  84.6  2.641
# snpSB00330           2  75.3  2.339
# S3_66095714          3 329.5  2.525
# S4_67712021          4 268.4  0.891
# S5_62754800          5  23.6  2.118
# Sbv3.1_06_40312464   6 119.3 21.291
# S7_64331637          7 289.5  2.446
# S8_4142831           8  97.8  0.740
# S9_49803067          9 192.3  1.175
# S10_7557386         10 189.8  4.017
# > summary(scan.cim_FLA_BY23, perms =  scan.cim.perm_FLA_BY23, alpha=0.1, pvalues=TRUE)
#                    chr pos   lod  pval
# Sbv3.1_06_40312464   6 119 21.29 0.000
# S10_7557386         10 190  4.02 0.036
# > plotPXG(phenogeno_cross, pheno.col =4 , marker = "S10_7557386", ylab = "FLA (BY23)")
# > plotPXG(phenogeno_cross, pheno.col =4 , marker = "Sbv3.1_06_40312464", ylab = "FLA (BY23)")
# > qtl <- makeqtl(phenogeno_cross, chr=6, pos= 119,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=4, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 1 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 285 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df       SS          MS      LOD     %var Pvalue(Chi2) Pvalue(F)
# Model   1 10544.92 10544.91624 37.03106 45.02909            0         0
# Error 283 12873.09    45.48796                                         
# Total 284 23418.01                                                     
# 
# 
# Estimated effects:
# -----------------
#               est      SE      t
# Intercept 80.5705  0.4058 198.56
# 6@118.9   -6.2356  0.4095 -15.23
# 
# > lodint(results =  scan.cim_FLA_BY23, chr = 6 , expandtomarkers = T)
#                    chr      pos        lod
# S6_42793946          6 109.6862  0.1129956
# Sbv3.1_06_40312464   6 119.3410 21.2911709
# S6_33754933          6 120.2554 18.0909838
# > qtl <- makeqtl(phenogeno_cross, chr=10, pos= 190,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=4, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 1 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 285 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS        MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  2134.228 2134.2284 5.913921 9.113621 1.802037e-07 2.038628e-07
# Error 283 21283.779   75.2077                                            
# Total 284 23418.008                                                      
# 
# 
# Estimated effects:
# -----------------
#               est      SE       t
# Intercept 82.4359  0.5344 154.267
# 10@189.8  -2.9620  0.5560  -5.327
# 
# > lodint(results =  scan.cim_FLA_BY23, chr = 10 , expandtomarkers = T)
#             chr      pos        lod
# S10_5865655  10 177.9316 0.01193006
# S10_7557386  10 189.8498 4.01729308
# snpSB00146   10 211.5789 2.20642465
# >


#FLA SM23-----
names(GrGn_qtl$pheno)
scan.cim_FLA_SM23 <- cim(phenogeno_cross, pheno.col=5, map.function="kosambi")

scan.cim.perm_FLA_SM23 <- cim(phenogeno_cross,pheno.col= 5, map.function="kosambi", n.perm=1000)

summary(scan.cim.perm_FLA_SM23)
summary(scan.cim_FLA_SM23)
plot(scan.cim_FLA_sm23, ylab = "LOD", main = "QTL profile in SM23 for FLA")
summary(scan.cim_FLA_SM23, perms =  scan.cim.perm_FLA_SM23, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =5 , marker = "S3_1239653", ylab = "FLA (SM23)")

qtl <- makeqtl(phenogeno_cross, chr=6, pos= 118,what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=5, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results =  scan.cim_FLA_SM23, chr = 6 , expandtomarkers = T)

# LOD thresholds (1000 permutations)
#     [,1]
# 5%  4.00
# 10% 3.37
# > summary(scan.cim_FLA_SM23)
#             chr   pos   lod
# S1_80660369   1 126.5  2.71
# S2_75994648   2  34.7  2.03
# S3_1239653    3 556.4  4.75
# S4_2151872    4  13.6  1.07
# S5_62754800   5  23.6  1.43
# S6_42122367   6 118.5 14.45
# S7_58333166   7 258.8  1.42
# S8_55843231   8  96.1  1.13
# S9_55093999   9 216.8  1.03
# S10_7557386  10 189.8  2.48
# > summary(scan.cim_FLA_SM23, perms =  scan.cim.perm_FLA_SM23, alpha=0.1, pvalues=TRUE)
#             chr pos   lod  pval
# S3_1239653    3 556  4.75 0.013
# S6_42122367   6 118 14.45 0.000
# > plotPXG(phenogeno_cross, pheno.col =4 ,infer = F, marker = "Sbv3.1_06_40312464", ylab = "FLA (BY23)")
# > plotPXG(phenogeno_cross, pheno.col =4 ,infer = F, marker = "S10_7557386", ylab = "FLA (BY23)")
# > qtl <- makeqtl(phenogeno_cross, chr=3, pos= 556,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=5, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 4 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 282 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS         MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  1794.563 1794.56255 6.663917 10.31121 3.029597e-08 3.484694e-08
# Error 280 15609.438   55.74799                                            
# Total 281 17404.000                                                       
# 
# 
# Estimated effects:
# -----------------
#               est      SE       t
# Intercept 75.2853  0.4501 167.258
# 3@556.4    2.6167  0.4612   5.674
# 
# > lodint(results =  scan.cim_FLA_SM23, chr = 3 , expandtomarkers = T)
#             chr      pos      lod
# S3_69405187   3 328.2474 2.788079
# S3_1239653    3 556.3851 4.746922
# S3_1486771    3 557.2074 3.935114
# > qtl <- makeqtl(phenogeno_cross, chr=6, pos= 118,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=5, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 4 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 282 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS         MS      LOD     %var Pvalue(Chi2) Pvalue(F)
# Model   1  6584.078 6584.07808 29.10592 37.83083            0         0
# Error 280 10819.922   38.64258                                         
# Total 281 17404.000                                                    
# 
# 
# Estimated effects:
# -----------------
#               est      SE      t
# Intercept 74.0802  0.3753 197.38
# 6@118.5   -4.9432  0.3787 -13.05
# 
# > lodint(results =  scan.cim_FLA_SM23, chr = 6 , expandtomarkers = T)
#             chr      pos        lod
# S6_42793946   6 109.6862  0.1150201
# S6_42122367   6 118.4694 14.4502524
# S6_33754933   6 120.2554 11.8180782

#FLA SD23----
names(GrGn_qtl$pheno)
scan.cim_FLA_SD23 <- cim(phenogeno_cross, pheno.col=6, map.function="kosambi")
scan.cim.perm_FLA_SD23 <- cim(phenogeno_cross, pheno.col= 6, map.function="kosambi", n.perm=1000)

summary(scan.cim_FLA_SD23)
summary(scan.cim.perm_FLA_SD23)
plot(scan.cim_FLA_intercept, ylab = "LOD", main = "QTL profile for FLA (intercept)")
summary(scan.cim_FLA_SD23, perms =  scan.cim.perm_FLA_SD23, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =6 ,infer = F, marker = "S6_1839684", ylab = "FLA (SD23)")

qtl <- makeqtl(phenogeno_cross, chr=6, pos= 132,what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=6, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results =  scan.cim_FLA_SD23, chr = 6 , expandtomarkers = T)

# > scan.cim_FLA_SD23 <- cim(phenogeno_cross, pheno.col=6, map.function="kosambi")
# > scan.cim.perm_FLA_SD23 <- cim(phenogeno_cross, pheno.col= 6, map.function="kosambi", n.perm=1000)
# > summary(scan.cim_FLA_SD23)
#             chr   pos    lod
# S1_71516550   1  70.0  1.192
# S2_6271581    2  43.4  1.365
# S3_52794964   3  33.0  1.904
# S4_12258225   4  92.8  0.878
# S5_69997060   5  57.3  0.943
# S6_1839684    6 131.8 11.355
# S7_53013212   7  60.6  2.472
# S8_61278748   8  86.0  1.169
# S9_3907892    9 148.0  0.440
# S10_7557386  10 189.8  1.127
# > summary(scan.cim.perm_FLA_SD23)
# LOD thresholds (1000 permutations)
#     [,1]
# 5%  3.81
# 10% 3.34
# > summary(scan.cim_FLA_SD23, perms =  scan.cim.perm_FLA_SD23, alpha=0.1, pvalues=TRUE)
#            chr pos  lod pval
# S6_1839684   6 132 11.4    0
# > plotPXG(phenogeno_cross, pheno.col =6 , marker = "S6_1839684", ylab = "FLA (SD23)")
# > plotPXG(phenogeno_cross, pheno.col =6 ,infer = F, marker = "S6_1839684", ylab = "FLA (SD23)")
# > qtl <- makeqtl(phenogeno_cross, chr=6, pos= 132,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=6, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 5 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 281 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS        MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  915.6321 915.63212 12.92466 19.08863 1.210143e-14 1.576517e-14
# Error 279 3881.1083  13.91078                                            
# Total 280 4796.7404                                                      
# 
# 
# Estimated effects:
# -----------------
#               est      SE       t
# Intercept 76.4685  0.2375 321.950
# 6@131.8   -1.9948  0.2459  -8.113
# 
# > lodint(results =  scan.cim_FLA_SD23, chr = 6 , expandtomarkers = T)
#             chr      pos       lod
# S6_2454378    6 128.0386  5.081265
# S6_1839684    6 131.8088 11.354811
# S6_1839684    6 131.8088 11.354811
# >

#Intercept ----
names(GrGn_qtl$pheno)
scan.cim_FLA_intcp <- cim(phenogeno_cross, pheno.col=7, map.function="kosambi")

scan.cim.perm_FLA_intcp <- cim(phenogeno_cross, pheno.col= 7, map.function="kosambi", n.perm=1000)

summary(scan.cim.perm_FLA_intcp)
summary(scan.cim_FLA_intcp)
plot(scan.cim_FLA_intcp, ylab = "LOD", main = "QTL profile for FLA (Intercept)")
summary(scan.cim_FLA_intcp, perms =  scan.cim.perm_FLA_intcp, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =7 ,infer = F, marker = "S10_13924090", ylab = "FLA (Intercept)")

qtl <- makeqtl(phenogeno_cross, chr=10, pos= 207,what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=7, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results =  scan.cim_FLA_intcp, drop= 2, chr = 10 , expandtomarkers = T)
??lodint


# > scan.cim.perm_FLA_intcp <- cim(phenogeno_cross, pheno.col= 7, map.function="kosambi", n.perm=1000)
# > summary(scan.cim.perm_FLA_intcp)
# LOD thresholds (1000 permutations)
# [,1]
# 5%  3.77
# 10% 3.32
# > summary(scan.cim_FLA_intcp)
# chr   pos   lod
# S1_20120607    1 140.2  2.51
# snpSB00330     2  75.3  2.72
# S3_1239653     3 556.4  3.61
# S4_67712021    4 268.4  1.51
# S5_62754800    5  23.6  3.00
# S6_42122367    6 118.5 16.85
# S7_12049713    7  56.0  3.29
# S8_57450092    8  41.7  1.38
# S9_49803067    9 192.3  1.57
# S10_13924090  10 206.7  4.41
# > summary(scan.cim_FLA_intcp, perms =  scan.cim.perm_FLA_intcp, alpha=0.1, pvalues=TRUE)
# chr pos   lod  pval
# S3_1239653     3 556  3.61 0.068
# S6_42122367    6 118 16.85 0.000
# S10_13924090  10 207  4.41 0.013
# > plotPXG(phenogeno_cross, pheno.col =7 , marker = "S3_1239653", ylab = "FLA (Intercept)")
# > plotPXG(phenogeno_cross, pheno.col =7 ,infer = F, marker = "S3_1239653", ylab = "FLA (Intercept)")
# > plotPXG(phenogeno_cross, pheno.col =7 ,infer = F, marker = "S6_42122367", ylab = "FLA (Intercept)")
# > plotPXG(phenogeno_cross, pheno.col =7 ,infer = F, marker = "S10_13924090", ylab = "FLA (Intercept)")
# > qtl <- makeqtl(phenogeno_cross, chr=3, pos= 556,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=7, 
#                    +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# > summary(fitqtl)
# 
# fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 286 
# 
# Full model result
# ----------------------------------  
#   Model formula: y ~ Q1 
# 
# df        SS         MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  1003.722 1003.72173 6.324701 9.682623 6.780926e-08 7.731454e-08
# Error 284  9362.496   32.96653                                            
# Total 285 10366.218                                                       
# 
# 
# Estimated effects:
#   -----------------
#   est      SE       t
# Intercept 77.1520  0.3432 224.798
# 3@556.4    1.9398  0.3515   5.518
# 
# > lodint(results =  scan.cim_FLA_intcp, chr = 3 , expandtomarkers = T)
# chr       pos       lod
# S3_8607931   3  20.11218 0.5177858
# S3_1239653   3 556.38508 3.6127079
# S3_1486771   3 557.20735 3.0056742
# > ??lodint
# > lodint(results =  scan.cim_FLA_intcp, drop= 2, chr = 3 , expandtomarkers = T)
# chr       pos       lod
# S3_8607931   3  20.11218 0.5177858
# S3_1239653   3 556.38508 3.6127079
# S3_1486771   3 557.20735 3.0056742
# > qtl <- makeqtl(phenogeno_cross, chr=6, pos= 118,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=7, 
#                    +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# > summary(fitqtl)
# 
# fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 286 
# 
# Full model result
# ----------------------------------  
#   Model formula: y ~ Q1 
# 
# df        SS         MS      LOD     %var Pvalue(Chi2) Pvalue(F)
# Model   1  4642.923 4642.92335 36.89032 44.78898            0         0
# Error 284  5723.294   20.15244                                         
# Total 285 10366.218                                                    
# 
# 
# Estimated effects:
#   -----------------
#   est      SE      t
# Intercept 76.2113  0.2690 283.29
# 6@118.5   -4.1197  0.2714 -15.18
# 
# > lodint(results =  scan.cim_FLA_intcp, drop= 2, chr = 6 , expandtomarkers = T)
# chr      pos          lod
# S6_42793946   6 109.6862 6.995289e-04
# S6_42122367   6 118.4694 1.684985e+01
# S6_1839684    6 131.8088 1.620132e+01
# > qtl <- makeqtl(phenogeno_cross, chr=10, pos= 207,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=7, 
#                    +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# > summary(fitqtl)
# 
# fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 286 
# 
# Full model result
# ----------------------------------  
#   Model formula: y ~ Q1 
# 
# df          SS       MS      LOD      %var Pvalue(Chi2) Pvalue(F)
# Model   1    81.86931 81.86931 0.492427 0.7897704     0.132095 0.1337964
# Error 284 10284.34832 36.21249                                          
# Total 285 10366.21763                                                   
# 
# 
# Estimated effects:
#   -----------------
#   est      SE       t
# Intercept 76.7764  0.3618 212.205
# 10@206.7   0.5577  0.3709   1.504
# 
# > lodint(results =  scan.cim_FLA_intcp, drop= 2, chr = 10 , expandtomarkers = T)
# chr      pos       lod
# S10_5865655   10 177.9316 0.1105332
# S10_13924090  10 206.6713 4.4142957
# S10_47071511  10 215.1906 1.7809632
# >


#Slope ----
names(GrGn_qtl$pheno)
scan.cim_FLA_slope<- cim(phenogeno_cross, pheno.col=8, map.function="kosambi")

scan.cim.perm_FLA_slope <- cim(phenogeno_cross, pheno.col= 8, map.function="kosambi", n.perm=1000)

summary(scan.cim.perm_FLA_slope)
summary(scan.cim_FLA_slope)
plot(scan.cim_FLA_intcp, ylab = "LOD", main = "QTL profile for FLA (Intercept)")
summary(scan.cim_FLA_slope, perms =  scan.cim.perm_FLA_slope, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =8, infer =F,  marker = "Sbv3.1_06_40312464", ylab = "FLA (slope)")

qtl <- makeqtl(phenogeno_cross, chr=6, pos= 119.3,what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=8, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results =  scan.cim_FLA_slope, chr = 6 , expandtomarkers = T)



# > summary(scan.cim_FLA_slope, perms =  scan.cim.perm_FLA_slope, alpha=0.1, pvalues=TRUE)
# chr   pos   lod  pval
# S1_66111434          1  83.8  6.19 0.001
# Sbv3.1_06_40312464   6 119.3 45.99 0.000
# > plotPXG(phenogeno_cross, pheno.col =8 , marker = "S1_66111434", ylab = "FLA (slope)")
# > plotPXG(phenogeno_cross, pheno.col =8, infer =F,  marker = "S1_66111434", ylab = "FLA (slope)")
# > plotPXG(phenogeno_cross, pheno.col =8, infer =F,  marker = "Sbv3.1_06_40312464", ylab = "FLA (slope)")
# > qtl <- makeqtl(phenogeno_cross, chr=1, pos= 83.8,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=8, 
#                    +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# > summary(fitqtl)
# 
# fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 286 
# 
# Full model result
# ----------------------------------  
#   Model formula: y ~ Q1 
# 
# df         SS        MS      LOD     %var Pvalue(Chi2) Pvalue(F)
# Model   1   315.6618 315.66180 1.007471 1.609142   0.03124273 0.0319899
# Error 284 19301.1166  67.96168                                         
# Total 285 19616.7784                                                   
# 
# 
# Estimated effects:
#   -----------------
#   est      SE      t
# Intercept 11.6839  0.4956 23.574
# 1@83.8     1.0763  0.4994  2.155
# 
# > lodint(results =  scan.cim_slope, chr = 1 , expandtomarkers = T)
# Error: object 'scan.cim_slope' not found
# 
# > lodint(results =  scan.cim_FLA_slope, chr = 1 , expandtomarkers = T)
# chr      pos      lod
# S1_67862324   1 81.36691 3.450823
# S1_66111434   1 83.84308 6.194297
# S1_65029310   1 84.91578 4.491301
# > qtl <- makeqtl(phenogeno_cross, chr=6, pos= 119.3,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=8, 
#                    +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# > summary(fitqtl)
# 
# fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 286 
# 
# Full model result
# ----------------------------------  
#   Model formula: y ~ Q1 
# 
# df        SS         MS      LOD     %var Pvalue(Chi2) Pvalue(F)
# Model   1  9014.353 9014.35297 38.21281 45.95226            0         0
# Error 284 10602.425   37.33248                                         
# Total 285 19616.778                                                    
# 
# 
# Estimated effects:
#   -----------------
#   est      SE      t
# Intercept 10.9356  0.3663  29.85
# 6@119.3   -5.7440  0.3697 -15.54
# 
# > lodint(results =  scan.cim_FLA_slope, chr = 6 , expandtomarkers = T)
# chr      pos      lod
# S6_40634627          6 118.9103 43.17125
# Sbv3.1_06_40312464   6 119.3410 45.98986
# S6_39929198          6 119.7650 44.30467
# >

#PGMR-----
#BY22----
names(GrGn_qtl$pheno)
scan.cim_PGMR_BY22<- cim(phenogeno_cross, pheno.col=9, map.function="kosambi")
perm_PGMR_BY22 <- cim(phenogeno_cross, pheno.col= 9, map.function="kosambi", n.perm=1000)
summary(scan.cim_PGMR_BY22)
# > summary(scan.cim_PGMR_BY22)
# chr   pos    lod
# S1_778962      1 115.8  1.732
# S2_3042598     2  83.7  2.074
# S3_4555653     3  38.1  1.805
# S4_59554213    4  67.2  5.612
# S5_63531462    5 109.9  0.444
# S6_40634627    6 118.9 25.420
# S7_55785125    7  70.8  1.345
# S8_58972170    8  55.4  1.994
# S9_53892157    9 111.6  1.212
# S10_56132414  10 249.5  2.933
summary(perm_PGMR_BY22)
# > summary(perm_PGMR_BY22)
# LOD thresholds (1000 permutations)
# [,1]
# 5%  3.93
# 10% 3.41
plot(scan.cim_PGMR_BY22, ylab = "LOD", main = "QTL profile PGMR (BY)")
summary(scan.cim_PGMR_BY22, perms =  perm_PGMR_BY22, alpha=0.1, pvalues=TRUE)
# > summary(scan.cim_PGMR_BY22, perms =  perm_PGMR_BY22, alpha=0.1, pvalues=TRUE)
# chr   pos   lod  pval
# S4_59554213   4  67.2  5.61 0.002
# S6_40634627   6 118.9 25.42 0.000
plotPXG(phenogeno_cross, pheno.col =9 ,infer = F,  marker = "S4_59554213", ylab = "PGMR (BY22)")
??plotPXG

qtl <- makeqtl(phenogeno_cross, chr=4, pos= 67.2,what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=9, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

effectplot(phenogeno_cross, mname1="4@67.2", mname2="6@118.9")

# fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 286 
# 
# Full model result
# ----------------------------------  
#   Model formula: y ~ Q1 
# 
# df       SS         MS      LOD     %var Pvalue(Chi2) Pvalue(F)
# Model   1 14.11948 14.1194771 17.66611 24.75793            0         0
# Error 284 42.91065  0.1510938                                         
# Total 285 57.03012                                                    
# 
# 
# Estimated effects:
#   -----------------
#   est      SE      t
# Intercept 2.19098 0.02336 93.802
# 6@118.9   0.22789 0.02357  9.667

lodint(results = scan.cim_PGMR_BY22, chr = 4 , expandtomarkers = T)

# > lodint(results = scan.cim_PGMR_BY22, chr = 6 , expandtomarkers = T)
# chr      pos        lod
# S6_42793946   6 109.6862  0.6242789
# S6_40634627   6 118.9103 25.4196591
# S6_39929198   6 119.7650 23.8345659
#______

# > lodint(results = scan.cim_PGMR_BY22, chr = 4 , expandtomarkers = T)
# chr      pos        lod
# S4_60890518   4 66.08522 0.08627806
# S4_59554213   4 67.18582 5.61165690
# S4_52124078   4 78.41874 0.58776115
# > 
# > summary(fitqtl)
# 
# fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 286 
# 
# Full model result
# ----------------------------------  
#   Model formula: y ~ Q1 
# 
# df         SS        MS      LOD      %var Pvalue(Chi2) Pvalue(F)
# Model   1  0.2612024 0.2612024 0.285095 0.4580078    0.2518682 0.2539508
# Error 284 56.7689216 0.1998906                                          
# Total 285 57.0301240                                                    
# 
# 
# Estimated effects:
#   -----------------
#   est       SE      t
# Intercept  2.15687  0.02696 79.993
# 4@67.2    -0.03116  0.02726 -1.143
# 
#SD22----
names(GrGn_qtl$pheno)
scan.cim_PGMR_SD22<- cim(phenogeno_cross, pheno.col=11, map.function="kosambi")

perm_PGMR_SD22 <- cim(phenogeno_cross, pheno.col= 11, map.function="kosambi", n.perm=1000)

summary(scan.cim_PGMR_SD22)
summary(perm_PGMR_SD22)
plot(scan.cim_PGMR_SD22, ylab = "LOD", main = "QTL profile PGMR (SD)")
summary(scan.cim_PGMR_SD22, perms =  perm_PGMR_SD22, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col = 11,infer = F,  marker = "S6_1839684", ylab = "PGMR (SD22)")

qtl <- makeqtl(phenogeno_cross, chr=6, pos= 132, what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=11, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results = scan.cim_PGMR_SD22, chr = 6 , expandtomarkers = T)

# > scan.cim_PGMR_SD22<- cim(phenogeno_cross, pheno.col=11, map.function="kosambi")
# > perm_PGMR_SD22 <- cim(phenogeno_cross, pheno.col= 11, map.function="kosambi", n.perm=1000)
# > summary(scan.cim_PGMR_SD22)
# chr   pos   lod
# S1_6666029    1 164.6 1.794
# S2_6271581    2  43.4 0.659
# S3_66095714   3 329.5 0.531
# S4_7147124    4  79.7 0.405
# S5_59574508   5  19.2 1.087
# S6_1839684    6 131.8 7.310
# S7_53916335   7  65.9 3.362
# S8_58332517   8  48.0 1.502
# S9_58176710   9   0.0 0.992
# S10_7557386  10 189.8 0.917
# > summary(perm_PGMR_SD22)
# LOD thresholds (1000 permutations)
# [,1]
# 5%  3.81
# 10% 3.42
# > plot(scan.cim_PGMR_SD22, ylab = "LOD", main = "QTL profile PGMR (SD)")
# > summary(scan.cim_PGMR_SD22, perms =  perm_PGMR_SD22, alpha=0.1, pvalues=TRUE)
# chr pos  lod pval
# S6_1839684   6 132 7.31    0
# > plotPXG(phenogeno_cross, pheno.col = 11,infer = F,  marker = "S6_1839684", ylab = "PGMR (SD22)")
# > qtl <- makeqtl(phenogeno_cross, chr=11, pos= 132,what=c("prob")) 
# Error in makeqtl(phenogeno_cross, chr = 11, pos = 132, what = c("prob")) : 
#   There's no chromosome number 11 in input cross object
# > qtl <- makeqtl(phenogeno_cross, chr=6, pos= 132, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=11, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 3 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 283 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS        MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  6.827152 6.8271521 9.858011 14.82109 1.608136e-11 1.972444e-11
# Error 281 39.236612 0.1396321                                            
# Total 282 46.063764                                                      
# 
# 
# Estimated effects:
# -----------------
#               est      SE       t
# Intercept 2.41691 0.02369 102.031
# 6@131.8   0.17142 0.02452   6.992
# 
# > # Total 285 57.0301240                                                    
# > # 
# > # 
# > # Estimated effects:
# > #   -----------------
# > #   est       SE      t
# > # Intercept  2.15687  0.02696 79.993
# > # 4@67.2    -0.03116  0.02726 -1.143
# > # 
# > #SD-------
# > names(GrGn_qtl$pheno)
#  [1] "FLA_BY22"  "FLA_SM22"  "FLA_SD22"  "FLA_BY23"  "FLA_SM23"  "FLA_SD23"  "Intcp_FLA"
#  [8] "slope_FLA" "PGMR_BY22" "PGMR_SM22" "PGMR_SD22" "PGMR_SM23" "PGMR_SD23" "PGMR_SM"  
# [15] "PGMR_SD"   "YLD_BY22"  "YLD_SM22"  "YLD_SD22"  "YLD_SM23"  "YLD_SD23"  "YLD_SM"   
# [22] "YLD_SD"    "TGMR_BY22" "TGMR_SD22" "TGMR_SM22" "TGMR_SD23" "TGMR_SM23" "TGMR_SM"  
# [29] "TGMR_SD"  
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=11, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 3 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 283 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS        MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  6.827152 6.8271521 9.858011 14.82109 1.608136e-11 1.972444e-11
# Error 281 39.236612 0.1396321                                            
# Total 282 46.063764                                                      
# 
# 
# Estimated effects:
# -----------------
#               est      SE       t
# Intercept 2.41691 0.02369 102.031
# 6@131.8   0.17142 0.02452   6.992
# 
# > summary(scan.cim_PGMR_SD22, perms =  perm_PGMR_SD22, alpha=0.1, pvalues=TRUE)
#            chr pos  lod pval
# S6_1839684   6 132 7.31    0
# > lodint(results = scan.cim_PGMR_SD22, chr = 6 , expandtomarkers = T)
#             chr      pos      lod
# S6_2454378    6 128.0386 4.550008
# S6_1839684    6 131.8088 7.309616
# S6_1839684    6 131.8088 7.309616

# > summary(scan.cim_PGMR_SD23)
# chr   pos   lod
# S1_61874941    1 100.6 0.726
# S2_3937345     2  29.5 0.683
# S3_3441687     3  39.2 1.470
# S4_61821462    4  64.7 0.608
# S5_66360856    5  39.1 1.116
# S6_5573833     6 126.4 0.452
# S7_1381634     7  86.0 1.363
# S8_4142831     8  97.8 0.574
# S9_6375482     9 158.3 0.575
# S10_52821474  10 228.1 0.978

#SM22----
names(GrGn_qtl$pheno)
scan.cim_PGMR_SM22<- cim(phenogeno_cross, pheno.col=10, map.function="kosambi")

perm_PGMR_SM22 <- cim(phenogeno_cross, pheno.col= 10, map.function="kosambi", n.perm=1000)

summary(scan.cim_PGMR_SM22)
summary(perm_PGMR_SM22)
plot(scan.cim_PGMR_SM23, ylab = "LOD", main = "QTL profile PGMR (SM)")
summary(scan.cim_PGMR_SM22, perms =  perm_PGMR_SM22, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col = 10,infer = F,  marker = "S6_33754933", ylab = "PGMR (SM22)")

qtl <- makeqtl(phenogeno_cross, chr=6, pos= 120,what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=10, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results = scan.cim_PGMR_SM22, chr = 6 , expandtomarkers = T)



# > summary(scan.cim_PGMR_SM22, perms =  perm_PGMR_SM22, alpha=0.1, pvalues=TRUE)
# chr pos  lod pval
# S6_33754933   6 120 8.15    0
# > summary(fitqtl)
# 
# fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 285 
# 
# Full model result
# ----------------------------------  
#   Model formula: y ~ Q1 
# 
# df       SS          MS     LOD     %var Pvalue(Chi2) Pvalue(F)
# Model   1 10.03046 10.03045708 20.0335 27.65408            0         0
# Error 283 26.24070  0.09272334                                        
# Total 284 36.27116                                                    
# 
# 
# Estimated effects:
#   -----------------
#   est      SE     t
# Intercept 2.37126 0.01833 129.4
# 6@119.8   0.19223 0.01848  10.4
# 
# > lodint(results = scan.cim_PGMR_SM22, chr = 6 , expandtomarkers = T)
# chr      pos        lod
# S6_42793946   6 109.6862 0.13585662
# S6_33754933   6 120.2554 8.14563992
# S6_6110410    6 126.0970 0.05824539
# >


#SM23----
names(GrGn_qtl$pheno)
scan.cim_PGMR_SM23<- cim(phenogeno_cross, pheno.col=12, map.function="kosambi")

perm_PGMR_SM23 <- cim(phenogeno_cross, pheno.col= 12, map.function="kosambi", n.perm=1000)

summary(scan.cim_PGMR_SM23)
summary(perm_PGMR_SM23)
plot(scan.cim_PGMR_SM23, ylab = "LOD", main = "QTL profile PGMR (SM)")
??plotPXG
summary(scan.cim_PGMR_SM23, perms =  perm_PGMR_SM23, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col = 12,infer = F,  marker = "S6_14495938", ylab = "PGMR (SM23)")

qtl <- makeqtl(phenogeno_cross, chr=9, pos= 111,what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=12, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results = scan.cim_PGMR_SM23, chr = 9 , expandtomarkers = T)

# 
# > summary(scan.cim_PGMR_SM23, perms =  perm_PGMR_SM23, alpha=0.1, pvalues=TRUE)
#             chr pos  lod  pval
# S6_14495938   6 120 7.07 0.000
# S7_59404098   7 265 3.74 0.059
# S9_1190784    9 111 3.75 0.057
# > qtl <- makeqtl(phenogeno_cross, chr=6, pos= 120,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=12, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 5 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 281 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df       SS        MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  4.65315 4.6531495 3.704383 5.890322 3.623025e-05 3.924164e-05
# Error 279 74.34337 0.2664637                                            
# Total 280 78.99652                                                      
# 
# 
# Estimated effects:
# -----------------
#                est       SE      t
# Intercept  2.36852  0.03131 75.645
# 6@119.8   -0.13191  0.03157 -4.179
# 
# > lodint(results = scan.cim_PGMR_SM23, chr = 6 , expandtomarkers = T)
#             chr      pos       lod
# S6_42793946   6 109.6862 0.1210560
# S6_14495938   6 120.2555 7.0709284
# S6_6110410    6 126.0970 0.3625298
# > qtl <- makeqtl(phenogeno_cross, chr=7, pos= 265,what=c("prob")) 
# > lodint(results = scan.cim_PGMR_SM23, chr = 7 , expandtomarkers = T)
#             chr      pos        lod
# S7_58325570   7 258.8445 0.09679025
# S7_59404098   7 264.9164 3.73910443
# S7_59597652   7 269.9756 0.40484008
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 281 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df       SS        MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  4.65315 4.6531495 3.704383 5.890322 3.623025e-05 3.924164e-05
# Error 279 74.34337 0.2664637                                            
# Total 280 78.99652                                                      
# 
# 
# Estimated effects:
# -----------------
#                est       SE      t
# Intercept  2.36852  0.03131 75.645
# 6@119.8   -0.13191  0.03157 -4.179
# 
# > qtl <- makeqtl(phenogeno_cross, chr=9, pos= 111,what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=12, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 5 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 281 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS        MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  4.421243 4.4212427 3.514338 5.596756 5.747373e-05 6.200838e-05
# Error 279 74.575276 0.2672949                                            
# Total 280 78.996518                                                      
# 
# 
# Estimated effects:
# -----------------
#                est       SE      t
# Intercept  2.25285  0.04610 48.870
# 9@111.4   -0.18888  0.04644 -4.067
# 
# > lodint(results = scan.cim_PGMR_SM23, chr = 9 , expandtomarkers = T)
#             chr      pos        lod
# S9_58176710   9   0.0000 0.08613330
# S9_1190784    9 111.3806 3.75060275
# S9_1646345    9 119.1402 0.03154426
# >

#Boucle plot effect PGMR------
-----------------------------
  # 1. Tableau : marqueur + colonne + label
  # -----------------------------
qtl_info <- data.frame(
  markers <- c("S4_59554213", "S6_40634627", "S6_1839684", "S6_33754933",
               "S6_14495938", "S7_59404098"),
  pheno_col = c(9, 9,
                11, 
                10,
                12, 12),
  ylab = c("PGMR (BY22)", "PGMR (BY22)",
           "PGMR (SD22)", 
           "PGMR (SM22)",
           "PGMR (SM23)",  "PGMR (SM23)")
)

# -----------------------------
# 2. Fichier haute résolution
# -----------------------------
png("Output/plot_PGMR.png", 
    width = 3000, height = 2500, res = 300)

# -----------------------------
# 3. Paramètres graphiques
# -----------------------------
par(mfrow = c(2, 3),
    mar = c(4.5, 4, 2, 1),
    oma = c(1, 1, 1, 1))
# -----------------------------
# 4. Boucle de tracé
# -----------------------------
# Boucle de dessin
for(i in 1:nrow(qtl_info)){
  plotPXG(
    phenogeno_cross,
    pheno.col = qtl_info$pheno_col[i],
    marker    = qtl_info$marker[i],
    infer     = FALSE,
    ylab      = qtl_info$ylab[i],
    cex.lab   = 1.3,  # taille des labels X/Y
    cex.axis  = 0.9  # titre (optionnel)
  )
}

dev.off()


#TGMR------
#BY22
names(GrGn_qtl$pheno)
scan.cim_TGMR_BY22<- cim(phenogeno_cross, pheno.col=19, map.function="kosambi")
perm_TGMR_BY22 <- cim(phenogeno_cross, pheno.col= 19, map.function="kosambi", n.perm=1000)
summary(scan.cim_TGMR_BY22)
summary(perm_TGMR_BY22)

plot(scan.cim_PGMR_BY22, ylab = "LOD", main = "QTL profile PGMR (BY)")
summary(scan.cim_TGMR_BY22, perms =  perm_TGMR_BY22, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =19 ,infer = F,  marker = "S4_59554213", ylab = "TGMR (BY22)")
??plotPXG

qtl <- makeqtl(phenogeno_cross, chr=6, pos=99.2, what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=23, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results = scan.cim_PGMR_BY22, chr = 4 , expandtomarkers = T)

# > summary(scan.cim_TGMR_BY22)
#              chr    pos   lod
# S1_54540069    1 126.88 1.796
# S2_57680557    2  97.75 3.099
# S3_12391410    3   9.85 0.766
# S4_4638445     4  16.79 0.700
# S5_3022709     5   3.85 0.433
# S6_45670817    6  99.23 3.708
# S7_63998483    7 288.19 1.555
# snpSB00309     8  78.76 0.813
# S9_5519208     9 156.77 0.813
# S10_52821474  10 228.15 0.480
# > summary(perm_TGMR_BY22)
# LOD thresholds (1000 permutations)
#     [,1]
# 5%  4.00
# 10% 3.51
# > summary(scan.cim_TGMR_BY22, perms =  perm_TGMR_BY22, alpha=0.1, pvalues=TRUE)
#             chr  pos  lod  pval
# S6_45670817   6 99.2 3.71 0.075
# 
# > plotPXG(phenogeno_cross, pheno.col =23 ,infer = F,  marker = "S6_45670817", ylab = "TGMR (BY22)")
# > qtl <- makeqtl(phenogeno_cross, chr=6, pos=99.2, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=23, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 2 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 284 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS        MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1   7.37571 7.3757098 4.078956 6.400192 1.463682e-05 1.596098e-05
# Error 282 107.86630 0.3825046                                            
# Total 283 115.24201                                                      
# 
# 
# Estimated effects:
# -----------------
#               est      SE      t
# Intercept 2.94250 0.03673 80.103
# 6@99.2    0.16629 0.03787  4.391
# 
# > lodint(results = scan.cim_PGMR_BY22, chr = 4 , expandtomarkers = T)
#             chr      pos        lod
# S4_60890518   4 66.08522 0.08627806
# S4_59554213   4 67.18582 5.61165690
# S4_52124078   4 78.41874 0.58776115


#SD22
names(GrGn_qtl$pheno)

scan.cim_TGMR_SD22<- cim(phenogeno_cross, pheno.col=20, map.function="kosambi")
perm_TGMR_SD22 <- cim(phenogeno_cross, pheno.col= 20, map.function="kosambi", n.perm=1000)
summary(scan.cim_TGMR_SD22)
summary(perm_TGMR_SD22)

plot(scan.cim_TGMR_SD22, ylab = "LOD", main = "QTL profile PGMR (BY)")
summary(scan.cim_TGMR_SD22, perms =  perm_TGMR_SD22, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =20 ,infer = F,  marker = "S7_54708458", ylab = "TGMR (SD22)")
??plotPXG

qtl <- makeqtl(phenogeno_cross, chr=7, pos=71, what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=24, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results = scan.cim_TGMR_SD22, chr = 7, expandtomarkers = T)

# > perm_TGMR_SD22 <- cim(phenogeno_cross, pheno.col= 24, map.function="kosambi", n.perm=1000)
# > summary(scan.cim_TGMR_SD22)
#             chr   pos   lod
# S1_6666029    1 164.6 2.178
# S2_4967260    2  47.6 2.253
# S3_55370875   3 106.8 0.673
# S4_466340     4   0.0 1.666
# S5_59574508   5  19.2 1.836
# S6_1839684    6 131.8 7.491
# S7_54708458   7  71.0 4.601
# S8_6118586    8 106.2 2.591
# S9_58176710   9   0.0 1.384
# S10_7557386  10 189.8 1.679
# > summary(perm_TGMR_SD22)
# LOD thresholds (1000 permutations)
#     [,1]
# 5%  3.74
# 10% 3.34
# > summary(scan.cim_TGMR_SD22, perms =  perm_TGMR_SD22, alpha=0.1, pvalues=TRUE)
#             chr pos  lod  pval
# S6_1839684    6 132 7.49 0.000
# S7_54708458   7  71 4.60 0.013
# > plotPXG(phenogeno_cross, pheno.col =24 ,infer = F,  marker = "S6_1839684", ylab = "TGMR (SD22)")
# > qtl <- makeqtl(phenogeno_cross, chr=6, pos=132, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=24, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 4 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 282 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS         MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  23.87872 23.8787193 7.687939 11.79857 2.678994e-09 3.146336e-09
# Error 280 178.50778  0.6375278                                            
# Total 281 202.38650                                                       
# 
# 
# Estimated effects:
# -----------------
#               est      SE     t
# Intercept 2.84200 0.05066 56.10
# 6@131.8   0.32090 0.05243  6.12
# 
# > lodint(results = scan.cim_TGMR_SD22, chr = 6 , expandtomarkers = T)
#             chr      pos      lod
# S6_2454378    6 128.0386 4.404602
# S6_1839684    6 131.8088 7.491365
# S6_1839684    6 131.8088 7.491365
# > plotPXG(phenogeno_cross, pheno.col =24 ,infer = F,  marker = "S7_54708458", ylab = "TGMR (SD22)")
# > qtl <- makeqtl(phenogeno_cross, chr=7, pos=71, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=24, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 4 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 282 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df       SS         MS      LOD     %var Pvalue(Chi2)   Pvalue(F)
# Model   1  17.3384 17.3383959 5.484463 8.566973 5.018012e-07 5.63472e-07
# Error 280 185.0481  0.6608861                                           
# Total 281 202.3865                                                      
# 
# 
# Estimated effects:
# -----------------
#                est       SE      t
# Intercept  2.73402  0.04841 56.475
# 7@71.0    -0.25024  0.04886 -5.122
# 
# > lodint(results = scan.cim_TGMR_SD22, chr = 7, expandtomarkers = T)
#             chr      pos         lod
# S7_53916335   7 65.92054 0.004738298
# S7_54708458   7 71.04548 4.600528891
# S7_54789016   7 75.81080 3.060081099
# > summary(scan.cim_TGMR_SD22, perms =  perm_TGMR_SD22, alpha=0.1, pvalues=TRUE)
#             chr pos  lod  pval
# S6_1839684    6 132 7.49 0.000
# S7_54708458   7  71 4.60 0.013


#SD23
names(GrGn_qtl$pheno)

scan.cim_TGMR_SD23<- cim(phenogeno_cross, pheno.col=22, map.function="kosambi")
perm_TGMR_SD23 <- cim(phenogeno_cross, pheno.col= 22, map.function="kosambi", n.perm=1000)
summary(scan.cim_TGMR_SD23)
summary(perm_TGMR_SD23)

plot(scan.cim_TGMR_SM22, ylab = "LOD", main = "QTL profile PGMR (BY)")
summary(scan.cim_TGMR_SD23, perms =  perm_TGMR_SD23, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =22 ,infer = F,  marker = "S10_7557386", ylab = "TGMR (SD23)")
??plotPXG

qtl <- makeqtl(phenogeno_cross, chr=10, pos=189.8, what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=22, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results = scan.cim_TGMR_SD23, chr = 10, expandtomarkers = T)

# > summary(scan.cim_TGMR_SM22)
# chr   pos   lod
# S1_61874941   1 100.6 7.297
# S2_57728127   2  97.5 2.197
# S3_232887     3 546.7 2.074
# S4_466340     4   0.0 3.816
# S5_66360856   5  39.1 1.005
# S6_42122367   6 118.5 9.752
# S7_54708717   7  71.7 1.335
# S8_59643388   8  60.0 1.138
# S9_49803067   9 192.3 0.331
# S10_5628553  10 175.0 0.625
# > summary(perm_TGMR_SM22)
# LOD thresholds (1000 permutations)
# [,1]
# 5%  3.71
# 10% 3.38
# > summary(scan.cim_TGMR_SM22, perms =  perm_TGMR_SM22, alpha=0.1, pvalues=TRUE)
# chr pos  lod  pval
# S1_61874941   1 101 7.30 0.000
# S4_466340     4   0 3.82 0.038
# S6_42122367   6 118 9.75 0.000

#SM22
names(GrGn_qtl$pheno)

scan.cim_TGMR_SM22<- cim(phenogeno_cross, pheno.col=21, map.function="kosambi")
perm_TGMR_SM22 <- cim(phenogeno_cross, pheno.col= 21, map.function="kosambi", n.perm=1000)
summary(scan.cim_TGMR_SM22)
summary(perm_TGMR_SM22)

plot(scan.cim_TGMR_SM22, ylab = "LOD", main = "QTL profile PGMR (SM22)")
summary(scan.cim_TGMR_SM22, perms =  perm_TGMR_SM22, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =21 ,infer = F,  marker = "S6_42122367", ylab = "TGMR (SM22)")
??plotPXG

qtl <- makeqtl(phenogeno_cross, chr=6, pos=118, what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=27, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results = scan.cim_TGMR_SM22, chr = 6, expandtomarkers = T)


# > summary(scan.cim_TGMR_SM22)
#              chr   pos    lod
# S1_61874941    1 100.6  4.767
# S2_57728127    2  97.5  4.403
# S3_11489823    3  14.0  1.441
# S4_63211680    4  58.3  1.755
# S5_66360856    5  39.1  2.347
# S6_42122367    6 118.5 12.005
# S7_62303916    7 277.4  1.457
# S8_6118586     8 106.2  2.262
# S9_49940810    9 114.6  1.628
# S10_55493385  10 241.6  0.741
# > summary(scan.cim_TGMR_SM22, perms =  perm_TGMR_SM22, alpha=0.1, pvalues=TRUE)
#             chr   pos   lod  pval
# S1_61874941   1 100.6  4.77 0.012
# S2_57728127   2  97.5  4.40 0.021
# S6_42122367   6 118.5 12.00 0.000
# > plotPXG(phenogeno_cross, pheno.col =25 ,infer = F,  marker = "S1_61874941", ylab = "TGMR (SM22)")
# > qtl <- makeqtl(phenogeno_cross, chr=1, pos=100.6, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=25, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 286 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS         MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  25.03393 25.0339340 8.304429 12.51631 6.245084e-10 7.409143e-10
# Error 284 174.97649  0.6161144                                            
# Total 285 200.01042                                                       
# 
# 
# Estimated effects:
# -----------------
#               est      SE      t
# Intercept 3.16604 0.04644 68.176
# 1@100.6   0.30607 0.04802  6.374
# 
# > lodint(results = scan.cim_TGMR_SM22, chr = 1, expandtomarkers = T)
#             chr       pos       lod
# S1_62943063   1  94.89423 0.1407407
# S1_61874941   1 100.55034 4.7671449
# S1_15315301   1 153.45119 2.8841936
# > plotPXG(phenogeno_cross, pheno.col =25 ,infer = F,  marker = "S2_57728127", ylab = "TGMR (SM22)")
# > plotPXG(phenogeno_cross, pheno.col =25 ,infer = F,  marker = "S2_57728127 ", ylab = "TGMR (SM22)")
# Error in plotPXG(phenogeno_cross, pheno.col = 25, infer = F, marker = "S2_57728127 ",  : 
#   Marker S2_57728127  not found
# > plotPXG(phenogeno_cross, pheno.col =25 ,infer = F,  marker = "S2_57728127", ylab = "TGMR (SM22)")
# > qtl <- makeqtl(phenogeno_cross, chr=2, pos=97.5, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=25, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 286 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS         MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  16.32362 16.3236216 5.287379 8.161385 8.035681e-07 8.973026e-07
# Error 284 183.68680  0.6467845                                            
# Total 285 200.01042                                                       
# 
# 
# Estimated effects:
# -----------------
#               est      SE      t
# Intercept 3.14419 0.04797 65.545
# 2@97.5    0.24392 0.04855  5.024
# 
# > lodint(results = scan.cim_TGMR_SM22, chr = 2, expandtomarkers = T)
#             chr      pos       lod
# S2_14749096   2 87.04723 0.1263962
# S2_57728127   2 97.54109 4.4030516
# S2_57981108   2 99.53119 2.7456908
# > plotPXG(phenogeno_cross, pheno.col =25 ,infer = F,  marker = "S6_42122367", ylab = "TGMR (SM22)")
# > qtl <- makeqtl(phenogeno_cross, chr=6, pos=118.5, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=25, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 286 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS         MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  35.52652 35.5265160 12.14488 17.76233  7.51621e-14 9.625634e-14
# Error 284 164.48391  0.5791687                                            
# Total 285 200.01042                                                       
# 
# 
# Estimated effects:
# -----------------
#               est      SE      t
# Intercept 3.23386 0.04561 70.907
# 6@118.5   0.36037 0.04601  7.832
# 
# > lodint(results = scan.cim_TGMR_SM22, chr = 6, expandtomarkers = T)
#             chr      pos        lod
# S6_42793946   6 109.6862  0.1185882
# S6_42122367   6 118.4694 12.0047269
# S6_6110410    6 126.0970  0.5963012
# >

#SM23
names(GrGn_qtl$pheno)

scan.cim_TGMR_SM23<- cim(phenogeno_cross, pheno.col=23, map.function="kosambi")
perm_TGMR_SM23 <- cim(phenogeno_cross, pheno.col= 23, map.function="kosambi", n.perm=1000)
summary(scan.cim_TGMR_SM23)
summary(perm_TGMR_SM23)

plot(scan.cim_TGMR_SM23, ylab = "LOD", main = "QTL profile PGMR (SM22)")
summary(scan.cim_TGMR_SM23, perms =  perm_TGMR_SM23, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =23 ,infer = F,  marker = "S2_57680557", ylab = "TGMR (SM23)")
??plotPXG


qtl <- makeqtl(phenogeno_cross, chr=7, pos=264.9, what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=23, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results = scan.cim_TGMR_SM23, chr = 7, expandtomarkers = T)


# > summary(scan.cim_TGMR_SM23, perms =  perm_TGMR_SM23, alpha=0.1, pvalues=TRUE)
#             chr   pos   lod  pval
# S1_17391798   1 116.7  3.44 0.080
# S2_57680557   2  97.8 12.00 0.000
# S3_4555653    3  38.1  5.25 0.003
# S7_59404098   7 264.9  3.88 0.039
# > qtl <- makeqtl(phenogeno_cross, chr=1, pos=116.7, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=23, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 8 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 278 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df         SS        MS      LOD     %var Pvalue(Chi2)  Pvalue(F)
# Model   1   2.883241 2.8832410 1.035995 1.701521   0.02894412 0.02967451
# Error 276 166.567545 0.6035056                                          
# Total 277 169.450786                                                    
# 
# 
# Estimated effects:
# -----------------
#                est       SE      t
# Intercept  3.76578  0.09231 40.795
# 1@116.7   -0.20190  0.09237 -2.186
# 
# > lodint(results = scan.cim_TGMR_SM23, chr = 1, expandtomarkers = T)
#             chr      pos      lod
# S1_64013631   1  91.1871 1.565268
# S1_17391798   1 116.6895 3.436284
# S1_59458791   1 116.8561 0.011898
# > qtl <- makeqtl(phenogeno_cross, chr=2, pos=97.8, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=23, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 8 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 278 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS         MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  24.71317 24.7131677 9.516241 14.58427 3.592893e-11 4.391942e-11
# Error 276 144.73762  0.5244117                                            
# Total 277 169.45079                                                       
# 
# 
# Estimated effects:
# -----------------
#               est      SE      t
# Intercept 3.89698 0.04388 88.807
# 2@97.8    0.30401 0.04428  6.865
# 
# > lodint(results = scan.cim_TGMR_SM23, chr = 2, expandtomarkers = T)
#             chr       pos        lod
# S2_56846223   2  94.64631  8.3819608
# S2_57680557   2  97.75126 12.0000191
# S2_58752330   2 104.32705  0.7931163
# > qtl <- makeqtl(phenogeno_cross, chr=3, pos=38.1, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=23, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 8 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 278 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS         MS     LOD    %var Pvalue(Chi2)    Pvalue(F)
# Model   1  11.73291 11.7329075 4.33161 6.92408 7.958545e-06 8.740411e-06
# Error 276 157.71788  0.5714416                                          
# Total 277 169.45079                                                     
# 
# 
# Estimated effects:
# -----------------
#               est      SE      t
# Intercept 3.97276 0.04591 86.530
# 3@38.1    0.21437 0.04731  4.531
# 
# > lodint(results = scan.cim_TGMR_SM23, chr = 23, expandtomarkers = T)
# Error in lodint(results = scan.cim_TGMR_SM23, chr = 23, expandtomarkers = T) : 
#   Chromosome 23 not found.
# 
# > lodint(results = scan.cim_TGMR_SM23, chr = 3, expandtomarkers = T)
#             chr      pos       lod
# S3_52787725   3 33.01659 0.1270354
# S3_4555653    3 38.08470 5.2467950
# S3_63217687   3 97.44606 0.1424754
# > qtl <- makeqtl(phenogeno_cross, chr=7, pos=264.9, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=23, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 8 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 278 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df         SS        MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1   6.726548 6.7265479 2.445196 3.969617 0.0007917194 0.0008360839
# Error 276 162.724238 0.5895806                                            
# Total 277 169.450786                                                      
# 
# 
# Estimated effects:
# -----------------
#                est       SE      t
# Intercept  3.90419  0.04725 82.621
# 7@264.9   -0.16196  0.04795 -3.378
# 
# > lodint(results = scan.cim_TGMR_SM23, chr = 7, expandtomarkers = T)
#             chr      pos      lod
# S7_56509169   7   0.0000 1.870454
# S7_59404098   7 264.9164 3.881728
# S7_60598717   7 272.8259 2.217134

#Boucle TGMR plot------

qtl_info_TGMR <- data.frame(
  marker = c("S4_59554213", 
             "S6_1839684", "S7_54708458", 
             "S7_50318042","S10_7557386",
             "S1_61874941", "S2_57728127", "S6_42122367",
             "S1_17391798", "S2_57680557", "S3_4555653"),
  
  pheno_col = c(19, 
                20, 20, 
                22, 22,
                21, 21, 21,
                23, 23, 23),
  
  ylab = c("TGMR (BY22)", 
           "TGMR (SD22)", "TGMR (SD22)", 
           "TGMR (SD23)", "TGMR (SD23)",
           "TGMR (SM22)", "TGMR (SM22)", "TGMR (SM22)",
           "TGMR (SM23)", "TGMR (SM23)", "TGMR (SM23)")
)


# -----------------------------
# 2. Fichier haute résolution
# -----------------------------
png("Output/plot_TGMR.png", 
    width = 2500, height = 3500, res = 300)

# -----------------------------
# 3. Paramètres graphiques
# -----------------------------
par(mfrow = c(4, 3),
    mar = c(4.5, 4, 2, 1),
    oma = c(1, 1, 1, 1))
# -----------------------------
# 4. Boucle de tracé
# -----------------------------
# Boucle de dessin
for(i in 1:nrow(qtl_info_TGMR)){
  plotPXG(
    phenogeno_cross,
    pheno.col = qtl_info_TGMR$pheno_col[i],
    marker    = qtl_info_TGMR$marker[i],
    infer     = FALSE,
    ylab      = qtl_info_TGMR$ylab[i],
    cex.lab   = 1.3,  # taille des labels X/Y
    cex.axis  = 0.9  # titre (optionnel)
  )
}

dev.off()

#YLD-----

#BY22
names(GrGn_qtl$pheno)

scan.cim_YLD_BY22<- cim(phenogeno_cross, pheno.col=14, map.function="kosambi")
perm_YLD_BY22 <- cim(phenogeno_cross, pheno.col= 14, map.function="kosambi", n.perm=1000)
summary(scan.cim_YLD_BY22)
summary(perm_YLD_BY22)

plot(scan.cim_YLD_BY22, ylab = "LOD", main = "QTL profile PGMR (SM22)")
summary(scan.cim_YLD_BY22, perms =  perm_YLD_BY22, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =16 ,infer = F,  marker = "S1_62943063", ylab = "YLD (BY22)")
??plotPXG


qtl <- makeqtl(phenogeno_cross, chr=1, pos=94.9, what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=16, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results = scan.cim_YLD_BY22, chr = 1, expandtomarkers = T)


# > summary(scan.cim_YLD_BY22)
# chr    pos   lod
# S1_62943063   1  94.89 4.220
# S2_57981108   2  99.53 1.521
# S3_12619687   3   9.85 0.763
# S4_61295445   4  65.69 0.363
# S5_62754800   5  23.59 0.462
# S6_40634627   6 118.91 2.185
# S7_1838966    7   0.00 0.826
# S8_54835462   8 105.18 1.192
# S9_9992707    9 174.52 0.839
# snpSB00146   10 211.58 0.550
# > summary(perm_YLD_BY22)
# LOD thresholds (1000 permutations)
# [,1]
# 5%  3.97
# 10% 3.53

# > summary(scan.cim_YLD_BY22, perms =  perm_YLD_BY22, alpha=0.1, pvalues=TRUE)
# chr  pos  lod  pval
# S1_62943063   1 94.9 4.22 0.033
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 283 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df           SS          MS          LOD         %var Pvalue(Chi2) Pvalue(F)
# Model   1 1.271381e+02    127.1381 5.744108e-05 9.347202e-05    0.9870236  0.987081
# Error 281 1.360171e+08 484046.6536                                                 
# Total 282 1.360172e+08                                                             
# 
# 
# Estimated effects:
# -----------------
#                est       SE      t
# Intercept 868.2686  42.2943 20.529
# 7@264.9     0.6954  42.9060  0.016
# 
# > lodint(results = scan.cim_YLD_BY22, chr = 1, expandtomarkers = T)
#             chr      pos      lod
# S1_64013631   1 91.18710 1.692929
# S1_62943063   1 94.89423 4.219967
# S1_61976676   1 99.16546 2.177315
# >
# 
# 

#YLD SD22
names(GrGn_qtl$pheno)

scan.cim_YLD_SD22<- cim(phenogeno_cross, pheno.col=16, map.function="kosambi")
perm_YLD_SD22 <- cim(phenogeno_cross, pheno.col= 16, map.function="kosambi", n.perm=1000)
summary(scan.cim_YLD_SD22)
summary(perm_YLD_SD22)

plot(scan.cim_YLD_SD22, ylab = "LOD", main = "QTL profile PGMR (SM22)")
summary(scan.cim_YLD_SD22, perms =  perm_YLD_SD22, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =16 ,infer = F,  marker = "S5_64396901", ylab = "YLD (SD22)")
??plotPXG

qtl <- makeqtl(phenogeno_cross, chr=5, pos=26.9, what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=16, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results = scan.cim_YLD_SD22, chr = 5, expandtomarkers = T)

# > scan.cim_YLD_SD22<- cim(phenogeno_cross, pheno.col=16, map.function="kosambi")
# > perm_YLD_SD22 <- cim(phenogeno_cross, pheno.col= 16, map.function="kosambi", n.perm=1000)
# > summary(scan.cim_YLD_SD22)
#             chr   pos   lod
# S1_16745318   1 149.0 2.378
# S2_2834191    2  15.5 1.166
# S3_69405187   3 328.2 0.921
# S4_59554213   4  67.2 1.179
# S5_64396901   5  26.9 4.352
# S6_59030350   6  44.6 1.715
# S7_6584892    7  35.6 1.180
# S8_14252601   8  95.8 0.722
# S9_49940810   9 114.6 1.971
# S10_1744527  10  54.9 0.812
# > summary(perm_YLD_SD22)
# LOD thresholds (1000 permutations)
#     [,1]
# 5%  4.13
# 10% 3.52
# > summary(scan.cim_YLD_SD22, perms =  perm_YLD_SD22, alpha=0.1, pvalues=TRUE)
#             chr  pos  lod pval
# S5_64396901   5 26.9 4.35 0.03
# > plotPXG(phenogeno_cross, pheno.col =18 ,infer = F,  marker = "S5_64396901", ylab = "YLD (SD22)")
# > qtl <- makeqtl(phenogeno_cross, chr=5, pos=26.9, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=16, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 4 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 282 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df       SS        MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  6203982 6203981.9 6.265604 9.725915 7.803765e-08 8.903486e-08
# Error 280 57584174  205657.8                                            
# Total 281 63788156                                                      
# 
# 
# Estimated effects:
# -----------------
#               est      SE      t
# Intercept  901.49   27.01 33.381
# 5@26.9    -151.43   27.57 -5.492
# 
# > lodint(results = scan.cim_YLD_SD22, chr = 5, expandtomarkers = T)
#             chr      pos        lod
# S5_59574508   5 19.18057 2.50486304
# S5_64396901   5 26.94659 4.35201425
# S5_64829267   5 29.58272 0.02373995
# >


#SM22
names(GrGn_qtl$pheno)

scan.cim_YLD_SM22<- cim(phenogeno_cross, pheno.col=15, map.function="kosambi")
perm_YLD_SM22 <- cim(phenogeno_cross, pheno.col= 15, map.function="kosambi", n.perm=1000)
summary(scan.cim_YLD_SM22)
summary(perm_YLD_SM22)

summary(scan.cim_YLD_SM22, perms =  perm_YLD_SM22, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =17 ,infer = F,  marker = "S6_47315279", ylab = "YLD (SM22)")
??plotPXG

qtl <- makeqtl(phenogeno_cross, chr=6, pos=6.6, what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=17, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results = scan.cim_YLD_SM22, chr = 6, expandtomarkers = T)

# > #SM22
# > names(GrGn_qtl$pheno)
#  [1] "FLA_BY22"  "FLA_SM22"  "FLA_SD22"  "FLA_BY23"  "FLA_SM23"  "FLA_SD23" 
#  [7] "Intcp_FLA" "slope_FLA" "PGMR_BY22" "PGMR_SM22" "PGMR_SD22" "PGMR_SM23"
# [13] "PGMR_SD23" "PGMR_SM"   "PGMR_SD"   "YLD_BY22"  "YLD_SM22"  "YLD_SD22" 
# [19] "YLD_SM23"  "YLD_SD23"  "YLD_SM"    "YLD_SD"    "TGMR_BY22" "TGMR_SD22"
# [25] "TGMR_SM22" "TGMR_SD23" "TGMR_SM23" "TGMR_SM"   "TGMR_SD"  
# > summary(scan.cim_YLD_SM22, perms =  perm_YLD_SM22, alpha=0.1, pvalues=TRUE)
#             chr   pos  lod  pval
# S1_64013631   1 91.19 4.28 0.035
# S3_71275540   3 99.28 4.25 0.037
# S6_47315279   6  6.62 6.53 0.006
# > qtl <- makeqtl(phenogeno_cross, chr=3, pos=99.28, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=17, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 1 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 285 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS       MS      LOD     %var Pvalue(Chi2)   Pvalue(F)
# Model   1  10036577 10036577 5.419747 8.384979  5.85668e-07 6.55988e-07
# Error 283 109660530   387493                                           
# Total 284 119697107                                                    
# 
# 
# Estimated effects:
# -----------------
#               est      SE      t
# Intercept 1514.96   68.11 22.243
# 3@99.3    -349.07   68.59 -5.089
# 
# > lodint(results = scan.cim_YLD_SM22, chr = 3, expandtomarkers = T)
#             chr        pos        lod
# S3_13359317   3   9.097768 2.62003185
# S3_71275540   3  99.278327 4.25401072
# S3_55432695   3 106.803546 0.01549337
# > qtl <- makeqtl(phenogeno_cross, chr=1, pos=99.19, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=17, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 1 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 285 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS        MS     LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1   5485701 5485700.8 2.90332 4.582985 0.0002556428 0.0002721277
# Error 283 114211406  403573.9                                           
# Total 284 119697107                                                     
# 
# 
# Estimated effects:
# -----------------
#               est      SE      t
# Intercept 1815.65   37.71 48.143
# 1@99.2    -143.54   38.93 -3.687
# 
# > lodint(results = scan.cim_YLD_SM22, chr = 1, expandtomarkers = T)
#             chr      pos       lod
# S1_65180762   1 84.91579 0.1468243
# S1_64013631   1 91.18710 4.2840540
# S1_62992453   1 94.89422 2.6884532
# > qtl <- makeqtl(phenogeno_cross, chr=6, pos=6.6, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=17, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 1 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 285 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS       MS      LOD     %var Pvalue(Chi2)   Pvalue(F)
# Model   1  13629834 13629834 7.481571 11.38694 4.364661e-09 5.09619e-09
# Error 283 106067273   374796                                           
# Total 284 119697107                                                    
# 
# 
# Estimated effects:
# -----------------
#               est      SE     t
# Intercept 1767.16   36.84 47.96
# 6@6.6      223.56   37.07  6.03
# 
# > lodint(results = scan.cim_YLD_SM22, chr = 1, expandtomarkers = T)
#             chr      pos       lod
# S1_65180762   1 84.91579 0.1468243
# S1_64013631   1 91.18710 4.2840540
# S1_62992453   1 94.89422 2.6884532
# > lodint(results = scan.cim_YLD_SM22, chr = 6, expandtomarkers = T)
#             chr      pos      lod
# S6_47055810   6 5.314446 4.568770
# S6_47315279   6 6.624349 6.528516
# S6_46064034   6 8.265517 3.765710
# > lodint(results = scan.cim_YLD_SM22, chr = 6, expandtomarkers = T)
#             chr      pos      lod
# S6_47055810   6 5.314446 4.568770
# S6_47315279   6 6.624349 6.528516
# S6_46064034   6 8.265517 3.765710
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 285 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS       MS      LOD     %var Pvalue(Chi2)   Pvalue(F)
# Model   1  13629834 13629834 7.481571 11.38694 4.364661e-09 5.09619e-09
# Error 283 106067273   374796                                           
# Total 284 119697107                                                    
# 
# 
# Estimated effects:
# -----------------
#               est      SE     t
# Intercept 1767.16   36.84 47.96
# 6@6.6      223.56   37.07  6.03

  #SM23
names(GrGn_qtl$pheno)

scan.cim_YLD_SM23<- cim(phenogeno_cross, pheno.col=17, map.function="kosambi")
perm_YLD_SM23 <- cim(phenogeno_cross, pheno.col= 17, map.function="kosambi", n.perm=1000)
summary(scan.cim_YLD_SM23)
summary(perm_YLD_SM23)

plot(scan.cim_YLD_SD23, ylab = "LOD", main = "QTL profile PGMR (SM22)")
summary(scan.cim_YLD_SM23, perms =  perm_YLD_SM23, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =19 ,infer = F,  marker = "Sbv3.1_06_40312464", ylab = "YLD (SM23)")
??plotPXG

qtl <- makeqtl(phenogeno_cross, chr=6, pos=119, what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=19, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results = scan.cim_YLD_SM23, chr = 6, expandtomarkers = T)

                   # chr   pos    lod
# S1_62992453          1  94.9  2.967
# S2_14749096          2  87.0  1.290
# S3_66095714          3 329.5  3.393
# S4_68260182          4 269.0  1.018
# S5_66360856          5  39.1  1.555
# Sbv3.1_06_40312464   6 119.3 14.026
# S7_62303916          7 277.4  1.363
# S8_4920840           8 102.3  2.190
# S9_51062242          9 199.1  0.421
# S10_8964377         10  88.9  2.886
# > summary(perm_YLD_SM23)
# LOD thresholds (1000 permutations)
#     [,1]
# 5%  3.83
# 10% 3.49
# > summary(scan.cim_YLD_SM23, perms =  perm_YLD_SM23, alpha=0.1, pvalues=TRUE)
#                    chr pos lod pval
# Sbv3.1_06_40312464   6 119  14    0
# > plotPXG(phenogeno_cross, pheno.col =17 ,infer = F,  marker = "Sbv3.1_06_40312464", ylab = "YLD (SM23)")
# > plotPXG(phenogeno_cross, pheno.col =19 ,infer = F,  marker = "Sbv3.1_06_40312464", ylab = "YLD (SM23)")
# > qtl <- makeqtl(phenogeno_cross, chr=6, pos=119, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=19, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 10 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 276 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS         MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  54010226 54010226.0 12.96156 19.44812 1.110223e-14 1.465494e-14
# Error 274 223704141   816438.5                                            
# Total 275 277714367                                                       
# 
# 
# Estimated effects:
# -----------------
#               est      SE      t
# Intercept 2232.05   55.26 40.390
# 6@118.9    453.79   55.79  8.133
# 
# > lodint(results = scan.cim_YLD_SM23, chr = 6, expandtomarkers = T)
#                    chr      pos         lod
# S6_42793946          6 109.6862  0.00548005
# Sbv3.1_06_40312464   6 119.3410 14.02599398
# S6_6110410           6 126.0970  0.02305628

#Boucle YLD plot------

  # 1. Tableau : marqueur + colonne + label
  # -----------------------------
qtl_info_YLD <- data.frame(
  markers <- c("S1_62943063", 
               "S5_64396901", 
               "S1_64013631", "S3_71275540", "S6_47315279",
               "Sbv3.1_06_40312464"),
  pheno_col = c(14,
                16, 
                15, 15,15,
                17),
  ylab = c("YLD (BY22)",
           "YLD (SD22)", 
           "YLD (SM22)", "YLD (SM22)","YLD (SM22)",
           "YLD (SM23)")
)

# -----------------------------
# 2. Fichier haute résolution
# -----------------------------
png("Output/plot_YLD.png", 
    width = 3000, height = 2500, res = 300)

# -----------------------------
# 3. Paramètres graphiques
# -----------------------------
par(mfrow = c(2, 3),
    mar = c(4.5, 4, 2, 1),
    oma = c(1, 1, 1, 1))
# -----------------------------
# 4. Boucle de tracé
# -----------------------------
# Boucle de dessin
for(i in 1:nrow(qtl_info_YLD)){
  plotPXG(
    phenogeno_cross,
    pheno.col = qtl_info_YLD$pheno_col[i],
    marker    = qtl_info_YLD$marker[i],
    infer     = FALSE,
    ylab      = qtl_info_YLD$ylab[i],
    cex.lab   = 1.3,  # taille des labels X/Y
    cex.axis  = 0.9  # titre (optionnel)
  )
}

dev.off()

#PDW-----
#BY22
names(GrGn_qtl$pheno)

scan.cim_PDW_SM23<- cim(phenogeno_cross, pheno.col=28, map.function="kosambi")
perm_PDW_SM23 <- cim(phenogeno_cross, pheno.col= 28, map.function="kosambi", n.perm=1000)
summary(scan.cim_PDW_SM23)
summary(perm_PDW_SM23)

plot(scan.cim_PDW_BY22, ylab = "LOD", main = "QTL profile PGMR (SM22)")
summary(scan.cim_PDW_SM23, perms =  perm_PDW_SM23, alpha=0.1, pvalues=TRUE)

plotPXG(phenogeno_cross, pheno.col =28 ,infer = F,  marker = "Sbv3.1_06_40312464", ylab = "PDW (SM23)")
??plotPXG

qtl <- makeqtl(phenogeno_cross, chr=6, pos=119.3, what=c("prob")) 
fitqtl <- fitqtl(GrGn_qtl, pheno.col=28, 
                 formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
summary(fitqtl)

lodint(results = scan.cim_PDW_SM23, chr = 6, expandtomarkers = T)

# > summary(scan.cim_PDW_SM23, perms =  perm_PDW_SM23, alpha=0.1, pvalues=TRUE)
#                    chr   pos  lod  pval
# S3_71605396          3  99.6 3.75 0.047
# Sbv3.1_06_40312464   6 119.3 8.38 0.000
# > qtl <- makeqtl(phenogeno_cross, chr=3, pos=99.6, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=28, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 10 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 276 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df        SS        MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
# Model   1  5563.786 5563.7856 6.824622 10.76274 2.069165e-08 2.395222e-08
# Error 274 46131.087  168.3616                                            
# Total 275 51694.873                                                      
# 
# 
# Estimated effects:
# -----------------
#              est     SE      t
# Intercept 39.953  1.451 27.531
# 3@99.6    -8.421  1.465 -5.749
# 
# > lodint(results = scan.cim_PDW_SM23, chr = 3, expandtomarkers = T)
#             chr       pos       lod
# S3_3441687    3  39.22099 1.7520284
# S3_71605396   3  99.59310 3.7478711
# S3_55432695   3 106.80355 0.2552379
# > summary(scan.cim_PDW_SM23, perms =  perm_PDW_SM23, alpha=0.1, pvalues=TRUE)
#                    chr   pos  lod  pval
# S3_71605396          3  99.6 3.75 0.047
# Sbv3.1_06_40312464   6 119.3 8.38 0.000
# > qtl <- makeqtl(phenogeno_cross, chr=6, pos=119.3, what=c("prob")) 
# > fitqtl <- fitqtl(GrGn_qtl, pheno.col=28, 
# +                  formula = y~Q1, qtl= qtl, method = "hk", get.ests = T)
# Warning message:
# In fitqtlengine(pheno = pheno, qtl = qtl, covar = covar, formula = formula,  :
#   Dropping 10 individuals with missing phenotypes.
# 
# > summary(fitqtl)
# 
# 		fitqtl summary
# 
# Method: Haley-Knott regression 
# Model:  normal phenotype
# Number of observations : 276 
# 
# Full model result
# ----------------------------------  
# Model formula: y ~ Q1 
# 
#        df       SS         MS      LOD    %var Pvalue(Chi2)    Pvalue(F)
# Model   1 10061.73 10061.7344 12.97315 19.4637 1.076916e-14 1.421085e-14
# Error 274 41633.14   151.9458                                           
# Total 275 51694.87                                                      
# 
# 
# Estimated effects:
# -----------------
#               est      SE      t
# Intercept 48.0810  0.7541 63.758
# 6@119.3    6.1942  0.7612  8.138
# 
# > lodint(results = scan.cim_PDW_SM23, chr = 6, expandtomarkers = T)
#                    chr      pos        lod
# S6_42793946          6 109.6862 0.22370050
# Sbv3.1_06_40312464   6 119.3410 8.38283144
# S6_6110410           6 126.0970 0.06812802
# >

#Interaction plot-----

png("Output/interaction_plot_FR.png", width = 2000, height = 3000, res = 300)

# Paramètres graphiques : 3 lignes, 1 colonne
par(mfrow = c(3,1),
    mar = c(4.5, 4, 2, 1),
    oma = c(1,1,1,1),
    cex.lab = 1.3,     # taille des labels X/Y
    cex.axis = 1.2)

effectplot(phenogeno, mname2="6@119.3", mname1="6@132")
title(main = "", cex.main = 1.3)

effectplot(phenogeno, mname2="6@119.3", mname1="6@118")
effectplot(phenogeno, mname2="6@132", mname1="6@118")

# Fin du fichier
dev.off()

#PLot chr 6-----

# Subset du chromosome 6
cim_chr6_FLA <- scan.cim_FLA_BY22[scan.cim_FLA_BY22$chr == 6, ]

# Plot
plot(cim_chr6_FLA, perm = perm_FLA_B22, col = "red")




# Fichier PNG
png("Output/interactionFR.png", width = 4000, height = 2300, res = 300)

# 4 lignes, 1 colonne
par(mfrow = c(2,2),
    mar = c(4.5, 6, 2, 1),
    oma = c(1,1,1,1),
    fg= "black")

cim_chr6_FLA <- scan.cim_FLA_BY22[scan.cim_FLA_BY22$chr == 6, ]
plot(cim_chr6_FLA, perm = perm_FLA_B22, col = "red", main = "LOD Chromosome 6",
     ylab= "LOD score", cex.main=2, cex.lab= 1.9, cex.axis= 1.9)

par(cex.lab = 1.9,      # taille du nom des axes
    cex.axis =1.9)
effectplot(phenogeno, mname2="6@119.3", mname1="6@132", main= "")
title(main = "Interaction Ma1 et qFLA6.2", cex.main = 1.9)

par(cex.lab = 1.9,      # taille du nom des axes
    cex.axis = 1.9)
effectplot(phenogeno, mname2="6@119.3", mname1="6@118", main= "")
title(main = "Interaction Ma1 et qFLA6.3", cex.main = 1.9)

par(cex.lab = 1.9,      # taille du nom des axes
    cex.axis = 1.9)
effectplot(phenogeno, mname2="6@132", mname1="6@118", main= "")
title(main = "Interaction qFLA6.2 et qFLA6.3", cex.main = 1.9)

dev.off()



cim_chr6_FLA <- scan.cim_FLA_BY22[scan.cim_FLA_BY22$chr == 6, ]

plot(cim_chr6_FLA, col = "red", main = "LOD Chromosome 6",
     ylab= "LOD score", cex.main=2, cex.lab= 1.9, cex.axis= 1.9)



# scanone sans covariable-----
out_yld <- scanone(cross, pheno.col="YLD", method="hk")
# scanone avec FLA comme covariate (additive)
out_yld_cond <- scanone(phenogeno_cross, pheno.col="YLD_SM23", addcovar=phenogeno_cross$pheno$FLA_BY22, method="hk")

summary(out_yld_cond)
# plot superposé
plot(out_yld, main="Scanone YLD: sans et avec FLA", ylim=c(0, max(out_yld$lod, out_yld_cond$lod, na.rm=TRUE)))
plot(out_yld_cond, add=TRUE, col="red")
legend("topright", legend=c("sans FLA","avec FLA"), col=c("black","red"), lwd=2)



