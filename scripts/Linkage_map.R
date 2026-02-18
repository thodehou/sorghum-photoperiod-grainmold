rm(list = ls())
setwd("/Users/thierry/Documents/ceraas/")
library(onemap)
library(vcfR)

vcfR.object <- read.vcfR("snp_filtered.vcf")

data_grinkan<- onemap_read_vcfR(vcfR.object = vcfR.object,
                                parent1 = "GRINKAN_202", 
                                parent2 = "Nganda_68", 
                                cross = "ri self") 


plot(data_grinkan, all= FALSE) 


marker_NA_filtered<- filter_missing(data_grinkan, threshold = 0.25) ##removing of missing data
marker_NA_filtered

##find redondant marker
#marker_bins <- find_bins(marker_NA_filtered, exact = FALSE)
#marker_bins

segrega_test2 <- test_segregation(marker_NA_filtered)# segregation with markers with any bins filtering
segrega_test2

sugges_alpha<-Bonferroni_alpha(segrega_test2) #suggestion of a alpha value that should be considered at 0.05 Bonferroni's correction


##recombination fractions----

LOD_sug <- suggest_lod(marker_NA_filtered)#calculate a suggested LOD score using all markers
LOD_sug

twopts_sug <- rf_2pts(marker_NA_filtered, LOD = LOD_sug, max.rf = 0.5)# by the suggested LOD
twopts_sug
print(twopts_sug, c(mrk1="snpSB00285", mrk2="snpSB00330"))# show an example



LGs<- group(mark_all)#forming the groups
LGs_upgma <- group_upgma(mark_all, expected.groups = 10, inter = F, 
                         comp.mat = T)#define a expected group number 
plot(LGs_upgma) #display the groups


##Estimate a genetic map----
maps<-vector("list", LGs_upgma$n.groups)
for(i in 1:LGs_upgma$n.groups)
  maps[[i]]<- make_seq(order_seq(input.seq= make_seq(LGs_upgma,i),twopt.alg =
                                   "rec"), "force")
draw_map(maps, names = FALSE, grid = TRUE, cex.mrk = 0.7)#Drow the map
draw_map2(maps,output="maps0.png")

file.out<-"map_grinkan.csv"
write_map(maps, file.out)#Extract the map into external file
