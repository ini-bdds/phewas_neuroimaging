# Takes as input the genotype data and outputs two columns: subjID and the genotype of interest, indicated by index_genotype, and only the data for which there are a higher than min_genetic percentage of subjects

genodata_handling <- function(geno_data,min_genetic,index_genotype) {
subj_col = which(names(geno_data)=="subjID")
aa=geno_data[,c(subj_col,index_genotype)]

uniq_geno = dim(unique(aa[2]))[1]
wkd <- unique(aa[,2])
wkd_pl = NULL
for (ww in 1:uniq_geno) {
    wkd_pl[ww]= length(aa[which (aa[2] == as.character(wkd[ww])),1])
}
wkd_tot = sum(wkd_pl)
geno_groups = 0
geno_excl = NULL
gexc = 1

for (ww in 1:uniq_geno) {
    # Exclude genotypes with less than % as specified by min_genetic
    if ((wkd_pl[ww]/wkd_tot * 100) > min_genetic) {
       geno_groups = geno_groups + 1
       }
     else {
    	 geno_excl[gexc]=as.character(wkd[ww])
	 gexc = gexc + 1
    }
}    
if (geno_groups < 2) {
   cat('WARNING: not enough data from the genotype in column ', 2, ' of ',genofile,' with at least ', min_genetic, '% of subjects in at least 2 groups \n')
   aa=NULL
   return(aa)
}

if (gexc > 1) {
   for (gg in 1:(gexc-1)) {
        aa<-aa[!(aa[,2]==geno_excl[gg]),]
   }
}	

return(aa)
}
