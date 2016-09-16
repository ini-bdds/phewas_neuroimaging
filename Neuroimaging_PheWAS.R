# Kristi Clark (kclark@ini.usc.edu)
# Sept 3, 2016

Neuroimaging_PheWAS <- function (genofile, imagfile, covfile, out_prefix, min_genetic, max_order, outlier_criteria) {

# Load libraries
library(heplots)
library(fgui)
library(nlstools)

# Load my functions
source("/ifs/loni/faculty/kclark/candl/kristi/PheWAS/genodata_handling.R")
source("/ifs/loni/faculty/kclark/candl/kristi/PheWAS/PheWAS_models.R")

# Read inputs
geno_data=read.csv(genofile)
imag_data = read.csv (imagfile)
cov_data = read.csv (covfile)

# Check and store variable names and indices
geno_names=names(geno_data)
my_genonames=NULL
index_genotype=NULL
a=1
for (aa in 1:length(geno_names)) {
    if (substr(geno_names[aa],1,2)=="rs") {
       index_genotype[a] = aa
       my_genonames[a] = geno_names[aa] 
       a=a+1
    }
}
num_geno = length(my_genonames)
if (num_geno == 0) {
   stop('ERROR: missing genotype data. looking for a column in ',genofile,' that starts with rs\n')
}

imag_names=names(imag_data)
my_imagnames=NULL
index_imag=NULL
index_value = 0
index_imag_subjID = 0
a=1
for (aa in 1:length(imag_names)) {
    if ((imag_names[aa]!="subjID") & (imag_names[aa]!="value")) {
       my_imagnames[a]=imag_names[aa]
       index_imag[a]=aa
       a=a+1
    }
    if (imag_names[aa]=="value") {
       index_value = aa
    }    
    if (imag_names[aa]=="subjID") {
       index_imag_subjID = aa
    }
}
if (index_value == 0) {
   stop('ERROR: missing imaging data. looking for a column in ',imagfile,' called value\n')
}
if (index_imag_subjID == 0) {
   stop('ERROR: looking for a column called subjID ',imagfile,' \n')
}
imag_variables=unique(imag_data[,c(my_imagnames)])
num_imag=dim(imag_variables)[1]

cov_names=names(cov_data)
my_covnames=NULL
index_cov=NULL
a=1
for (aa in 1:length(cov_names)) {
    if (cov_names[aa]!="subjID") {
       index_cov[a] = aa
       my_covnames[a] = cov_names[aa] 
       a = a+1
    }
}
num_cov=length(my_covnames)
pairs_cov=NULL
a=1
if (max_order > 1) {
   for (z in 1:num_cov) {
       for (y in 1:num_cov) {
           if (y<z) {
              pairs_cov[a]=paste(my_covnames[z],my_covnames[y],sep="*")
	      a=a+1
	  }
	}
    }	  
}
num_pairs=length(pairs_cov)

# Initialize output files
sorted_outINT=NULL
output_figureINT=NULL
sorted_outINT2=NULL
output_figureINT2=NULL
resultsALL=NULL

for (g in 1:length(num_geno)) {
    geno_subset_data=genodata_handling(geno_data,min_genetic,index_genotype[g])
    if (!is.null(geno_subset_data)) {
       # create output names for main effects
       outputfileALL = paste(out_prefix,my_genonames[g],"ALL.csv",sep="_")
       sorted_outME = paste(out_prefix,my_genonames[g],"sortedME.csv",sep="_")
       output_figureME = paste(out_prefix,my_genonames[g],"ME.png",sep="")

       for (i in 1:num_imag) {
       	   # Select a vector of data
	   k = 1
	   imag_subset_data = imag_data[which(imag_data[,c(names(imag_variables[k]))] == imag_variables[i,k]),]	   
	   for (k in 2:dim(imag_variables)[2]) {
	       imag_subset_data = imag_subset_data[which(imag_subset_data[,c(names(imag_variables[k]))] == imag_variables[i,k]),]
	   }
	   imag_subset_data=imag_subset_data[,c(index_imag_subjID,index_value)]

	   if ((!(is.nan(mean(imag_subset_data[,2])))) & ((mean(imag_subset_data[,2]))!=0)) {
	   # remove outliers from image data
	   imag_clean_subset = imag_subset_data[(imag_subset_data$value > quantile(imag_subset_data$value, 0.25) - outlier_criteria*IQR(imag_subset_data$value)) & (imag_subset_data$value < quantile(imag_subset_data$value, 0.75) + outlier_criteria*IQR(imag_subset_data$value)),]

       	   # Estimate the main effects, i.e. max_order 0
	   resultsALL=PheWAS_models_0(imag_clean_subset,geno_subset_data,cov_data,my_genonames,g,imag_variables,i,resultsALL)

           if (max_order > 0) {
       	      for (c in 1:num_cov) {
       	      	  sorted_outINT[c] = paste(out_prefix,my_genonames[g],"sortedINT",my_covnames[c],".csv",sep="")
		  output_figureINT[c] = paste(out_prefix,my_genonames[g],"INT",my_covnames[c],".png",sep="")
 	       }
	      resultsALL=PheWAS_models_1(imag_clean_subset,geno_subset_data,cov_data,my_genonames,g,imag_variables,i,resultsALL)
	   }

	   if (max_order > 1) {
   	      for (z in 1:num_cov) {
       	      	  for (y in 1:num_cov) {
           	      if (y<z) {
		      	 sorted_outINT2[c] = paste(out_prefix,my_genonames[g],"sortedINT_2ndord",(paste(my_covnames[z],my_covnames[y],sep="*")),".csv",sep="")
	  	    	 output_figureINT2[c] = paste(out_prefix,my_genonames[g],"INT_2ndord",(paste(my_covnames[z],my_covnames[y],sep="*")),".png",sep="")       	              }
	          }
	       }
	      resultsALL=PheWAS_models_2(imag_clean_subset,geno_subset_data,cov_data,my_genonames,g,imag_variables,i,resultsALL)
	   } # ends the max_order 2 covariates loop
	   } # ends conditional statement that values exist
	   write.csv(resultsALL,outputfileALL,quote=FALSE,row.names=FALSE)
	} # ends the image data loop
	save_resultsME <- resultsALL[order(resultsALL$MElog_pvals, decreasing=TRUE),]
	write.csv(save_resultsME,sorted_outME,quote=FALSE,row.names=FALSE)
	
	# save manhattan plot
	png(filename=output_figureME)
	geno_only=resultsALL[which(resultsALL$model=="genotype"),]
	plot(geno_only$MElog_pvals)
	abline(h=(-log10(0.05/dim(geno_only)[1])), col = "blue")
	dev.off()

    } # ends the is.null statement for genotype data
} # ends loop through genotypes

}