# Estimate the models for PheWAS

PheWAS_models_0 <- function (imag_subset_data,geno_subset_data,cov_data,my_genonames,g,imag_variables,i,resultsALL) {

SNP=my_genonames[g]
INTlog_pvals=0
INTeffect_sizes=0

# Merge genotype data with imaging data
# Merge genotype data with imaging data
geno_imag_data <- merge(x=imag_subset_data,y=geno_subset_data,by.x="subjID",by.y="subjID")

# Merge geno_imag data with covariate data
geno_cov_imag_data <- merge(x=geno_imag_data,y=cov_data,by.x="subjID",by.y="subjID")
geno_index=which(names(geno_cov_imag_data)==SNP)

# Simplest model, genotype only
model="genotype"
simple<-lm(geno_cov_imag_data$value ~ geno_cov_imag_data[,geno_index])
MElog_pvals=-log10(anova(simple)[1,5])
MEeffect_sizes=etasq(simple,anova=FALSE,partial=TRUE)[1,]
modelBIC=BIC(simple)
numsubjincl=dim(simple$model)[1]

addrow <- cbind(SNP,model,imag_variables[i,],MElog_pvals,MEeffect_sizes,INTlog_pvals,INTeffect_sizes,modelBIC,numsubjincl)	  
resultsALL = rbind(resultsALL,addrow)
return(resultsALL)
}

PheWAS_models_1 <- function (imag_subset_data,geno_subset_data,cov_data,my_genonames,g,imag_variables,i,resultsALL) {

SNP=my_genonames[g]
INTlog_pvals=0
INTeffect_sizes=0

# Merge genotype data with imaging data
geno_imag_data <- merge(x=imag_subset_data,y=geno_subset_data,by.x="subjID",by.y="subjID")

# Merge geno_imag data with covariate data
geno_cov_imag_data <- merge(x=geno_imag_data,y=cov_data,by.x="subjID",by.y="subjID")
geno_index=which(names(geno_cov_imag_data)==SNP)

my_covnames=NULL
index_cov=NULL
a=1
cov_names=names(cov_data)
for (aa in 1:length(cov_names)) {
    if (cov_names[aa]!="subjID") {
       index_cov[a] = aa
       my_covnames[a] = cov_names[aa] 
       a = a+1
    }
}

# Add each covariate one at a time to estimate the effect of the genotype but including the covariates
for (c in 1:length(my_covnames)) {
    model=paste("genotype",my_covnames[c],sep="+")
    cov_index=which(names(geno_cov_imag_data)==my_covnames[c])
    one_cov<-lm(geno_cov_imag_data$value ~ geno_cov_imag_data[,geno_index] + geno_cov_imag_data[,cov_index])
    MElog_pvals=-log10(anova(one_cov)[1,5])
    MEeffect_sizes=etasq(one_cov,anova=FALSE,partial=TRUE)[1,]
    modelBIC=BIC(one_cov)
    numsubjincl=dim(one_cov$model)[1]

    addrow <- cbind(SNP,model,imag_variables[i,],MElog_pvals,MEeffect_sizes,INTlog_pvals,INTeffect_sizes,modelBIC,numsubjincl)	  
    resultsALL = rbind(resultsALL,addrow)
}

# Add pairs of covariates, still estimating genotype effects only
for (z in 1:length(my_covnames)) {
    for (y in 1:length(my_covnames)) {
    	if (y<z) {
    	   model=paste("genotype",my_covnames[z],my_covnames[y],sep="+")
    	   cov_indexA=which(names(geno_cov_imag_data)==my_covnames[z])
	   cov_indexB=which(names(geno_cov_imag_data)==my_covnames[y])
	   
    	   two_cov<-lm(geno_cov_imag_data$value ~ geno_cov_imag_data[,geno_index] + geno_cov_imag_data[,cov_indexA] + geno_cov_imag_data[,cov_indexB])
    
            MElog_pvals=-log10(anova(two_cov)[1,5])
    	    MEeffect_sizes=etasq(two_cov,anova=FALSE,partial=TRUE)[1,]
    	    modelBIC=BIC(two_cov)
    	    numsubjincl=dim(two_cov$model)[1]

    	    addrow <- cbind(SNP,model,imag_variables[i,],MElog_pvals,MEeffect_sizes,INTlog_pvals,INTeffect_sizes,modelBIC,numsubjincl)	  
    	    resultsALL = rbind(resultsALL,addrow)
	}
     }
}

# Add triplets of covariates, still estimating genotype effects only
for (z in 1:length(my_covnames)) {
    for (y in 1:length(my_covnames)) {
    	for (x in 1:length(my_covnames)) {
    	    if ((x<y) & (y<z)) {
    	       model=paste("genotype",my_covnames[z],my_covnames[y],my_covnames[x],sep="+")
    	       cov_indexA=which(names(geno_cov_imag_data)==my_covnames[z])
	       cov_indexB=which(names(geno_cov_imag_data)==my_covnames[y])
    	       cov_indexC=which(names(geno_cov_imag_data)==my_covnames[x])
    	       three_cov<-lm(geno_cov_imag_data$value ~ geno_cov_imag_data[,geno_index] + geno_cov_imag_data[,cov_indexA] + geno_cov_imag_data[,cov_indexB] + geno_cov_imag_data[,cov_indexC])
    
		MElog_pvals=-log10(anova(three_cov)[1,5])
    	    	MEeffect_sizes=etasq(three_cov,anova=FALSE,partial=TRUE)[1,]
    	    	modelBIC=BIC(three_cov)
    	    	numsubjincl=dim(three_cov$model)[1]

    	    	addrow <- cbind(SNP,model,imag_variables[i,],MElog_pvals,MEeffect_sizes,INTlog_pvals,INTeffect_sizes,modelBIC,numsubjincl)
    	    	resultsALL = rbind(resultsALL,addrow)
	    }
	}
     }
}

return(resultsALL) 	   
}

PheWAS_models_2 <- function (imag_subset_data,geno_subset_data,cov_data,my_genonames,g,imag_variables,i,resultsALL) {

SNP=my_genonames[g]

# Merge genotype data with imaging data
geno_imag_data <- merge(x=imag_subset_data,y=geno_subset_data,by.x="subjID",by.y="subjID")

# Merge geno_imag data with covariate data
geno_cov_imag_data <- merge(x=geno_imag_data,y=cov_data,by.x="subjID",by.y="subjID")
geno_index=which(names(geno_cov_imag_data)==SNP)

my_covnames=NULL
index_cov=NULL
a=1
cov_names=names(cov_data)
for (aa in 1:length(cov_names)) {
    if (cov_names[aa]!="subjID") {
       index_cov[a] = aa
       my_covnames[a] = cov_names[aa] 
       a = a+1
    }
}

# Add each covariate one at a time to estimate the effect of the genotype AND the interaction of the genotype with the first covariate in the list
for (c in 1:length(my_covnames)) {
    interaction=paste("genotype",my_covnames[c],sep="*")
    model=paste("genotype",my_covnames[c],interaction,sep="+")
    cov_index=which(names(geno_cov_imag_data)==my_covnames[c])
    one_cov<-lm(geno_cov_imag_data$value ~ geno_cov_imag_data[,geno_index] + geno_cov_imag_data[,cov_index] + geno_cov_imag_data[,geno_index]*geno_cov_imag_data[,cov_index])
    MElog_pvals=-log10(anova(one_cov)[1,5])
    MEeffect_sizes=etasq(one_cov,anova=FALSE,partial=TRUE)[1,]
    INTlog_pvals=-log10(anova(one_cov)[3,5])
    INTeffect_sizes=etasq(one_cov,anova=FALSE,partial=TRUE)[3,]
    modelBIC=BIC(one_cov)
    numsubjincl=dim(one_cov$model)[1]

    addrow <- cbind(SNP,model,imag_variables[i,],MElog_pvals,MEeffect_sizes,INTlog_pvals,INTeffect_sizes,modelBIC,numsubjincl)	  
    resultsALL = rbind(resultsALL,addrow)
}

# Add pairs of covariates, estimating genotype and interaction with each covariate separately
for (z in 1:length(my_covnames)) {
    for (y in 1:length(my_covnames)) {
    	if (y<z) {
	   interaction=paste("genotype",my_covnames[z],sep="*")
    	   model=paste("genotype",my_covnames[z],my_covnames[y],interaction,sep="+")
    	   cov_indexA=which(names(geno_cov_imag_data)==my_covnames[z])
	   cov_indexB=which(names(geno_cov_imag_data)==my_covnames[y])
	   
    	   two_cov<-lm(geno_cov_imag_data$value ~ geno_cov_imag_data[,geno_index] + geno_cov_imag_data[,cov_indexA] + geno_cov_imag_data[,cov_indexB] + geno_cov_imag_data[,geno_index]*geno_cov_imag_data[,cov_indexA])
    
            MElog_pvals=-log10(anova(two_cov)[1,5])
    	    MEeffect_sizes=etasq(two_cov,anova=FALSE,partial=TRUE)[1,]
	    INTlog_pvals=-log10(anova(two_cov)[4,5])
    	    INTeffect_sizes=etasq(two_cov,anova=FALSE,partial=TRUE)[4,]
    	    modelBIC=BIC(two_cov)
    	    numsubjincl=dim(two_cov$model)[1]

    	    addrow <- cbind(SNP,model,imag_variables[i,],MElog_pvals,MEeffect_sizes,INTlog_pvals,INTeffect_sizes,modelBIC,numsubjincl)	  
    	    resultsALL = rbind(resultsALL,addrow)

	   interaction=paste("genotype",my_covnames[y],sep="*")
    	   model=paste("genotype",my_covnames[z],my_covnames[y],interaction,sep="+")
    	   two_cov<-lm(geno_cov_imag_data$value ~ geno_cov_imag_data[,geno_index] + geno_cov_imag_data[,cov_indexA] + geno_cov_imag_data[,cov_indexB] + geno_cov_imag_data[,geno_index]*geno_cov_imag_data[,cov_indexB])
    
            MElog_pvals=-log10(anova(two_cov)[1,5])
    	    MEeffect_sizes=etasq(two_cov,anova=FALSE,partial=TRUE)[1,]
	    INTlog_pvals=-log10(anova(two_cov)[4,5])
    	    INTeffect_sizes=etasq(two_cov,anova=FALSE,partial=TRUE)[4,]
    	    modelBIC=BIC(two_cov)
    	    numsubjincl=dim(two_cov$model)[1]

    	    addrow <- cbind(SNP,model,imag_variables[i,],MElog_pvals,MEeffect_sizes,INTlog_pvals,INTeffect_sizes,modelBIC,numsubjincl)	  
    	    resultsALL = rbind(resultsALL,addrow)
	}
     }
}

# Add triplets of covariates, still estimating genotype effects only
for (z in 1:length(my_covnames)) {
    for (y in 1:length(my_covnames)) {
    	for (x in 1:length(my_covnames)) {
    	    if ((x<y) & (y<z)) {
	       interaction=paste("genotype",my_covnames[z],sep="*")
    	       model=paste("genotype",my_covnames[z],my_covnames[y],my_covnames[x],interaction,sep="+")
    	       cov_indexA=which(names(geno_cov_imag_data)==my_covnames[z])
	       cov_indexB=which(names(geno_cov_imag_data)==my_covnames[y])
    	       cov_indexC=which(names(geno_cov_imag_data)==my_covnames[x])
    	       three_cov<-lm(geno_cov_imag_data$value ~ geno_cov_imag_data[,geno_index] + geno_cov_imag_data[,cov_indexA] + geno_cov_imag_data[,cov_indexB] + geno_cov_imag_data[,cov_indexC] + geno_cov_imag_data[,geno_index]*geno_cov_imag_data[,cov_indexA])
    
		MElog_pvals=-log10(anova(three_cov)[1,5])
    	    	MEeffect_sizes=etasq(three_cov,anova=FALSE,partial=TRUE)[1,]
	    	INTlog_pvals=-log10(anova(three_cov)[5,5])
    	    	INTeffect_sizes=etasq(three_cov,anova=FALSE,partial=TRUE)[5,]
    	    	modelBIC=BIC(three_cov)
    	    	numsubjincl=dim(three_cov$model)[1]

    	    	addrow <- cbind(SNP,model,imag_variables[i,],MElog_pvals,MEeffect_sizes,INTlog_pvals,INTeffect_sizes,modelBIC,numsubjincl)
    	    	resultsALL = rbind(resultsALL,addrow)

	       interaction=paste("genotype",my_covnames[y],sep="*")
    	       model=paste("genotype",my_covnames[z],my_covnames[y],my_covnames[x],interaction,sep="+")
    	       three_cov<-lm(geno_cov_imag_data$value ~ geno_cov_imag_data[,geno_index] + geno_cov_imag_data[,cov_indexA] + geno_cov_imag_data[,cov_indexB] + geno_cov_imag_data[,cov_indexC] + geno_cov_imag_data[,geno_index]*geno_cov_imag_data[,cov_indexB])
    
		MElog_pvals=-log10(anova(three_cov)[1,5])
    	    	MEeffect_sizes=etasq(three_cov,anova=FALSE,partial=TRUE)[1,]
	    	INTlog_pvals=-log10(anova(three_cov)[5,5])
    	    	INTeffect_sizes=etasq(three_cov,anova=FALSE,partial=TRUE)[5,]
    	    	modelBIC=BIC(three_cov)
    	    	numsubjincl=dim(three_cov$model)[1]

    	    	addrow <- cbind(SNP,model,imag_variables[i,],MElog_pvals,MEeffect_sizes,INTlog_pvals,INTeffect_sizes,modelBIC,numsubjincl)
    	    	resultsALL = rbind(resultsALL,addrow)

	       interaction=paste("genotype",my_covnames[x],sep="*")
    	       model=paste("genotype",my_covnames[z],my_covnames[y],my_covnames[x],interaction,sep="+")
    	       three_cov<-lm(geno_cov_imag_data$value ~ geno_cov_imag_data[,geno_index] + geno_cov_imag_data[,cov_indexA] + geno_cov_imag_data[,cov_indexB] + geno_cov_imag_data[,cov_indexC] + geno_cov_imag_data[,geno_index]*geno_cov_imag_data[,cov_indexC])
    
		MElog_pvals=-log10(anova(three_cov)[1,5])
    	    	MEeffect_sizes=etasq(three_cov,anova=FALSE,partial=TRUE)[1,]
	    	INTlog_pvals=-log10(anova(three_cov)[5,5])
    	    	INTeffect_sizes=etasq(three_cov,anova=FALSE,partial=TRUE)[5,]
    	    	modelBIC=BIC(three_cov)
    	    	numsubjincl=dim(three_cov$model)[1]

    	    	addrow <- cbind(SNP,model,imag_variables[i,],MElog_pvals,MEeffect_sizes,INTlog_pvals,INTeffect_sizes,modelBIC,numsubjincl)
    	    	resultsALL = rbind(resultsALL,addrow)



	    }
	}
     }
}

return(resultsALL) 	   
}