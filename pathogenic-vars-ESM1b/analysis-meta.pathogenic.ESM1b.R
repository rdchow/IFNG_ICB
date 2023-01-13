#!/usr/bin/Rscript
library(metafor)
library(MASS)
library(stringr)
library(logistf)

data=read.table("FINAL-DATASET.pretreat.pathogenic.ESM1b.v2.txt",sep="\t",header=TRUE)

# remove patients with non-pathogenic other alterations
data = data[!(data$PredPathogenicAlt_ESM1b == 0 & data$AnyAlt == 1),] #removes 31 patients

cohorts = unique(data$CohortName)
data$BooleanResponse = data$BinaryResponse
data[data$BooleanResponse == "NR","BooleanResponse"] = 0
data[data$BooleanResponse == "R","BooleanResponse"] = 1
data$BooleanResponse = as.numeric(data$BooleanResponse)
data$AnyAlt = factor(data$AnyAlt,levels=c("0","1"))
data$PredPathogenicAlt_ESM1b = factor(data$PredPathogenicAlt_ESM1b,levels=c("0","1"))

# for solo datasets, run multivariable logistic regression with TMB as covariate
# for composite datasets, run simple logistic regression without TMB covariate (since variation in how they are calculated)

# cohort | tissue | Alteration_Estimate | Alteration_Error | Alteration_95CI | Alteration_Pvalue | TMB_Estimate | TMB_Error | TMB_95CI | TMB_Pvalue 

mydf = matrix(nrow=length(cohorts),ncol=12)
for (i in 1:length(cohorts)){
    #c = "Composite_Bladder"
    #c = "Van Allen_Melanoma"
    #c = "Gao_Melanoma"
    c = cohorts[i]
    mydat = data[data$CohortName == c,]
    mytissue = strsplit(c,"_")[[1]][2]
    mydf[i,1] = c
    mydf[i,2] = mytissue

    # if solo
    if (all(mydat$Type == "solo")){ 
        if (sum(is.na(mydat$TMB)) == 0){ # if TMB available for all tumors 

            # binarize TMB
            mydat$TMBbin = ifelse(mydat$TMB >= median(mydat$TMB), 1, 0)
            mydat$TMBbin = factor(mydat$TMBbin,levels=c("0","1"))

            mylm = logistf(BooleanResponse ~ PredPathogenicAlt_ESM1b + TMBbin,data=mydat)
            mycoefs = as.data.frame(t(coef(mylm)))
            mydf[i,3] = mycoefs$PredPathogenicAlt_ESM1b1
            mydf[i,8] = mycoefs$TMBbin1

            myerrors= as.data.frame(t(sqrt(diag(vcov(mylm)))))
            mydf[i,4] = myerrors$PredPathogenicAlt_ESM1b1
            mydf[i,9] = myerrors$TMBbin1

            myconfs = as.data.frame(t(confint(mylm)))
            mydf[i,5] = round(myconfs$PredPathogenicAlt_ESM1b1[1],digits=3)
            mydf[i,6] = round(myconfs$PredPathogenicAlt_ESM1b1[2],digits=3)
            mydf[i,10] = round(myconfs$TMBbin1[1],digits=3)
            mydf[i,11] = round(myconfs$TMBbin1[2],digits=3)

            mypvals = as.data.frame(t(mylm$prob))
            mydf[i,7] = mypvals$PredPathogenicAlt_ESM1b1
            mydf[i,12] = mypvals$TMBbin1

        } else { # if missing TMB data (only Gao cohort)
            mylm = logistf(BooleanResponse ~ PredPathogenicAlt_ESM1b,data=mydat)
            mycoefs = as.data.frame(t(coef(mylm)))
            mydf[i,3] = mycoefs$PredPathogenicAlt_ESM1b1
            mydf[i,7] = NA

            myerrors= as.data.frame(t(sqrt(diag(vcov(mylm)))))
            mydf[i,4] = myerrors$PredPathogenicAlt_ESM1b1
            mydf[i,8] = NA

            myconfs = as.data.frame(t(confint(mylm)))
            mydf[i,5] = round(myconfs$PredPathogenicAlt_ESM1b1[1],digits=3)
            mydf[i,6] = round(myconfs$PredPathogenicAlt_ESM1b1[2],digits=3)
            mydf[i,10] = NA
            mydf[i,11] = NA

            mypvals = as.data.frame(t(mylm$prob))
            mydf[i,7] = mypvals$PredPathogenicAlt_ESM1b1
            mydf[i,12] = NA

        }
    } else if (all(mydat$Type == "composite")){ # if composite
        mylm = logistf(BooleanResponse ~ PredPathogenicAlt_ESM1b,data=mydat)
        mycoefs = as.data.frame(t(coef(mylm)))
        mydf[i,3] = mycoefs$PredPathogenicAlt_ESM1b1
        mydf[i,7] = NA

        myerrors= as.data.frame(t(sqrt(diag(vcov(mylm)))))
        mydf[i,4] = myerrors$PredPathogenicAlt_ESM1b1
        mydf[i,8] = NA

        myconfs = as.data.frame(t(confint(mylm)))
        mydf[i,5] = round(myconfs$PredPathogenicAlt_ESM1b1[1],digits=3)
        mydf[i,6] = round(myconfs$PredPathogenicAlt_ESM1b1[2],digits=3)
        mydf[i,10] = NA
        mydf[i,11] = NA

        mypvals = as.data.frame(t(mylm$prob))
        mydf[i,7] = mypvals$PredPathogenicAlt_ESM1b1
        mydf[i,12] = NA

    }
}

mydf = as.data.frame(mydf)
colnames(mydf) = c("Cohort","Tissue","Alteration_Estimate","Alteration_StdErr","Alteration_95CI.low","Alteration_95CI.hi","Alteration_Pval","TMB_Estimate","TMB_StdErr","TMB_95CI.low","TMB_95CI.hi","TMB_Pval")

mydf$Alteration_Estimate = as.numeric(mydf$Alteration_Estimate)
mydf$Alteration_StdErr = as.numeric(mydf$Alteration_StdErr)
mydf$Alteration_Variance = mydf$Alteration_StdErr^2
mydf = mydf[with(mydf,order(Tissue,Cohort)),]
write.table(mydf,"regression-results.pretreat.pathogenic.removeNonpathMutated.ESM1b.v2.txt",sep="\t",row.names=FALSE)


# run meta-analysis

###############################################
# multi-level regression
mydf = read.table("regression-results.pretreat.pathogenic.removeNonpathMutated.ESM1b.v2.txt",sep="\t",header=TRUE)
so = read.table("sample-order.ESM1b.txt",sep="\t",header=TRUE)
so = so[match(mydf$Cohort,so$Cohort),]
mydf$Order = so$Order
mydf = mydf[mydf$Tissue != "Breast",]

regres = rma.mv(yi=Alteration_Estimate,V=Alteration_Variance,slab=Cohort,data=mydf,method="REML",random= ~1|Tissue/Cohort)
myres=rma(yi=Alteration_Estimate,vi=Alteration_Variance,data=mydf)

W <- diag(1/regres$vi)
X <- model.matrix(regres)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(regres$sigma2) / (sum(regres$sigma2) + (regres$k-regres$p)/sum(diag(P)))

# lung subgroup
lungdf = mydf[mydf$Tissue == "Lung",]
lungres = rma(yi=Alteration_Estimate,vi=Alteration_Variance,data=lungdf)

# melanoma subgroup
meldf = mydf[mydf$Tissue == "Melanoma",]
melres = rma(yi=Alteration_Estimate,vi=Alteration_Variance,data=meldf)

##############################################
# overall forest
mycols = mydf$Tissue
mycols[mycols=='Bladder'] = "#000C95"
mycols[mycols=='Lung'] = "#7876B1FF"
mycols[mycols=='Melanoma'] = "#20854EFF"
mycols[mycols=='Breast'] = "#CC347C"
mycols[mycols=='CRC'] = "#2FACCF"
mycols[mycols=='Endometrial'] = "#E18727FF"
mycols[mycols=='Gastric'] = "#6D3000"
mycols[mycols=='RCC'] = "#999933"

pdf("forest-plot-multilevelReg.pretreat.pathogenic.removeNonpathMutated.ESM1b.v2.pdf",height=7.5,width=8,useDingbats=FALSE)
forest(mydf$Alteration_Estimate, ci.lb=mydf$Alteration_95CI.low, ci.ub=mydf$Alteration_95CI.hi, slab=mydf$Cohort, header=TRUE, ilab=mydf[,c("Tissue")],ilab.xpos=c(-7), order=mydf$Order,col=mycols,alim=c(-4,4),steps=5,plim=c(1,2),ylim=c(-1,22))

addpoly(regres, row=-0.8, mlab="Mixed Effects",col="black",addpred=TRUE,cex=1.25)
addpoly(melres, row=0.5, mlab="Melanoma subgroup",col="#20854EFF",addpred=TRUE,cex=1.25)
addpoly(lungres, row= 8.5, mlab="Lung subgroup",col="#7876B1FF",addpred=TRUE,cex=1.25)
dev.off()