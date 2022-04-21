# /usr/local/pkg/R/3.4.2/bin/R
# qsub -cwd -l lmem,s_vmem=32G,mem_req=32G -N ewas -b y '/usr/local/pkg/R/3.4.2/bin/Rscript --vanilla script.R'

# load packages
library(data.table)
library(ggplot2)

# load beta
setwd("path/to/workind/directory")
beta <- as.data.frame(fread("path/to/DNAm/data(bed format)"),stringsAsFactors=F)

# load covariates
baseformula <- "m ~ group + age + sex + CD8T + CD4T + NK + Bcell + Mono + Neu"
cov <- read.table("path/to/phenotype/data",header=T,stringsAsFactors=F,row.names=1)
cell <- read.csv("path/to/cellCount/data",header=T,row.names=1)

# get ID
beta_id <- colnames(beta)[5:ncol(beta)]
cov_id <- row.names(cov)
cell_id <- row.names(cell)
id <- intersect(beta_id,intersect(cov_id,cell_id))

# reshape
beta2 <- beta[,c("chr","end","cr",id)]
cov2 <- cov[id,]
cell2 <- cell[id,]

# get covariates
CD8T  <- cell2[,"CD8T"]
CD4T  <- cell2[,"CD4T"]
NK    <- cell2[,"NK"]
Bcell <- cell2[,"Bcell"]
Mono  <- cell2[,"Mono"]
Neu  <- cell2[,"Neu"]
group <- as.factor(cov2[,"cancer"])
age   <- as.numeric(cov2[,"age"])
sex   <- as.factor(cov2[,"sex"])

# define function
lmfunc <- function(x){
	# browser()
	counter + 1 ->> counter
	#cat(counter)
	chr <- as.character(x[1])
	pos <- as.character(x[2])
	cr <- as.numeric(x[3])
	m <- as.numeric(x[4:length(x)])
	model <- as.formula(baseformula)
	res <- c(chr,pos)
	ri <- (as.numeric(as.character(quantile(na.exclude(m),0.95)))) - (as.numeric(as.character(quantile(na.exclude(m),0.05))))
	res <- c(res,cr,ri)
	if(cr < 0.5) {res <- c(res, as.numeric(rep(NA,4)))
	} else {
		res.lm <- NA
		tryCatch({res.lm <- lm(model, na.action=na.omit)},error = function(e) {res.lm <<- NA})		
		if(exists("res.lm")){
			if(!is.na(res.lm)){
				res0 <- summary(res.lm)
				res1 <- res0$coefficients
				coe_group    <- res1[2,1]
				se_group     <- res1[2,2]
				Pval_group   <- res1[2,4]
				rsq <- summary(res.lm)$r.squared
				res <- c(res,coe_group,se_group,Pval_group,rsq)
			} else {res <- c(res,as.numeric(rep(NA,4)))
			}
		} else {res <- c(res,as.numeric(rep(NA,4)))
		}
	}
}

counter <- 0
# lm
lmn <- 0
for(i in 1:nrow(beta2)){
	lmn <- lmn + 1
	startline <- lmn
	endline <- lmn + 1000
	if(endline >= nrow(beta2)){
		endline <- nrow(beta2)
	}
	beta2split <- beta2[c(startline:endline),]
	lm.res.split <- apply(beta2split,1,lmfunc)
	lm.res.split <- as.data.frame(t(lm.res.split))
	colnames(lm.res.split) <- c("chr", "pos",
		"cr","ri","coef","se","pval","rsq"
	)
	outfn <- "lm.res"
	if(lmn == 1){
		write.table(lm.res.split,outfn,row.names=F,quote=F,sep="\t")
	} else {
		write.table(lm.res.split,outfn,row.names=F,col.names=F,quote=F,sep="\t",append=T)
	}
	lmn <- endline
	if(endline == nrow(beta2)){
		break
	}
}

