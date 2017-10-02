############
#	R script for random forest ORF delineation
############

# Library
suppressMessages(library(ggplot2))
suppressMessages(library(randomForest))

run_randomforest <- function(train,pos,neg,ntree){
	set.seed(12345)	# set seed
	rf.output <- randomForest(class~., data=train, mtry=3,ntree=ntree,classwt = c(P=pos, N=neg))
	return(rf.output)
}

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

feat <- c("class","start_coverage","start_rpkm","coverage","accu_prop","stop_rpkm","SD_score")


###################################
### Ecoli data
###################################

ORFs.ecoli <- read.table(file="/data/elvis/REPARATION/Ecoli_sample/Ribo_Rich_define_media_rep_1/work_dir/all_valid_ORFs.txt", sep="\t", h=T)
prodigal.ecoli <- read.table(file="/data/elvis/REPARATION/Ecoli_sample/Ribo_Rich_define_media_rep_1/work_dir/positive_set.txt", sep="\t", h=T)

ORFs.ecoli <- ORFs.ecoli[which(ORFs.ecoli$start_rpkm > 0),]
ORFs.ecoli$log_rpkm <- log(ORFs.ecoli$rpkm)
ORFs.ecoli <- ORFs.ecoli[which(ORFs.ecoli$log_rpkm >= 0),]
ORFs.ecoli <- ORFs.ecoli[which(ORFs.ecoli$coverage >= 0.1),]

positive.ecoli <- ORFs.ecoli[which(as.vector(ORFs.ecoli$orf_id) %in% as.vector(prodigal.ecoli$orf_id)),]
positive.ecoli <- positive.ecoli[which(positive.ecoli$start_codon == "ATG"),]
positive.ecoli$class <- "P"
positive.ecoli$class <- as.factor(positive.ecoli$class)

negative.ecoli <- ORFs.ecoli[which(ORFs.ecoli$start_codon == 'CTG' & ORFs.ecoli$length >= 90),]
negative.ecoli$class <- "N"
negative.ecoli$class <- as.factor(negative.ecoli$class)

# keep only longest ORF
family.neg.ecoli <- unique(as.vector(negative.ecoli$gene))
family.neg.vector.ecoli <- character(length = length(family.neg.ecoli))
for (i in 1:length(family.neg.ecoli)) {
	g <- family.neg.ecoli[i]
	family.ecoli <- negative.ecoli[which(negative.ecoli$gene == g),]
	idx.ecoli <- sample(dim(family.ecoli)[1],1)
	family.neg.vector.ecoli[i] <- as.vector(family.ecoli[idx.ecoli,]$orf_id)
}
negative.ecoli <- negative.ecoli[which(as.vector(negative.ecoli$orf_id) %in% family.neg.vector.ecoli),]


###################################
### SALT_mono data
###################################

ORFs.salt_mono <- read.table(file="/data/elvis/REPARATION/SALT_sample/Mono/work_dir/all_valid_ORFs.txt", sep="\t", h=T)
prodigal.salt_mono <- read.table(file="/data/elvis/REPARATION/SALT_sample/Mono/work_dir/positive_set.txt", sep="\t", h=T)

ORFs.salt_mono <- ORFs.salt_mono[which(ORFs.salt_mono$start_rpkm > 0),]
ORFs.salt_mono$log_rpkm <- log(ORFs.salt_mono$rpkm)
ORFs.salt_mono <- ORFs.salt_mono[which(ORFs.salt_mono$log_rpkm >= 0),]
ORFs.salt_mono <- ORFs.salt_mono[which(ORFs.salt_mono$coverage >= 0.1),]

positive.salt_mono <- ORFs.salt_mono[which(as.vector(ORFs.salt_mono$orf_id) %in% as.vector(prodigal.salt_mono$orf_id)),]
positive.salt_mono$class <- "P"
positive.salt_mono <- positive.salt_mono[which(positive.salt_mono$start_codon == "ATG"),]
positive.salt_mono$class <- as.factor(positive.salt_mono$class)

negative.salt_mono <- ORFs.salt_mono[which(ORFs.salt_mono$start_codon == 'CTG' & ORFs.salt_mono$length >= 90),]
negative.salt_mono$class <- "N"
negative.salt_mono$class <- as.factor(negative.salt_mono$class)

# keep only longest ORF
family.neg.salt_mono <- unique(as.vector(negative.salt_mono$gene))
family.neg.vector.salt_mono <- character(length = length(family.neg.salt_mono))
for (i in 1:length(family.neg.salt_mono)) {
	g <- family.neg.salt_mono[i]
	family.salt_mono <- negative.salt_mono[which(negative.salt_mono$gene == g),]
	idx.salt_mono <- sample(dim(family.salt_mono)[1],1)
	family.neg.vector.salt_mono[i] <- as.vector(family.salt_mono[idx.salt_mono,]$orf_id)
}
negative.salt_mono <- negative.salt_mono[which(as.vector(negative.salt_mono$orf_id) %in% family.neg.vector.salt_mono),]



###################################
### Manual Search tuning
###################################

# ntree
nfolds <- 10
set.seed(12345)
folds.ecoli.p <- sample(1:nfolds, nrow(positive.ecoli), replace = TRUE)
set.seed(12345)
folds.ecoli.n <- sample(1:nfolds, nrow(negative.ecoli), replace = TRUE)
set.seed(12345)
folds.salt_mono.p <- sample(1:nfolds, nrow(positive.salt_mono), replace = TRUE)
set.seed(12345)
folds.salt_mono.n <- sample(1:nfolds, nrow(negative.salt_mono), replace = TRUE)


ecoli <- vector()
salt_mono <- vector()

nodesize <- 6
mtry <- 3
ntree <- c(101,251,seq(501,6001,500))
for (j in 1:length(ntree)) {

	salt_mono_cv_error <- vector()
	ecoli_cv_error <- vector()

	cat("Processing ntree ",ntree[j],"\n", sep=" ")

	for(i in 1:nfolds) {

		# ecoli
	    testIndexes.ecoli.p <- which(folds.ecoli.p==i,arr.ind=TRUE)
	    testIndexes.ecoli.n <- which(folds.ecoli.n==i,arr.ind=TRUE)

	    testset.ecoli.p <- positive.ecoli[testIndexes.ecoli.p, feat]
		trainset.ecoli.p <- positive.ecoli[-testIndexes.ecoli.p, feat]
	    testset.ecoli.n <- negative.ecoli[testIndexes.ecoli.n, feat]
		trainset.ecoli.n <- negative.ecoli[-testIndexes.ecoli.n, feat]

		testset.ecoli <- rbind(testset.ecoli.p,testset.ecoli.n)
		trainset.ecoli <- rbind(trainset.ecoli.p,trainset.ecoli.n)

		pos.prob.ecoli <- dim(trainset.ecoli.p)[1]/dim(trainset.ecoli.n)[1]
		neg.prob.ecoli <- 1 - pos.prob.ecoli

		rf.ecoli <- run_randomforest(trainset.ecoli,pos.prob.ecoli,neg.prob.ecoli,ntree[j])
		testset.ecoli$pred <- predict(rf.ecoli, testset.ecoli, type="response")
		ecoli_conf.matrix <- table(testset.ecoli$pred, testset.ecoli$class)
		ecoli_TP <- ecoli_conf.matrix[1,1]
		ecoli_FP <- ecoli_conf.matrix[1,2]
		ecoli_FN <- ecoli_conf.matrix[2,1]
		ecoli_cv_error[i] <- ecoli_TP/(ecoli_TP + ecoli_FP)

		#Salt
	    testIndexes.salt_mono.p <- which(folds.salt_mono.p==i,arr.ind=TRUE)
	    testIndexes.salt_mono.n <- which(folds.salt_mono.n==i,arr.ind=TRUE)

	    testset.salt_mono.p <- positive.salt_mono[testIndexes.salt_mono.p, feat]
		trainset.salt_mono.p <- positive.salt_mono[-testIndexes.salt_mono.p, feat]
	    testset.salt_mono.n <- negative.salt_mono[testIndexes.salt_mono.n, feat]
		trainset.salt_mono.n <- negative.salt_mono[-testIndexes.salt_mono.n, feat]

		trainset.salt_mono <- rbind(trainset.salt_mono.p,trainset.salt_mono.n)
		testset.salt_mono <- rbind(testset.salt_mono.p,testset.salt_mono.n)

		pos.prob.salt_mono <- dim(trainset.salt_mono.p)[1]/dim(trainset.salt_mono.n)[1]
		neg.prob.salt_mono <- 1 - pos.prob.salt_mono

		rf.salt_mono <- run_randomforest(trainset.salt_mono,pos.prob.salt_mono,neg.prob.salt_mono,ntree[j])
		testset.salt_mono$pred <- predict(rf.salt_mono, testset.salt_mono, type="response")
		salt_mono_conf.matrix <- table(testset.salt_mono$pred, testset.salt_mono$class)
		salt_mono_TP <- salt_mono_conf.matrix[1,1]
		salt_mono_FP <- salt_mono_conf.matrix[1,2]
		salt_mono_FN <- salt_mono_conf.matrix[2,1]
        salt_mono_cv_error[i] <- salt_mono_TP/(salt_mono_TP + salt_mono_FP)
	}

	ecoli[j] <- mean(ecoli_cv_error)
	salt_mono[j] <- mean(salt_mono_cv_error)

}

df <- data.frame(ntree=ntree,ecoli=ecoli,salt_mono=salt_mono)
write.table(df, file="ntree.txt", sep = "\t",col.names = TRUE,row.names = F, quote = FALSE)

min.ecoli <- paste(df[which(df$ecoli == min(df$ecoli)),]$ntree, collapse =" ")
min.salt <- paste(df[which(df$salt == min(df$salt)),]$ntree, collapse =" ")

cat("Ecoli ", min.ecoli, "\n")
cat("Salt ", min.salt, "\n")

pdf(file="ntree_tuning.pdf")
ymax <- max(df[,-1])
ymin <- min(df[,-1])
plot(df[,1],df[,2], xlab="ntree", ylab="Precision", col=2, pch=16, ylim=c(ymin,ymax))
lines(df[,1],df[,2], col=2)
for (i in 3:ncol(df)) {
	points(df[,1],df[,i], pch=16,col=i)
	lines(df[,1],df[,i], pch=16,col=i)
}

cex <- 1.25
colnames <- colnames(df)
add_legend("topright", legend=colnames[2:ncol(df)],pch=20,col=c(2:ncol(df)),horiz=TRUE, bty='n', cex=cex)
dev.off()



