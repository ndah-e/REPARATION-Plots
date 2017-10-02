############
#	R script for random forest ORF delineation
############

# Library
suppressMessages(library(ggplot2))
suppressMessages(library(randomForest))
suppressMessages(library(DMwR))

run_randomforest <- function(train){
	set.seed(12345)	# set seed
	rf.output <- randomForest(class~., data=train)
	return(rf.output)
}

run_randomforest_2 <- function(train,pos,neg){
	set.seed(12345)	# set seed
	rf.output <- randomForest(class~., data=train,classwt = c(P=pos, N=neg))
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

over_sampling_ecoli <- vector()
under_sampling_ecoli <- vector()
smote_ecoli <- vector()
classwt_ecoli <- vector()

over_sampling_salt_mono <- vector()
under_sampling_salt_mono <- vector()
smote_salt_mono <- vector()
classwt_salt_mono <- vector()



for(i in 1:nfolds) {

	print(paste("processing fold ", i))

# ----------- ecoli ----------- #
	testIndexes.ecoli.p <- which(folds.ecoli.p==i,arr.ind=TRUE)
	testIndexes.ecoli.n <- which(folds.ecoli.n==i,arr.ind=TRUE)

	testset.ecoli.p <- positive.ecoli[testIndexes.ecoli.p, feat]
	trainset.ecoli.p <- positive.ecoli[-testIndexes.ecoli.p, feat]
	testset.ecoli.n <- negative.ecoli[testIndexes.ecoli.n, feat]
	trainset.ecoli.n <- negative.ecoli[-testIndexes.ecoli.n, feat]

	testset.ecoli <- rbind(testset.ecoli.p,testset.ecoli.n)
	trainset.ecoli <- rbind(trainset.ecoli.p,trainset.ecoli.n)

	# UNDER SAMPLING
	set.seed(12345)
	under_samp.ecoli.n <- sample(1:nrow(trainset.ecoli.n),nrow(trainset.ecoli.p), replace = FALSE)
	trainset.ecoli.under.samp <- trainset.ecoli.n[under_samp.ecoli.n,]
	trainset.ecoli.under.samp <- rbind(trainset.ecoli.under.samp,trainset.ecoli.p)
	
	rf.ecoli <- run_randomforest(trainset.ecoli.under.samp)
	testset.ecoli$pred <- predict(rf.ecoli, testset.ecoli, type="response")
	ecoli_conf.matrix <- table(testset.ecoli$pred, testset.ecoli$class)
	ecoli_TP <- ecoli_conf.matrix[1,1]
	ecoli_FP <- ecoli_conf.matrix[1,2]
	ecoli_FN <- ecoli_conf.matrix[2,1]
    under_sampling_ecoli[i] <- ecoli_TP/(ecoli_TP + ecoli_FP)

	# OVER SAMPLING
	set.seed(12345)
	over_samp.ecoli.p <- sample(1:nrow(trainset.ecoli.p),nrow(trainset.ecoli.n), replace = TRUE)
	trainset.ecoli.over.samp <- trainset.ecoli.p[over_samp.ecoli.p,]
	trainset.ecoli.over.samp <- rbind(trainset.ecoli.over.samp,trainset.ecoli.n)
	
	rf.ecoli <- run_randomforest(trainset.ecoli.over.samp)
	testset.ecoli$pred <- predict(rf.ecoli, testset.ecoli, type="response")
	ecoli_conf.matrix <- table(testset.ecoli$pred, testset.ecoli$class)
	ecoli_TP <- ecoli_conf.matrix[1,1]
	ecoli_FP <- ecoli_conf.matrix[1,2]
	ecoli_FN <- ecoli_conf.matrix[2,1]
    over_sampling_ecoli[i] <- ecoli_TP/(ecoli_TP + ecoli_FP)

	# SMOTE
	trainset.ecoli.SMOTE <- SMOTE(class ~ ., trainset.ecoli, perc.over = 200, perc.under=100)
	rf.ecoli <- run_randomforest(trainset.ecoli.SMOTE)
	testset.ecoli$pred <- predict(rf.ecoli, testset.ecoli, type="response")
	ecoli_conf.matrix <- table(testset.ecoli$pred, testset.ecoli$class)
	ecoli_TP <- ecoli_conf.matrix[1,1]
	ecoli_FP <- ecoli_conf.matrix[1,2]
	ecoli_FN <- ecoli_conf.matrix[2,1]
    smote_ecoli[i] <- ecoli_TP/(ecoli_TP + ecoli_FP)

	#CLASSWT_ecoli
	pos.ecoli <- nrow(trainset.ecoli.p)/nrow(trainset.ecoli.n)
	neg.ecoli <- 1 - pos.ecoli

	rf.ecoli <- run_randomforest_2(trainset.ecoli,pos.ecoli,neg.ecoli)
	testset.ecoli$pred <- predict(rf.ecoli, testset.ecoli, type="response")
	ecoli_conf.matrix <- table(testset.ecoli$pred, testset.ecoli$class)
	ecoli_TP <- ecoli_conf.matrix[1,1]
	ecoli_FP <- ecoli_conf.matrix[1,2]
	ecoli_FN <- ecoli_conf.matrix[2,1]
    classwt_ecoli[i] <- ecoli_TP/(ecoli_TP + ecoli_FP)

# ----------- Salt ----------- #
	testIndexes.salt_mono.p <- which(folds.salt_mono.p==i,arr.ind=TRUE)
	testIndexes.salt_mono.n <- which(folds.salt_mono.n==i,arr.ind=TRUE)

	testset.salt_mono.p <- positive.salt_mono[testIndexes.salt_mono.p, feat]
	trainset.salt_mono.p <- positive.salt_mono[-testIndexes.salt_mono.p, feat]
	testset.salt_mono.n <- negative.salt_mono[testIndexes.salt_mono.n, feat]
	trainset.salt_mono.n <- negative.salt_mono[-testIndexes.salt_mono.n, feat]

	trainset.salt_mono <- rbind(trainset.salt_mono.p,trainset.salt_mono.n)
	testset.salt_mono <- rbind(testset.salt_mono.p,testset.salt_mono.n)

	# UNDER SAMPLING
	set.seed(12345)
	under_samp.salt_mono.n <- sample(1:nrow(trainset.salt_mono.n),nrow(trainset.salt_mono.p), replace = FALSE)
	trainset.salt_mono.under.samp <- trainset.salt_mono.n[under_samp.salt_mono.n,]
	trainset.salt_mono.under.samp <- rbind(trainset.salt_mono.under.samp,trainset.salt_mono.p)
	
	rf.salt_mono <- run_randomforest(trainset.salt_mono.under.samp)
	testset.salt_mono$pred <- predict(rf.salt_mono, testset.salt_mono, type="response")
	salt_mono_conf.matrix <- table(testset.salt_mono$pred, testset.salt_mono$class)
	salt_mono_TP <- salt_mono_conf.matrix[1,1]
	salt_mono_FP <- salt_mono_conf.matrix[1,2]
	salt_mono_FN <- salt_mono_conf.matrix[2,1]
    under_sampling_salt_mono[i] <- salt_mono_TP/(salt_mono_TP + salt_mono_FP)

	# OVER SAMPLING
	set.seed(12345)
	over_samp.salt_mono.p <- sample(1:nrow(trainset.salt_mono.p),nrow(trainset.salt_mono.n), replace = TRUE)
	trainset.salt_mono.over.samp <- trainset.salt_mono.p[over_samp.salt_mono.p,]
	trainset.salt_mono.over.samp <- rbind(trainset.salt_mono.over.samp,trainset.salt_mono.n)
	
	rf.salt_mono <- run_randomforest(trainset.salt_mono.over.samp)
	testset.salt_mono$pred <- predict(rf.salt_mono, testset.salt_mono, type="response")
	salt_mono_conf.matrix <- table(testset.salt_mono$pred, testset.salt_mono$class)
	salt_mono_TP <- salt_mono_conf.matrix[1,1]
	salt_mono_FP <- salt_mono_conf.matrix[1,2]
	salt_mono_FN <- salt_mono_conf.matrix[2,1]
    over_sampling_salt_mono[i] <- salt_mono_TP/(salt_mono_TP + salt_mono_FP)

	# SMOTE
	trainset.salt_mono.SMOTE <- SMOTE(class ~ ., trainset.salt_mono, perc.over = 200, perc.under=100)
	rf.salt_mono <- run_randomforest(trainset.salt_mono.SMOTE)
	testset.salt_mono$pred <- predict(rf.salt_mono, testset.salt_mono, type="response")
	salt_mono_conf.matrix <- table(testset.salt_mono$pred, testset.salt_mono$class)
	salt_mono_TP <- salt_mono_conf.matrix[1,1]
	salt_mono_FP <- salt_mono_conf.matrix[1,2]
	salt_mono_FN <- salt_mono_conf.matrix[2,1]
    smote_salt_mono[i] <- salt_mono_TP/(salt_mono_TP + salt_mono_FP)

	#CLASSWT_salt_mono
	pos.salt_mono <- nrow(trainset.salt_mono.p)/nrow(trainset.salt_mono.n)
	neg.salt_mono <- 1 - pos.salt_mono

	rf.salt_mono <- run_randomforest_2(trainset.salt_mono,pos.salt_mono,neg.salt_mono)
	testset.salt_mono$pred <- predict(rf.salt_mono, testset.salt_mono, type="response")
	salt_mono_conf.matrix <- table(testset.salt_mono$pred, testset.salt_mono$class)
	salt_mono_TP <- salt_mono_conf.matrix[1,1]
	salt_mono_FP <- salt_mono_conf.matrix[1,2]
	salt_mono_FN <- salt_mono_conf.matrix[2,1]
    classwt_salt_mono[i] <- salt_mono_TP/(salt_mono_TP + salt_mono_FP)

}

head <- c("over_sampling_ecoli","under_sampling_ecoli","smote_ecoli","classwt_ecoli","over_sampling_salt_mono","under_sampling_salt_mono","smote_salt_mono","classwt_salt_mono")
precision <- c(mean(over_sampling_ecoli),mean(under_sampling_ecoli),mean(smote_ecoli),mean(classwt_ecoli),mean(over_sampling_salt_mono),mean(under_sampling_salt_mono),mean(smote_salt_mono),mean(classwt_salt_mono))

df <- data.frame(type=head,precision=precision)
write.table(df, file="imbalance_precision.txt", sep = "\t",col.names = TRUE,row.names = F, quote = FALSE)

pdf(file="imbalance_precision.pdf", width=9)
#barplot (df$precision, las=1, col=rainbow(nrow(df)), names=df$type, main="CV precision For Imbalance Analysis")
plot (df, las=1, col=rainbow(nrow(df)), main="CV precision For Imbalance Analysis")
dev.off()



