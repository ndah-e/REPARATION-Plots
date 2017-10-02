############
#	R script for random forest ORF delineation
############

# Library
suppressMessages(library(ggplot2))
suppressMessages(library(randomForest))
suppressMessages(library(SiZer))


# Functions
run_randomforest <- function(train,pos,neg,mtry){
	set.seed(131)	# set seed
	rf.output <- randomForest(class~., data=train, mtry=4, ntree=3001, nodesize=4,classwt = c(P=pos, N=neg))
	return(rf.output)
}

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

getROC_AUC = function(probs, true_Y){
    probsSort = sort(probs, decreasing = TRUE, index.return = TRUE)
    val = unlist(probsSort$x)
    idx = unlist(probsSort$ix)  

    roc_y = true_Y[idx];
    stack_x = cumsum(roc_y == 2)/sum(roc_y == 2)
    stack_y = cumsum(roc_y == 1)/sum(roc_y == 1)    

    auc = sum((stack_x[2:length(roc_y)]-stack_x[1:length(roc_y)-1])*stack_y[2:length(roc_y)])
	print(auc)
    return(list(stack_x=stack_x, stack_y=stack_y, auc=auc))
}

feat <- c("class","start_coverage","start_rpkm","coverage","accu_prop","stop_rpkm","SD_score")


# Input varibles
args <- commandArgs(TRUE)
ORFs_file <- as.character(args[1])
Positive_set <- as.character(args[2])




###################################
### Ecoli data
###################################

ORFs <- read.table(file=ORFs_file, sep="\t", h=T)
prodigal <- read.table(file=Positive_set, sep="\t", h=T)

ORFs <- ORFs[which(ORFs$start_rpkm > 0),]
ORFs$log_rpkm <- log(ORFs$rpkm)

positive <- ORFs[which(as.vector(ORFs$orf_id) %in% as.vector(prodigal$orf_id)),]
positive$class <- "P"
positive$class <- as.factor(positive$class)

# S-CURVE Threshold Estimation
model <- nls(coverage ~ SSfpl(log_rpkm, D, A, C, B), data=positive)
x = fitted(model)
y <- predict(model)

# SiZer bent point prediction
bent <- bent.cable(x, y, grid.size=100)
MINRPKM <- round(bent$alpha, digits = 2)
MINCOV <- round(predict(model, newdata = data.frame(log_rpkm = MINRPKM))[1], digits = 2)

#positive set
positive <- positive[which(positive$log_rpkm >= MINRPKM),]
positive <- positive[which(positive$coverage >= MINCOV),]

# keep all ORFs above threshold
ORFs <- ORFs[which(ORFs$log_rpkm >= MINRPKM),]
ORFs <- ORFs[which(ORFs$coverage >= MINCOV),]

negative <- ORFs[which(ORFs$start_codon == 'CTG' & ORFs$length >= min(positive$length)),]
negative$class <- "N"
negative$class <- as.factor(negative$class)

# keep only longest ORF
family.neg <- unique(as.vector(negative$gene))
family.neg.vector <- character(length = length(family.neg))
for (i in 1:length(family.neg)) {
	g <- family.neg[i]
	family <- negative[which(negative$gene == g),]
	idx <- sample(dim(family)[1],1)
	family.neg.vector[i] <- as.vector(family[idx,]$orf_id)
}
negative <- negative[which(as.vector(negative$orf_id) %in% family.neg.vector),]


###################################
### Manual Search tuning
###################################

# ntree
nfolds <- 10
set.seed(131)
folds.p <- sample(1:nfolds, nrow(positive), replace = TRUE)
set.seed(131)
folds.n <- sample(1:nfolds, nrow(negative), replace = TRUE)

accuracy <- vector()
precision <- vector()
sensitivity <- vector()
auc <- vector()

#aList = getROC_AUC(probs, true_Y) 
#auc = unlist(aList$auc)
for(i in 1:nfolds) {

	# ecoli
	testIndexes.p <- which(folds.p==i,arr.ind=TRUE)
	testIndexes.n <- which(folds.n==i,arr.ind=TRUE)

	testset.p <- positive[testIndexes.p, feat]
	trainset.p <- positive[-testIndexes.p, feat]
	testset.n <- negative[testIndexes.n, feat]
	trainset.n <- negative[-testIndexes.n, feat]

	testset <- rbind(testset.p,testset.n)
	trainset <- rbind(trainset.p,trainset.n)

	pos.prob <- dim(trainset.p)[1]/dim(trainset.n)[1]
	neg.prob <- 1 - pos.prob

	rf <- run_randomforest(trainset,pos.prob,neg.prob,mtry[j])
	testset$prob <- predict(rf, testset, type="prob")[,2]
	testset$pred <- predict(rf, testset, type="response")

	conf.matrix <- table(testset$pred, testset$class)
	TP <- conf.matrix[1,1]
	FP <- conf.matrix[1,2]
	FN <- conf.matrix[2,1]
	TN <- conf.matrix[2,1]

	precision[i] <- TP/(TP + FP)
	sensitivity[i] <- TP/(TP + FN)

	aList = getROC_AUC(testset$prob, testset$class)
	print(unlist(aList$auc))
	auc[i] <- unlist(aList$auc)

    testset$rightPred <- testset$pred == testset$class
    accuracy[i] <- sum(testset$rightPred)/nrow(testset)
}

cat("Accuracy ", round(mean(accuracy)*100,digit=2), "\n")
cat("Precision ", round(mean(precision)*100,digit=2), "\n")
cat("Sensitivity ", round(mean(sensitivity)*100,digit=2), "\n")
cat("AUC ", round(mean(auc)*100,digit=2), "\n")





