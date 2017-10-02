########################
#	Analysis
########################


# Library and Functions
suppressMessages(library(ggplot2))
#library(plotrix)

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


########################
#	Functions
########################

# Input varibles
args <- commandArgs(TRUE)
prefix <- as.character(args[1])
offset <- as.numeric(args[2])
annotated <- as.character(args[3])
extended <- as.character(args[4])
truncated <- as.character(args[5])
novel <- as.character(args[6])


# read files
annotated <- read.table(file=annotated, sep="\t", h=T)
extended <- read.table(file=extended, sep="\t", h=T)
truncated <- read.table(file=truncated, sep="\t", h=T)
novel <- read.table(file=novel, sep="\t", h=T)

offset <- 25

# Annotation
anno.0 <- which(annotated[,1] == 0)
start.anno <- anno.0 - offset
if (start.anno < 0) {start.anno <- 0}
stop.anno <- anno.0 + offset
row.anno <- c(start.anno:stop.anno)
annotated <- annotated[row.anno,]
no.orf <- strsplit(colnames(annotated)[1], "A")[[1]]
anno.no <- no.orf[2]

# Extensions
ext.0 <- which(extended[,1] == 0)
start.ext <- ext.0 - offset
if (start.ext < 0) {start.ext <- 0}
stop.ext <- ext.0 + offset
row.ext <- c(start.ext:stop.ext)
extended <- extended[row.ext,]
no.orf <- strsplit(colnames(extended)[1], "A")[[1]]
ext.no <- no.orf[2]


# Truncations
trunc.0 <- which(truncated[,1] == 0)
start.trunc <- trunc.0 - offset
if (start.trunc < 0) {start.trunc <- 0}
stop.trunc <- trunc.0 + offset
row.trunc <- c(start.trunc:stop.trunc)
truncated <- truncated[row.trunc,]
no.orf <- strsplit(colnames(truncated)[1], "A")[[1]]
trunc.no <- no.orf[2]

# Novel
novel.0 <- which(novel[,1] == 0)
start.novel <- novel.0 - offset
if (start.novel < 0) {start.novel <- 0}
stop.novel <- novel.0 + offset
row.novel <- c(start.novel:stop.novel)
novel <- novel[row.novel,]
no.orf <- strsplit(colnames(novel)[1], "A")[[1]]
novel.no <- no.orf[2]

y.max <- max(c(annotated$average_score, extended$average_score, truncated$average_score, novel$average_score))
y.min <- min(c(annotated$average_score, extended$average_score, truncated$average_score, novel$average_score))

ylim <- c(y.min, y.max)

# Plot
lables <- c(paste(paste("Annotated (",anno.no,sep=""),")", sep=""),paste(paste("Extension (",ext.no,sep=""),")", sep=""),paste(paste("Novel (",novel.no,sep=""),")", sep=""),paste(paste("Truncation (",trunc.no,sep=""),")", sep=""))
colors <- c(1,2,4,3)


cex <- 2.75
pdf(file=paste(prefix,"_conservation.pdf",sep=""), width=14, height=9)
mar.default = c(5, 4, 4, 2) + 0.1
par(mar=mar.default + c(0,4,0,0)) 
plot(annotated, type="l", ylab="Average Conservation score", col=colors[1], ann = FALSE, ylim=ylim, lwd=1.5, cex.axis=cex, cex.lab=3, xaxt="n")
axis(1, at = seq(-offset, offset, by = 5), las=0, cex.axis=cex)
title(ylab = "Average Conservation score", xlab="Distance from start codon", cex.lab =3.5,line = 4)
lines(extended, col=colors[2], lwd=1.5)
lines(truncated, col=colors[4], lwd=1.5)
lines(novel, col=colors[3], lwd=1.5)
add_legend("topright", legend=lables, pch=20,col=colors,horiz=TRUE, bty='n', cex=1)
dev.off()


