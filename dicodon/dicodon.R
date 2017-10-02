########################
#	Analysis
########################


# Library and Functions
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
library(RColorBrewer)
require(lattice)
library(plotrix)

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
file <- as.character(args[1])
prefix <- as.character(args[2])

# data
dicodon <- read.table(file=file, sep="\t", h=T)


# PLOT ORF categories
codon <- paste(prefix,"_dicodon_bar.pdf")
pdf(file=codon, width=10)
#barplot(as.matrix(dicodon[,c(2:4)]), col= c("blue", "red", "green"), legend=c("Annotated","Novel","Random"), cex.lab = 2.5, cex.axis = 3,main="", ylab="Frequency", xlab="Amino acids", names.arg=dicodon$aa, beside=T)
#text(labels=dicodon$aa, cex= 2.5, col= "blue")

trellis.par.set(fontsize=list(text=30))
dt <- data.frame(Category=c(rep("Annotated", 20),rep("Novel",20),rep("Random",20)), 
                aa=c(as.vector(dicodon$aa),as.vector(dicodon$aa),as.vector(dicodon$aa)),Frequency=c(dicodon$Annotated, dicodon$Novel, dicodon$Random))
#barchart(Frequency ~ aa, groups=Category, data=dt,scales=list(x=list(cex=1)), col=c("darkblue","red", "green"),auto.key=TRUE,par.settings = simpleTheme(col=c("darkblue","red", "green")))
barchart(Frequency ~ aa, groups=Category, data=dt,scales=list(x=list(cex=1)), col=c("darkblue","red", "green"))

#axis.default(labels = as.vector(dicodon$aa))

dev.off()


# PLOT ORF categories
#codon <- paste(prefix,"_dicodon.pdf")
#pdf(file=codon, width=10)
#plot(dicodon$Annotated, dicodon$Random, col= "blue", pch = 19, cex=0.01, cex.lab = 2.5, cex.axis = 3,main=paste(prefix, cor(dicodon$Annotated, dicodon$Random), sep=" correlation="), ylab="Random Sequence", xlab="Annotated ORFs")
#text(dicodon$Annotated, dicodon$Random, labels=dicodon$aa, cex= 2.5, col= "blue")
#plot(dicodon$Novel, dicodon$Random, col= "blue", pch = 19, cex=0.01, cex.lab = 2.5,cex.axis = 3, main=paste(prefix, cor(dicodon$Novel, dicodon$Random), sep=" correlation="), ylab="Random sequences", xlab="Novel ORFs")
#text(dicodon$Novel, dicodon$Random, labels=dicodon$aa, cex= 2.5, col= "blue")
#plot(dicodon$Annotated, dicodon$Novel, col= "blue", pch = 19, cex=0.01, cex.lab=2.5, cex.axis=3,main=paste(prefix, cor(dicodon$Novel, dicodon$Annotated), sep=" correlation="), ylab="Novel ORFs", xlab="Annotated ORFs")
#text(dicodon$Annotated, dicodon$Novel, labels=dicodon$aa, cex= 2.5, col= "blue")
#dev.off()



