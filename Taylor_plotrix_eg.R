# Taylor diagram example
# From plotrix package

library(datasets)
library(ncdf4)
library(plotrix)

#use 'add=TRUE' after plotting the first dataset

taylor.diagram(as.vector(volcano), # makes a vector
               as.vector(volcano), # makes a vector
               add=FALSE,
               col="red",
               pch=4,
               pos.cor=TRUE,
               xlab="MERRA SD (Normalised)",
               ylab="RCA4 runs SD (normalised)",
               main="Taylor Diagram",
               show.gamma=TRUE,
               ngamma=3,
               sd.arcs=1,
               ref.sd=TRUE,
               grad.corr.lines=c(0.2,0.4,0.6,0.8,0.9),
               pcex=1,cex.axis=1,
               normalize=TRUE,
               mar=c(5,4,6,6),
               lwd=10,
               font=5,
               lty=3)

taylor.diagram(as.vector(volcano),
               as.vector(volcano),
               add=FALSE,
               col="red",
               pch=4,
               pos.cor=TRUE,
               xlab="MERRA SD (Normalised)",
               ylab="RCA4 runs SD (normalised)",
               main="Taylor Diagram",
               show.gamma=TRUE,
               ngamma=3,
               sd.arcs=1,
               ref.sd=TRUE,
               grad.corr.lines=c(0.2,0.4,0.6,0.8,0.9),
               pcex=1,cex.axis=1,
               normalize=TRUE,
               mar=c(5,4,6,6),
               lwd=10,
               font=5,
               lty=3)

legend(1.5,1.5,cex=1.2,pt.cex=1.2,legend=c("volcano"),pch=4,col=c("red"))



