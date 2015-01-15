#!/usr/bin/Rscript
source("sim-patchy-selection-fns.R")

infile <- "sims/89011-r101-101-sb0.1-sm-0.01-N100-pophistory-run.Rdata"
outfile <- gsub(".R[dD]ata",".png",infile)
png( file=outfile, width=5*288, height=3*288, res=288, pointsize=10 )
layout(matrix(1:6,nrow=2))
par(mar=c(1,1,1,1)+.1)
do.times <- c(6,8,10,12,14,16)
stepsize <- 50
for (k in do.times) {
    plotpop( pophist$pophist[,,,k], params=pophist$pop$params, plotlegend=FALSE, xlab='', ylab='', main=paste("t=",stepsize*k) )
}
dev.off()

if (FALSE) {

infiles <- if (!interactive()) {commandArgs(TRUE)} else { scan(what='char') }

for (infile in infiles) {
    (load(infile))
    # outfile <- gsub(".R[dD]ata",".pdf",infile)
    # if (!interactive()) { pdf( file=outfile, width=5, height=3, pointsize=10 ) }
    outfile <- gsub(".R[dD]ata",".png",infile)
    if (!interactive()) { png( file=outfile, width=5*144, height=3*144, res=144, pointsize=10 ) }
    layout(matrix(1:6,nrow=2))
    par(mar=c(1,1,1,1)+.1)

    do.times <- c(6,8,10,12,14,16)
    if (interactive()) { layout(1) }
    for (k in do.times) {
        plotpop( pophist$pophist[,,,k], params=pophist$pop$params, plotlegend=FALSE, xlab='', ylab='' )
        if ( interactive() && length(locator(1))==0) { break; }
    }
    if (!interactive()) { dev.off() }
}
}
