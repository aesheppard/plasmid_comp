#!/usr/bin/Rscript

library(gplots)
library(argparse)

parser = ArgumentParser()
parser$add_argument('--matfile', help = 'Data file containing matrix of values to be plotted, in csv format.', required = TRUE)
parser$add_argument('--pdffile', help = 'Output pdf file name,', required = TRUE)
parser$add_argument('--reorder', help = 'Flag indicating whether to reorder samples. Default is to plot in the same order as the input file.', action = 'store_true')
parser$add_argument('--orderfile', help = 'File to write list of samples, in the order plotted.')
parser$add_argument('--scale', help = 'Multiplication factor for scale labels, i.e. bin size that was used.', default = 1)
parser$add_argument('--scaleinc', help = 'Scale bar increments, in kb', default = 10)
parser$add_argument('--contigs', help = 'File containing list of contig names and lengths. Default assumes single contig.')
parser$add_argument('--contigcol', help = 'Colour to use for indicating contig breaks.', default = 'black')
parser$add_argument('--singlescale', help = 'Use a single scale instead of restarting for each contig', action = 'store_true')
parser$add_argument('--cexRow', help = 'Character size adjustment for row labels.', default = 1)
parser$add_argument('--pdfExp', help = 'Expansion factor for pdf dimensions.', default = 1)
parser$add_argument('--rainbowstart', help = 'Value at which to start rainbow colours, or lower threshold if single colour (gradient) used', default = 80)
parser$add_argument('--singlecol', help = 'Single colour to use instead of rainbow, e.g. blue, black, etc. Must be a valid R colour. May or may not be used as gradient (see below).')
parser$add_argument('--nogradient', help = 'Single colour indicating presence above the given threshold, rather than a true heatmap. Ignored unless --singlecol also used.', action = 'store_true')
args = parser$parse_args()

blasts = read.csv(args$matfile, header = TRUE)
rnames = blasts[,1]
mat_data = data.matrix(blasts[,2:ncol(blasts)])
rownames(mat_data) = rnames

numvalues = 1000
proprain = 1 - (as.numeric(args$rainbowstart) / 100.0)
numrain = proprain * numvalues
if (!is.null(args$singlecol)) {
	belowthresh = rep('white', numvalues - numrain)
	if (args$nogradient) {
		rainb = rep(args$singlecol, numrain)
	} else {
		rainb = colorRampPalette(colors = c('white', args$singlecol))(numrain)
	}
} else {
	rainb = rainbow(numrain * 1.25)[1:numrain]
	belowthresh = c(colorRampPalette(colors = c('white', rainb[1]))(numvalues - numrain))
}
heatcols = c(belowthresh, rainb)

pdfdims = as.numeric(args$pdfExp) * 7
pdf(args$pdffile, width = pdfdims, height = pdfdims)

labspace = as.numeric(args$scaleinc) * 1000 / as.numeric(args$scale)

contigs = read.delim(args$contigs, header = FALSE, col.names = c('name', 'size'))
sizes = as.integer(as.numeric(contigs$size) / as.numeric(args$scale))

indices = c()
total = 0
labels = c()
for (size in sizes) {
    if (size > 0) {
    	for (i in 1:size) {
	    currpos = i
	    if (args$singlescale) {
	    	currpos = currpos + total
	    }
    	    if (currpos %% labspace == 0) {
	   	newval = currpos * as.numeric(args$scale) / 1000
	    } else {
	   	newval = ''
	    } 
    	    labels = c(labels, newval)
    	}
        total = total + size
    	indices = c(indices, total)
    }
}


dend = FALSE
if (args$reorder) {
   dend = TRUE
}

hm = heatmap.2(mat_data, dendrogram="none", col = heatcols, trace = "none", density.info = "none", Rowv = dend, Colv = FALSE, labCol = labels, srtCol = 0, adjCol = c(0.5,0), colsep = indices, sepcolor = args$contigcol, cexRow = args$cexRow, offsetCol = -0.2, offsetRow = -0.2)

if (!is.null(args$orderfile)) {
	write(rev(rownames(mat_data)[hm$rowInd]), file = args$orderfile)
}

dev.off()
