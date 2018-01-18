# rnaseq.R
# plots related to RNAseq data

library(GenomicAlignments)
library(GenomicRanges)
library(ggplot2)
library(grid)
#library(gtable) # for combining plots vertically

#' Make a 'sashimi plot': coverage plus splice events
#' 
#' @param aln a named list of \code{GAlignments} objects
#' @param tx a \code{GRanges} object, preferably as returned by \code{import.gff3()}
#' @param meta a dataframe of extra metadata; linked to samples by column \code{"iid"}, and column
#' 	\code{"reads"} is used to normalize coverage, if available
#' @param smooth integer; if >0, size of windows in which to calculate smoothed coverage estimate
#' @param min.coverage only show regions with coverage strictly greater than this
#' @param min.splices only show splice events with multiplicity strictly greater than this
#' @param max.coverage truncate coverage at this value
#' @param log.coverage logical; if \code{TRUE}, show coverage in log10 scale
#' @param splice.scale numeric vector of length 2 used for drawing splice-junction arcs; just play with it to find good values
#' @param colours a named vector of colours used for coverage plots; grey/black is default
#' 
#' @value a \code{grid} object with the completed plot
#' 
sashimiplot <- function(aln, tx, meta = NULL, smooth = 0, min.coverage = 0, min.splices = 0, max.coverage = Inf, log.coverage = FALSE,
						splice.scale = c(1.5e3, 5), colours = NULL, colour.by = NULL, title = NULL, ...) {
	
	if (!inherits(aln[[1]], "GAlignments"))
		stop("Must supply a list of GAlignments objects.")
	
	## set boundaries of region to plot: approx. the range of all transcripts
	message("Calculating boundaries of region of interest...")
	zoom <- reduce(tx, min.gapwidth = 100e4)
	
	if (smooth > 0) {
		## calculated smoothed coverage
		message("Smoothing coverage in windows of size ", smooth, " bp...")
		bins <- seq( start(zoom), end(zoom), smooth )
		bins.gr <- GRanges(seqnames = as.vector(seqnames(zoom))[1],
						ranges = IRanges(start = bins[ -length(bins) ], end = bins[-1]-1))
		cvg <- ldply(aln, function(x) {
			transform(as.data.frame(bins.gr), score = countOverlaps(bins.gr, x))
		})
		#cvg <- ldply(bins, as.data.frame)
		colnames(cvg)[1] <- "iid"
		cvg$panel <- cvg$iid
	}
	else {
		## calculate raw coverage
		message("Estimating read coverage...")
		cvg <- lapply(lapply(aln, coverage), as, "GRanges")
		cvg <- ldply(cvg, as.data.frame)
		colnames(cvg)[1] <- "iid"
		cvg$panel <- cvg$iid
	}
	
	## discover splice junctions
	message("Discovering splice junctions...")
	sj <- lapply(aln, summarizeJunctions)
	## ignore splices to out-of-range places
	sj <- lapply(sj, subsetByOverlaps, zoom, type = "within")
	## make swoop shapes representing splicing events
	sj <- ldply(sj, as.data.frame)
	swoops <- make.swoop(sj, scale = splice.scale)
	swoops$panel <- swoops$.id
	
	## add metadata, if any
	if (!is.null(meta))
		cvg <- merge(cvg, meta, all.x = TRUE)
	
	## rescale by total coverage, if provided
	if (!is.null(meta$reads))
		cvg$depth <- with(cvg, score/(reads/1e7))
	else
		cvg$depth <- cvg$score
	
	## set colours for coverage panels, if not provided
	if (is.null(colours)) {
		colours <- rep("darkgrey", length(aln))
		names(colours) <- names(aln)
	}
	
	## set x-axis limits to cover region of interest
	xlims <- range(as.vector(ranges(zoom)))
	
	## draw coverage plot with splice events
	cvg <- subset(cvg, depth > min.coverage)
	cvg$depth <- pmin(cvg$depth, max.coverage)
	swoops <- subset(swoops, score > min.splices)
	#swoops$scolor <- swoops$row-rep(c(-1, cumsum(table(swoops$panel)/100)[-4])+1, table(swoops$panel))
	
	#return(swoops)
	#return(cvg)
	
	#NFO <- swoops %>% group_by(.id, seqnames, seqnames, start, end) %>% summarise_all(mean)
	message("Rendering plots...")
	p0 <- ggplot(cvg) +
	    geom_rect(aes(xmin = start-1, xmax = end+1, ymin = 0, ymax = depth, fill = iid)) +
	    
		scale_x_continuous(limits = xlims) +
		scale_fill_manual(values = colours) +
		#scale_colour_manual(values = colours) +
		guides(fill = FALSE, colour = FALSE) +
		facet_grid(panel ~ .) +
		ylab("coverage (reads per 10 million)\n") + xlab("") + theme_bw() +
		theme(axis.text.x = element_blank(), axis.title.x = element_blank()) #+ theme_gbrowse() 
	if(nrow(swoops)>0) {
	    p0 <- p0 + geom_line(
	        data = swoops, 
	        aes(x = x, y = y, group = row, colour = as.factor(row %% 5) ), 
	        size=1.1
	    ) +
	        geom_label_repel(
	            data = swoops[round(seq(1, nrow(swoops)-100, length.out = nrow(swoops)/100)+50),], 
	            aes(x = x, y = y, colour = as.factor(row %% 5), label =score)
	        ) 
	}
	if (log.coverage)
		p0 <- p0 + scale_y_continuous("coverage", trans = "log10")
	
	if (!is.null(title))
	    p0 <- p0 + ggtitle(title)
	
	#return(p0)
	## draw transcript ideograms
	p1 <- ggplot() +
		geom_tx(tx[ tx$type == "exon" ], at = 0, height = 5, fill = "grey50", colour.by = colour.by, ...) +
		scale_x_continuous("\nposition (Mb)", limits = xlims, labels = function(x) x/1e6) +
		scale_fill_manual(values = c("grey50", "black", "blue")) +
		facet_grid(panel ~ .) +
		guides(fill = FALSE) +
		ylab("") + theme_bw() +
		 theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
	
	## combine plots
	rez <- gtable:::rbind_gtable( ggplotGrob(p0), ggplotGrob(p1), size = "first" )
	
	## done
	return(rez)
	
}

#' Use Bezier curves to make 'swoop' shapes represnting splice junctions
make.swoop <- function(df, scale = c(1.5e3, 5), ...) {
	
	df$row <- seq_len(nrow(df))
	ddply(df, .(row), function(d) {
		w <- with(d, end-start)[1]
		mid <- c(0.2,0.8)*w+d$start[1]
		x <- c(d$start[1], mid, d$end[1])
		y <- c(0, -1*rep(pmax(w/scale[1], scale[2]), 2), 0)
		rez <- as.data.frame(Hmisc::bezier(x, y))
		for (c in colnames(d)) {
			rez[,c] <- d[1,c]
		}
		return(rez)
	})
	
}

#' Render a GRanges of exons into a transcript ideogram
#' 
#' @param exons a \code{GRanges} containing exons; metadata column \code{"Parent"} groups them into distinct transcripts
#' @param at vertical offset for first transcript
#' @param height vertical size of each transcript
#' @param arrows if arrows desired along intron connectors (to show strand orientation), this should be a call
#' 	to \code{grid::arrow()} returning the desired sort of arrow
#' @param do.introns if \code{TRUE}, draw lines through introns to connect exons
#' @param colour line colour for intron connectors
#' @param colour.by name of metadata column to use for exon boxes
#' @param fill fill colour for exon boxes, if \code{colour.by} not specified
#' @param stroke line colour for borders of exon boxes
#' 
#' @value a list of \code{ggplot2} geoms, which can be added to an existing plot with \code{`+`}.
#' 
#' @details Designed for use with the output of \code{rtracklayer::import.gff3()} with GFFs from Ensembl.  For
#' 	other use cases, mileage may vary.
#' 
geom_tx <- function(exons, at = 0, height = 1, arrows = grid::arrow(length = unit(4, "pt"), type = "closed"),
					do.introns = TRUE, colour.by = NULL, stroke = NA, colour = "black", fill = "black",
					label = FALSE, label.size = 3, drop = TRUE, ...) {
	
	.make.tx.df <- function(ex) {
		
		if (!length(ex))
			return(NULL)
		
		if (drop)
			values(ex)$at <- as.numeric(ex$parent)*height + at
		else
			values(ex)$at <- as.numeric(ex$parent)
		values(ex)$height <- height
		if (drop)
			exx <- droplevels(as.data.frame(ex))
		else
			exx <- as.data.frame(ex)
		exx$panel <- "transcripts"
		if (is.null(colour.by) || is.na(colour.by)) {
			rez <- list( geom_rect(data = exx, aes(xmin = start, xmax = end, ymin = at-height*0.4, ymax = at+height*0.4),
								   fill = fill, colour = stroke, ...) )
		}
		else {
			exx$fill <- exx[ ,colour.by ]
			rez <- list( geom_rect(data = exx, aes(xmin = start, xmax = end, ymin = at-height*0.4, ymax = at+height*0.4, fill = fill),
								   colour = stroke, ...) )
		}
		
		if (label) {
			labels <- data.frame(xpos = max(exx$end) + 1000, ypos = exx$at[1]+exx$height[1]*0.4, label = exx$parent[1], panel = exx$panel[1])
			rez[[3]] <- geom_text(data = labels, aes(x = xpos, y = ypos, label = label), size = label.size, hjust = 0)
		}
		
		if (length(ex) > 1) {
			## get exons from introns
			if (min(start(ex)) > 0)
				.introns <- GenomicRanges::gaps(ex)[-1]
			else
				.introns <-  GenomicRanges::gaps(ex)
			if (drop)
				introns <- droplevels(as.data.frame(.introns))
			else
				introns <- as.data.frame(.introns)
			#message(nrow(introns), " introns...")
			
			if (any(strand(ex) == "-")) {
				x <- introns$start
				introns$start <- introns$end
				introns$end <- x
			}
			
			#print(introns)
			introns$at <- ex$at[1]
			introns$height <- height
			introns$panel <- "transcripts"
			#if (drop)
			#	introns$middle <- with(introns, at+height*0.8/2)
			#else
			#	introns$middle <- with(introns, at + as.numeric(parent))
			rez <- c(rez, geom_segment(data = introns, aes(x = start, xend = end, y = at, yend = at),
									   arrow = arrows, colour = colour, ...))
			
		}
		
		return(rez)
		
	}
	
	if (!inherits(exons, "GRanges"))
		exons <- makeGRangesFromDataFrame(exons, keep.extra.columns = TRUE)
	
	## if input was GFF3, do everything on a per-tx basis
	if (!is.null(exons$Parent)) {
		
		if (is.null(exons$parent)) {
			## rtracklayer::import.gff3() uses list format for Parent field; assume 1 parent, listed first
			if (inherits(exons$Parent,"List"))
				exons$parent <- sapply(exons$Parent, "[", 1)
			else
				exons$parent <- as.character(exons$parent)
		}
		else {
			## respect existing 'parent' column
		}
		## 'parent' needs to be a factor; force this to be so
		if (!is.factor(exons$parent))
			exons$parent <- reorder(factor(exons$parent), start(exons), min)
		
	}
	else {
		exons$parent <- 1
	}
	
	exons <- exons[ order(exons$parent) ]
	exl <- split(exons, exons$parent)
	lapply(exl, .make.tx.df)
	
}