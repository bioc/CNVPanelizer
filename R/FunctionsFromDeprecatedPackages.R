# Dependencies from packages that got deprecated

countBamInGRanges <- function (bam.file, granges, min.mapq = 1, read.width = 1, stranded.start = FALSE, get.width = FALSE, remove.dup = FALSE) {
  rds.counts <- integer(length(granges))
  seq.names <- unique(as.character(seqnames(granges)))
  seq.names.in.bam <- names(scanBamHeader(bam.file)[[1]]$targets)
  for (seq.name in seq.names) {
    if (seq.name %in% seq.names.in.bam) {
      granges.subset <- granges[seqnames(granges) == seq.name]
      strand(granges.subset) <- "*"
      scan.what <- c("pos","mapq")
      if (stranded.start | get.width) {
        scan.what <- c(scan.what, "qwidth")
      }
      if (stranded.start | remove.dup) {
        scan.what <- c(scan.what, "strand")
      }
      rds <- scanBam(bam.file, param = ScanBamParam(what = scan.what, which = range(granges.subset)))
      if (min.mapq > 0) {
        mapq.test <- rds[[1]]$mapq >= min.mapq & !is.na(rds[[1]]$mapq)
      } else {
        mapq.test <- rep(TRUE,length(rds[[1]]$mapq))
      }
      if (remove.dup) {
        if (get.width) {
          # this check is fast and accurate, assuming that read widths
          # are less than 1e4 bp and positions are less than 1e12 bp
          # the precision for doubles is 16 significant digits
          mapq.test <- mapq.test & !duplicated(rds[[1]]$pos + as.numeric(rds[[1]]$strand)/10 + rds[[1]]$qwidth/1e5)
        } else {
          mapq.test <- mapq.test & !duplicated(rds[[1]]$pos + as.numeric(rds[[1]]$strand)/10)
        }
      }
      if (sum(mapq.test) > 0) {
        if (stranded.start) {
          rds.ranges <- GRanges(seq.name, IRanges(start = ifelse(rds[[1]]$strand[mapq.test] %in% c("+","*"), rds[[1]]$pos[mapq.test] , rds[[1]]$pos[mapq.test] + rds[[1]]$qwidth[mapq.test] - read.width), end = ifelse(rds[[1]]$strand[mapq.test] %in% c("+","*"), rds[[1]]$pos[mapq.test] + read.width - 1, rds[[1]]$pos[mapq.test] + rds[[1]]$qwidth[mapq.test] - 1)))
        }
        else if (get.width) {
          rds.ranges <- GRanges(seq.name, IRanges(start = rds[[1]]$pos[mapq.test], width = rds[[1]]$qwidth[mapq.test]))
        }
        else {
          rds.ranges <- GRanges(seq.name, IRanges(start = rds[[1]]$pos[mapq.test], width = read.width))
        }
        rds.counts.seq.name <- countOverlaps(granges.subset, rds.ranges)
        rds.counts[as.logical(seqnames(granges) == seq.name)] <- rds.counts.seq.name
      }
      else {
        rds.counts[as.logical(seqnames(granges) == seq.name)] <- 0
      }
    }
    else {
      rds.counts[as.logical(seqnames(granges) == seq.name)] <- 0
    }
  }
  if (sum(rds.counts) == 0) {
    warning("No reads found with minimum mapping quality")
  }
  rds.counts
}

