# This script calculates the full LD matrix within a given region or for a given variant

# gds.file <- "/Users/tmajaria/Documents/projects/public_workflows/genesis_wdl/inputs/1KG_phase3_subset.gds"
# sample.ids.file <- NULL
# ref.var <- "1:96922279"
# interval <- NULL
# half.interval <- 1000000
# min.mac <- 5
# max.mac <- 1000
# min.maf <- 0.00001
# max.maf <- 0.5
# ld.method <- "r"
# out.pref <- "ld_interval_testing"

# Load packages
packages <- c("data.table","SeqArray", "GenomicRanges", "SNPRelate", "dplyr", "tidyr", "SeqVarTools")
lapply(packages, library, character.only = TRUE)

# Parse inputs
input_args <- commandArgs(trailingOnly=T)

gds.file <- input_args[1]
sample.ids.file <- if (input_args[2] == "NA") NULL else input_args[2]
ref.var <- if (input_args[3] == "NA") NULL else input_args[3]
interval <- if (input_args[4] == "NA") NULL else unlist(strsplit(input_args[4], "[[:punct:]]"))
half.interval <- as.numeric(input_args[5])
min.mac <- as.numeric(input_args[6])
max.mac <- as.numeric(input_args[7])
min.maf <- as.numeric(input_args[8])
max.maf <- as.numeric(input_args[9])
ld.method <- tolower(input_args[10])
out.pref <- input_args[11]
visual.bool <- ifelse(startsWith(tolower(input_args[12]), "t"), TRUE, FALSE)


## # these are from the DCC pipeline, credit -> S. Gogarten https://github.com/UW-GAC/analysis_pipeline/blob/master/TopmedPipeline/R/filterVariants.R

.calcMAF <- function(gds) {
  ref.freq <- alleleFrequency(gds)
  pmin(ref.freq, 1-ref.freq)
}

.calcMAC <- function(gds) {
  ref.cnt <- alleleCount(gds)
  n.obs <- SeqVarTools:::.nSampObserved(gds)
  round(pmin(ref.cnt, 2*n.obs - ref.cnt))
}

.variantDF <- function(gds) {
  data.frame(variant.id=seqGetData(gds, "variant.id"),
             chr=seqGetData(gds, "chromosome"),
             pos=seqGetData(gds, "position"),
             ref=refChar(gds),
             alt=altChar(gds),
             nAlleles=seqNumAllele(gds),
             maf=.calcMAF(gds),
             mac=.calcMAC(gds),
             stringsAsFactors=FALSE)
}
.expandAlleles <- function(gds) {
  .variantDF(gds) %>%
    separate_rows_("alt", sep=",") %>%
    group_by_("variant.id") %>%
    mutate_(allele.index=~1:n()) %>%
    as.data.frame()
}

.ldPair <- function(var1, gds, samples, var2, method){
  ld <- snpgdsLDMat(gds, method = method, slide = 0, sample.id = samples, snp.id = c(var1, var2))
  if (method != 'dprime'){
    ld.df <- ld$LD * ld$LD  
  } else {
    ld.df <- ld$LD
  }
  
  row.names(ld.df) <- colnames(ld.df) <- ld$snp.id
  return(ld.df)
}

.errorDict <<- list(
  "samples" = "No samples remain after filtering by sample list. Check that sample list matches samples in GDS file.",
  "refvar" = "Reference variant not found in GDS file.",
  "interval" = "Interval not found in GDS file.",
  "variants" = "No variants remain after filtering by mac and maf thresholds."
  )

.exitError <- function(why, out.file){
  seqClose(gds.data)
  write.csv(.errorDict[[why]], file = out.file, row.names = T, col.names = NA, sep = ",", quote = F) 
  stop(.errorDict[[why]])
}

.ldMethodName <- function(ld.method){
  n <- list(
    'r' = expression(R^2),
    'composite' = "Composite",
    'corr' = "Correlation",
    'dprime' = "D\'",
    'cov' = "Cov"
    )
  n[[ld.method]]
}

plotLD <- function(ld, ld.method, single.var, out.file){
  library(ggplot2)
  library(reshape2)
  ld.melt <- melt(as.matrix(ld))
  ld.melt$Var1 <- factor(ld.melt$Var1, levels = rev(row.names(ld)))
  png.dims <- c(max(2, 0.1*ncol(ld)), max(2, 0.1*nrow(ld)))
  if (single.var){
    plot <- ggplot() +
      geom_raster(data = ld.melt, aes(x=Var2, y=Var1, fill=value)) +
      scale_fill_gradient2(low = 'white', high = 'red') +
      scale_x_discrete(position = "top") +
      theme(axis.title.x = element_blank()) +
      theme(axis.title.y = element_blank()) +
      theme(axis.text.x = element_text(angle = 90, size = 6)) +
      theme(axis.text.y = element_text(size = 6)) +
      # theme(legend.position = "right") +
      labs(fill = .ldMethodName(ld.method)) +
      guides(fill = guide_colourbar(barwidth = 0.5, barheight = 3)) + 
      coord_fixed()
      # theme(legend.position = "none")

  
    png(filename = out.file, width = png.dims[1], height = 2.5, units = "in", res=400, type = "cairo")
  } else {
    plot <- ggplot() +
      geom_raster(data = ld.melt, aes(x=Var2, y=Var1, fill=value)) +
      scale_fill_gradient2(low = 'white', high = 'red') +
      scale_x_discrete(position = "top") +
      theme(axis.title.x = element_blank()) +
      theme(axis.title.y = element_blank()) +
      theme(axis.text.x = element_text(angle = 90, size = 6)) +
      theme(axis.text.y = element_text(size = 6)) +
      labs(fill = .ldMethodName(ld.method)) +
      guides(fill = guide_colourbar(barwidth = 0.5, barheight = 3)) +
      coord_fixed()
      # theme(legend.position = "none")
    png(filename = out.file, width = png.dims[1], height = png.dims[2]+2, units = "in", res=400, type = "cairo")
  }
  
  print(plot)
  dev.off()
}

##

# Show the input arguments
print(paste0("gds.file = ", gds.file, " (", input_args[1], ")"))
print(paste0("sample.ids.file = ", sample.ids.file, " (", input_args[2], ")"))
print(paste0("ref.var = ", ref.var, " (", input_args[3], ")"))
print(paste0("interval = ", interval, " (", input_args[4], ")"))
print(paste0("half.interval = ", half.interval, " (", input_args[5], ")"))
print(paste0("min.mac = ", ref.var, " (", input_args[6], ")"))
print(paste0("max.mac = ", max.mac, " (", input_args[7], ")"))
print(paste0("min.maf = ", min.maf, " (", input_args[8], ")"))
print(paste0("max.maf = ", max.maf, " (", input_args[9], ")"))
print(paste0("ld.method = ", ld.method, " (", input_args[10], ")"))
print(paste0("out.pref = ", out.pref, " (", input_args[11], ")"))


### Processing inputs ###
# check inputs
if (all(is.null(c(ref.var, interval)))){
  stop("Need to set either reference variant or interval")
} 

# Open gds file and get sample ids
gds.data <<- seqOpen(gds.file)
gds.samples <- seqGetData(gds.data, "sample.id")

# Load sample ids
if (!is.null(sample.ids.file)){
  sample.ids <- read.table(sample.ids.file, stringsAsFactors = F, header = F)$V1
  sample.ids <- sample.ids[sample.ids %in% gds.samples]
  if (length(sample.ids) == 0) .exitError("samples", paste0(out.pref, ".ERROR.csv"))
  seqSetFilter(gds.data, sample.id=sample.ids)
  gds.samples <- seqGetData(gds.data, "sample.id")
}

# parse ref var
if (!is.null(ref.var)){
  ref.c <- unlist(strsplit(ref.var, "[[:punct:]]"))
  ref.chr <- ref.c[1]
  ref.pos <- as.numeric(ref.c[2])
}

# parse interval
if (is.null(interval)){
  interval.chr <- ref.chr
  interval.start <- max(0,ref.pos-half.interval)
  interval.end <- ref.pos+half.interval
} else {
  interval.chr <- interval[1]
  interval.start <- as.numeric(as.character(interval[2]))
  interval.end <- as.numeric(as.character(interval[3]))
}

# make sure that we get the right ld method
if (!(tolower(ld.method) %in% c("comp","composite","r","d","dprime","d'","corr","correlation","cov"))){
  warning("LD method not recognised, defaulting to correlation. LD method must be one of 'comp','r','d','corr','cov'. See SNPRelate documentation for snpgdsLDMat.")
  ld.method <- "corr"
}

# fix ld method name
if (ld.method == 'comp') ld.method <- 'composite'
if (ld.method %in% c("d","dprime","d'")) ld.method <- 'dprime'
if (ld.method == 'correlation') ld.method <- 'corr'

# create output prefix
out.file <- paste0(paste(out.pref, interval.chr, interval.start, interval.end, sep = "."), ".csv")
print(paste0("out.file = ", out.file))

# make sure that interval and gds have the same chr format
gds.chrs <- unique(seqGetData(gds.data,"chromosome"))
gds.chr <- gds.chrs[1]

if (startsWith(gds.chr,"chr") && !(startsWith(interval.chr, "chr"))){
  interval.chr <- sub("^","chr",interval.chr)
} else if (!(startsWith(gds.chr,"chr")) && startsWith(interval.chr, "chr")){
  interval.chr <- sub("chr","",interval.chr)
}

# check that reference variant chromosome in gds file
if (!is.null(ref.var)){
  if (!(ref.chr %in% gds.chrs)) .exitError("refvar", paste0(out.pref, ".ERROR.csv"))
}

# check that interval chr is in gds file
if (!(interval.chr %in% gds.chrs)) .exitError("interval", paste0(out.pref, ".ERROR.csv"))

# filter gds to only those variants in our interval
seqSetFilterChrom(gds.data, interval.chr, from.bp=interval.start , to.bp=interval.end)

# check that we still have variants
if (length(seqGetData(gds.data, "variant.id")) == 0) .exitError("variants", paste0(out.pref, ".ERROR.csv"))

# filter gds by mac and maf
seqSetFilterCond(gds.data, mac=c(min.mac,max.mac), maf=c(min.maf,max.maf))

# check that we still have variants
if (length(seqGetData(gds.data, "variant.id")) == 0) .exitError("variants", paste0(out.pref, ".ERROR.csv"))

# prepare marker output
markers <- .expandAlleles(gds.data)
markers$MarkerName <- paste(markers$chr, markers$pos, markers$ref, markers$alt, sep=":")

# make sure that the reference variant is in markers
if (!is.null(ref.var)){
  if (!any(grepl(ref.var, markers$MarkerName))){
    seqResetFilter(gds.data)
    seqSetFilterChrom(gds.data, ref.chr, from.bp=ref.pos-1 , to.bp=ref.pos+1)
    
    # check that the reference variant is in gds file
    if (length(seqGetData(gds.data, "variant.id")) == 0) .exitError("refvar", paste0(out.pref, ".ERROR.csv"))

    ref.df <- .expandAlleles(gds.data)
    ref.df$MarkerName <- paste(ref.df$chr, ref.df$pos, ref.df$ref, ref.df$alt, sep=":")
    ref.df <- ref.df[grep(ref.var, ref.df$MarkerName),][1,]
    markers <- rbind(markers, ref.df)
  }
}

markers <- markers[!duplicated(markers$variant.id),]
row.names(markers) <- markers$variant.id

#########################

### Generate LD Matrix ###
if (is.null(ref.var)){
  ld <- snpgdsLDMat(gds.data, method = ld.method, slide = 0, sample.id = gds.samples, snp.id = markers$variant.id)
  seqClose(gds.data)
  
  if (ld.method != 'dprime'){
    ld.df <- data.frame(ld$LD * ld$LD)
    # ld.df <- data.frame(ld$LD)
  } else {
    ld.df <- data.frame(ld$LD)
  }
  row.names(ld.df) <- colnames(ld.df) <- markers[as.character(ld$snp.id),]$MarkerName
} else {
  ref.id <- row.names(markers[grep(ref.var, markers$MarkerName),])[1]
  ld.list <- lapply(markers$variant.id[markers$variant.id != ref.id], .ldPair, gds.data, gds.samples, ref.id, ld.method)
  ld.list <- unlist(lapply(ld.list, function(x) x[as.character(ref.id), ]))
  ld.list <- ld.list[!duplicated(names(ld.list))]
  ld.df <- t(data.frame(ld.list[order(names(ld.list))]))
  row.names(ld.df) <- markers[ref.id,]$MarkerName
  colnames(ld.df) <- markers[colnames(ld.df),]$MarkerName
}

# make visulaization
viz.file <- sub(".csv$", ".png", out.file)
if (visual.bool) {
  single.var <- ifelse(is.null(ref.var), FALSE, TRUE)
  plotLD(ld.df, ld.method, single.var, viz.file)
} else {
  file.create(viz.file)
}

# write only vector of LD values if given reference variant
# otherwise write entire matrix
write.csv(ld.df, file = out.file, row.names = T, col.names = NA, sep = ",", quote = F)  




