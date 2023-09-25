# Calculate LD from a SeqArray GDS file.
ldCalc <- function(gds,
                   regions=NULL, variant_ids=NULL, sample_ids=NULL,
                   ref_alleles=NULL,
                   maf_thresh=NaN, mac_thresh=1L, variant_missingness_thresh=NaN,
                   slide=0, method="corr",
                   n_threads=parallel::detectCores(), yes_really=F) {
  #' ldCalc
  #'
  #' Calculate LD from `SeqArray` GDS files.
  #'
  #' @param gds A `SeqArray` GDS file handle (obtained from `seqOpen("example.gds")`)
  #' @param regions Character vector of regions of the format chr:start-end, or a data.frame of bed-like format (1st column chromosome, 2nd start position (0-indexed), 3rd end, other cols unused)
  #' @param variant_ids Character vector of variant IDs. If given, the resulting LD matrix's colnames will be variant IDs. Otherwise, the colnames will be of the format chr:pos:ref:alt.
  #' @param sample_ids Character vector of sample IDs, to subset the VCF/BCF(s) by.
  #' @param ref_alleles Character vector of variant reference alleles corresponding to variant_ids.
  #' @param maf_thresh MAF threshold
  #' @param mac_thresh MAC threshold
  #' @param variant_missingness Missingness threshold
  #' @param slide Size of sliding window. Full LD matrix will be computed by default.
  #' @param method Measure of LD to use. "corr" by default, referring to Pearson r. See `?snpgdsLDMat` from the `SNPRelate` package for more details.
  #' @param n_threads Number of threads for parallel processing.
  #' @param yes_really Safety measure. LD for a large number of variants can consume a lot of hard drive space. If calculating LD for a very large number of variants (>20k), pass yes_really=T to confirm you wish to do this.

  if(class(gds) != "SeqVarGDSClass")               stop("gds is class ", class(gds), ", not SeqVarGDSClass as it should be.")
  if(is.null(ref_alleles) & !is.null(variant_ids)) stop("ldCalc warning: if providing variant_ids, it is highly recommended to also provide ref_alleles. If the reference alleles (i.e. non-effect alleles) in your GWAS data =/= those in the reference panel, the LD will be incorrect.")

  if(length( sample_ids)==1 && file.exists( sample_ids))   sample_ids <- readLines( sample_ids)
  if(length(variant_ids)==1 && file.exists(variant_ids))  variant_ids <- readLines(variant_ids)
  if(length(    regions)==1 && file.exists(    regions))      regions <- read.table(   regions)
  if(class(regions) %in% c("data.frame", "data.table" ))      regions <- paste0(regions[,1],':',regions[,2],'-',regions[,3])
  # TODO: add GRanges support

  if(!is.null(variant_ids)) if(anyDuplicated(seqGetData(gds,"annotation/id")) > 0) stop("Error: detected duplicates in gds file's \"annotation/id\" field.")
  if(anyDuplicated(seqGetData(gds,    "sample.id")) > 0) stop("Error: detected duplicates in gds file's \"sample.id\" field.")
  if(anyDuplicated(variant_ids)                     > 0) stop("Error: detected duplicates in variant_ids.")
  if(anyDuplicated( sample_ids)                     > 0) stop("Error: detected duplicates in sample_ids.")


  seqFilterPush(gds) # Save user's pre-existing filters to restore it afterwards. 

  # Filter by regions
  if(!is.null(regions)) {
    chrs <-            sapply(regions, function(reg) strsplit(reg,":|-")[[1]])[1,]
    begs <- as.numeric(sapply(regions, function(reg) strsplit(reg,":|-")[[1]])[2,])
    ends <- as.numeric(sapply(regions, function(reg) strsplit(reg,":|-")[[1]])[3,])
    seqSetFilterChrom(gds, chrs, from.bp=begs, to.bp=ends, intersect=T)
  }

  # Filter by sample_ids
  if(!is.null(sample_ids)) {
    message(length(sample_ids) - sum(seqGetData(gds,"sample.id") %in% sample_ids), "/", length(sample_ids), " sample_ids will be omitted because they were not in the gds file.")
    seqSetFilter(gds, sample.id = sample_ids, action="intersect")
  }

  # Filter by MAF/MAC/missingness
  # Note seqSetFilterCond implicitly intersects existing filters. (Interestingly seqSetFilterChrom does not.)
  seqSetFilterCond(gds, maf=maf_thresh, mac=mac_thresh, missing.rate=variant_missingness_thresh)

  # Filter by variant_ids
  if(!is.null(variant_ids)) {
    message(length(variant_ids) - sum(seqGetData(gds,"annotation/id") %in% variant_ids), "/", length(variant_ids), " variant_ids will be omitted because they were not in the gds file.")
    seqSetFilter(gds, variant.sel = seqGetData(gds,"annotation/id") %in% variant_ids, action="intersect")
  }


  # Must convert to SNP GDS format (an older version of SeqArray GDS format) to use snpgdsAlleleSwitch.
  tmp_filename <- "/tmp/snpgds_file.gds"
  seqGDS2SNP(gds, tmp_filename, compress.geno="", compress.annotation="")

  seqFilterPop(gds) # Restore user's pre-existing filter


  snpgds <- snpgdsOpen(tmp_filename, readonly=F)
  tryCatch({ # So that the user is not left with a dangling snpgds file pointer
      snpgdsGetData <- function(snpgds,field) read.gdsn(index.gdsn(snpgds,field))

      n_variants <- length(snpgdsGetData(snpgds,"snp.id"))
      if(n_variants > 20000 && !yes_really) stop("You're about to compute an LD matrix for ", n_variants, " variants, are you sure? Run ldCalc again with yes_really=TRUE if so.") # TODO: estimate final file size for style poitns

      # Allele switching
      if(!is.null(ref_alleles) & !is.null(variant_ids)) {
        # Given IDs/alleles and those in the GDS file not in same order. Match the ref alleles to match the GDS file.
        ref_alleles <- ref_alleles[match(snpgdsGetData(snpgds,"snp.rs.id"), variant_ids)]
        snpgdsAlleleSwitch(snpgds, toupper(ref_alleles))
      }

      # LD
      ld <- snpgdsLDMat(snpgds, slide=slide, method=method, num.thread=n_threads)

      # If user filtered using variant_ids, they probably want IDs as the colnames. Otherwise, chr:pos:ref:alt.
      if(!is.null(variant_ids)) {
        colnames(ld$LD) <- snpgdsGetData(snpgds, "snp.rs.id") # Not actually necessarily rs IDs. "snp.rs.id" is just what "annotation/id" becomes after SeqArray -> SNP GDS conversion.
      } else {
        colnames(ld$LD) <- mapply( USE.NAMES=F,
          snpgdsGetData(snpgds, "snp.chromosome"),
          snpgdsGetData(snpgds, "snp.position"),
          strsplit(snpgdsGetData(snpgds, "snp.allele"),'/'),
          FUN = function(chr,pos,alleles) paste0(chr,':',pos,':',alleles[[1]],':',alleles[[2]])
        )
      }
      
    },

    error = function(err) {
      snpgdsClose(snpgds)
      unlink(tmp_filename)
      stop()
  })
  snpgdsClose(snpgds)
  unlink(tmp_filename)

  return(ld$LD)
}
