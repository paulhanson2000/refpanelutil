# TODO: Add built-in mac/maf/miss, useful if s/o else w/ restricted data is doing LD for you

# Calculate LD from a SeqArray GDS file.

# This function can accept files as input, and be run standalone on the command line:
  # R -q -e "source(ldCalc.R); ldCalc("myfile.gds", sample_ids="my_samples.txt", ...)
ldCalc <- function(gds,
                   sample_ids=NULL, variant_ids=NULL, regions=NULL,
                   slide=0, method="corr",
                   n_threads=parallel::detectCores(), yes_really=F) {

  #if(!dir.exists(dirname(output_name))) stop("The directory of output_name ",output_name," does not exist, please create it so output can be written there.")
  gds_was_filename <- F
  if(is.character(gds) && file.exists(gds)) { gds_was_filename <- T; gds <- seqOpen(gds, allow.duplicate=T) }
  if(class(gds) != "SeqVarGDSClass") stop("gds is class ", class(gds), " and not SeqVarGDSClass as it should be.")
  seqFilterPush(gds) # Save user's prexisting filter if they have one (good manners)
  # TODO: tryCatch to avoid hidden side effect of filter being kept on error? Here's a graceful exit function: die <- function(...) { seqFilterPop(gds); if(gds_was_filename) seqClose(gds); stop(...) } # Exit gracefully if error, popping the filter.

  if(length( sample_ids)==1 && file.exists(    sample_ids))   sample_ids <- readLines(    sample_ids)
  if(length(variant_ids)==1 && file.exists(   variant_ids))  variant_ids <- readLines(   variant_ids)
  if(length(    regions)==1 && file.exists(       regions)) {    regions <- read.table(      regions); regions <- paste0(regions[,1],":",regions[,2],"-",regions[,3]) }
  # TODO: add GRanges support

  if(anyDuplicated(variant_ids) > 0) stop("Error: detected duplicates in variant_ids.")
  is_chr_pos_ids <- any(grepl(":",variant_ids))
  is_rs_ids <- any(grepl("^rs",variant_ids))
  if(is_chr_pos_ids && is_rs_ids) stop("Error: detected a mix of rs and chr:pos:allele IDs in variant_ids. Please use one type or the other.")

  # Filter by regions
  if(!is.null(regions)) {
    chrs <- sub("chr","",strsplit(regions, ":|-")[[1]][1])
    begs <-   as.numeric(strsplit(regions, ":|-")[[1]][2])
    ends <-   as.numeric(strsplit(regions, ":|-")[[1]][3])
    seqSetFilterChrom(gds, chrs, from.bp=begs, to.bp=ends, intersect=T)
  }

  # Filter by variant_ids
  if(!is.null(variant_ids)) {
    chrpos_id_sel <- gsub("_",":", seqGetData(gds, "$chrom_pos_allele")) %in% variant_ids # Because seqGetData "$chrom_pos_allele" returns chr:pos_ref_alt not chr:pos:ref:alt
        rs_id_sel <-               seqGetData(gds,   "annotation/id"  )  %in% variant_ids
           id_sel <- rs_id_sel | chrpos_id_sel

    message(length(variant_ids) - sum(id_sel), "/", length(variant_ids), " variant_ids will be omitted because they were not in common with the gds file.")
    seqSetFilter(gds, variant.sel = id_sel, action="intersect")
  }

  n_variants <- sum(seqGetFilter(gds)$variant.sel)
  if(n_variants > 20000 && !yes_really) stop("You're about to compute an LD matrix for ", n_variants, " variants, are you sure? Run ldCalc again with yes_really=TRUE if so.") # TODO: estimate final file size for style poitns

  # Filter by sample_ids
  if(!is.null(sample_ids)) seqSetFilter(gds, sample.id = sample_ids, action="intersect")

  # TODO: Add allele switching
    # Need to sort the input variants (specifically, effect_alleles) by whatever order it is in the gds file before alleleSwitch
       # ^ Look to what you did with match() in the original ldCalc function for help
    # After sorted and stuff, seqGDS2SNP(), seqClose and unlink if made temp copy, snpgdsAlleleSwitch(), snpgdsLDMat(),

  ld <- snpgdsLDMat(gds,
                    sample.id = seqGetData(gds,  "sample.id"),
                       snp.id = seqGetData(gds, "variant.id"),
                    slide = slide, method = method,
                    num.thread = n_threads)
  colnames(ld$LD) <- seqGetData(gds,"$chrom_pos_allele")

  seqFilterPop(gds)
  if(gds_was_filename) seqClose(gds)
  return(ld$LD)
}
