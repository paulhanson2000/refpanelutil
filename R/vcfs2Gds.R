vcfs2Gds <- function(files, output_name,
                     regions=NULL, variant_ids=NULL, sample_ids=NULL,
                     chr_nm_map=NULL, exclude_annos=NULL, extra_cmds=NULL,
                     scratch_dir="/tmp/", bcftools_bin="bcftools") {
  #' Vcfs2Gds
  #'
  #' Convert a bunch of VCF files into one SeqArray GDS file.
  #'
  #' @param files VCF file path(s) or URL(s).
  #' @param output_name Filename of output GDS file. 
  #' @param storage.option <TODO, not implemented> Storage/Compression options. Passed directly to seqVCF2GDS. See ?seqVCF2GDS and ?seqStorageOption for details.
  #' @param sample_ids Character vector of sample IDs, to subset the VCFs by.
  #' @param variant_ids Character vector of variant IDs, to subset the VCFs by.
  #' @param regions A data.frame of bed-like format (1st column chrom, 2nd chromStart, 3rd chromEnd, other cols unused)
  #' @param chr_nm_map A data.frame with two columns. Used to convert the VCF files' chromosome/contig names from the names in the first column to the second. E.g. to convert from "chr1" to "1".
  #' @param exclude_annos <TODO, see bcftoolsPipableCmd(), maybe redundant w/ info_fields and format_fields params, choose one>
  #' @param extra_cmds <TODO, document>
  #' @param bcftools_bin Path to the program "bcftools".
  #' @section A tip to increase speed: Even if you already provide variant_ids, also providing regions will speed things up a lot. It is easier for computers to compare chromosome and position numbers than ID strings.
  #' @section TODO: Support and chr:pos(:ref:alt) IDs? Support mixed types of IDs? Infer regions if all IDs are chr:pos(:ref:alt), and maybe mention in these docs that that optimization is made once implemented. Deal with "chr#" vs. "#" syntax for the user of not matching (I think bcftools norm or  bcftools annotate --rename-chrs can help). Allow GRanges regions input? Support BCF files too?--Note would have to use pipe(..., "rb") maybe.

  bcftools_cmd <- bcftoolsPipableCmd(files, output_type="v",
                                     regions=regions, variant_ids=variant_ids, sample_ids=sample_ids,
                                     chr_nm_map=chr_nm_map, exclude_annos=exclude_annos, extra_cmds=extra_cmds,
                                     scratch_dir=scratch_dir, bcftools_bin=bcftools_bin)
  vcf_connection <- pipe(bcftools_cmd, "rt")

  # TODO: add a "safe" option to download and convert one VCF at a time, seqMerge()ing them. This would prevent having to start completely over again if the download corrupts in one place.

  f <- SeqArray::seqVCF2GDS(vcf_connection, output_name, storage.option="ZIP_RA", ignore.chr.prefix="")
  close(vcf_connection)
  unlink(c(paste0(basename(files),".tbi"), paste0(basename(files),".csi"))) # Index files clutter the working directory when accessing over the internet.
  
  # TODO: Commented out b/c it takes really long and ZIP is small enough.
  #seqRecompress(f, "LZMA") # See ?seqRecompress -- seqVCF2GDS uses a LOT of memory when given a high compression option. Better to use a low compression, and then recompress later.
  f
}
