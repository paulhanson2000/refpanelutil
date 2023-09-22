liftOver <- function(to_lift, from_build, to_build, out_file=NULL, liftover_bin="liftOver") {
  #' liftOver
  #'
  #' Simple wrapper for UCSC's liftOver, to convert genome coordinates between reference genomes. Required chain files are automatically downloaded. BED file specification: http://www.genome.ucsc.edu/FAQ/FAQformat.html#format1 
  #'
  #' @param to_lift A BED format file, or a BED-like data.frame. NOTE: UCSC-style chromosome naming must be used (e.g. "chr1", not just "1"), or liftOver won't work.
  #' @param from_build A string, e.g. "hg19"
  #' @param to_build A string, e.g. "hg38"
  #' @param out_file Desired output filename. Results are also returned as an R data.frame whether this option is given or not. 
  #' @param liftover_bin Path to the program liftOver. Can be downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
  #' @return A BED-like data frame with lifted genome coordinates.
  #' @section Command-line use: Since this function can accept files, you could run it from the command line like so: R -q -e source("liftOver.R"); liftOver("my_file.bed", "hg19", "hg38", liftover_bin="path/to/liftOver"). At that point you should probably just use liftOver by itself though.

  # Find or download liftOver executable
  if(suppressWarnings(system(liftover_bin, ignore.stderr=T, ignore.stdout=T))==127) liftover_bin <- paste0("./", liftover_bin)
  if(suppressWarnings(system(liftover_bin, ignore.stderr=T, ignore.stdout=T))==127) {
    message("Could not find liftOver at path: \"",liftover_bin,"\". Please edit the liftover_bin argument to point to your liftOver executable.\n(WINDOWS USERS: there is no liftOver program for Windows, so this function will not work.)")
    yes <- "y" == readline("Automatically download liftOver from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64 and put it in the current directory? [y/n].\n")
    if(yes) {
      download.file("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver", "liftOver")
      system("chmod +x liftOver")
      liftover_bin <- "./liftOver"
    } else { stop("Not downloading liftOver.") }
  }

  # Get chain file
  chain_file_path <- dirname(liftover_bin)
  chain_file_name <- paste0(from_build, "To", tools::toTitleCase(as.character(to_build)), ".over.chain.gz")
  chain_file <- paste0(chain_file_path,"/",chain_file_name)

  if(from_build == to_build) stop("Using liftOver to convert from one reference genome to the same reference genome!")
  if(!file.exists(chain_file)) {
    err_code <- download.file(paste0("https://hgdownload.soe.ucsc.edu/goldenPath/",from_build,"/liftOver/",chain_file_name), destfile = chain_file)
    if(err_code > 0) stop("Please download the chain file https://hgdownload.soe.ucsc.edu/goldenPath/",from_build,"/liftOver/",chain_file_name," and add it to ",chain_file_path,"/.")
  }

  wrote_tmp_file <- F
  if(!is.character(to_lift)) { # User passed in a data.frame
    wrote_tmp_file <- T
    write.table(to_lift, "/tmp/to_lift.bed", sep=' ', row.names=F, col.names=F, quote=F)
    to_lift <- "/tmp/to_lift.bed"
  }
  if(!file.exists(to_lift)) stop("File ",to_lift," does not exist!")

  if(is.null(out_file)) out_file <- "/tmp/lifted.bed"

  system(paste(liftover_bin, "-multiple -noSerial", to_lift, chain_file, out_file, "/tmp/unmapped.tmp"))

  n_bed_rows <- as.integer(system(paste("wc -l", to_lift, "| awk '{print $1}'"), intern=T))
  n_unmapped <- tryCatch(nrow(read.table("/tmp/unmapped.tmp")), error=function(e) 0)
  if(n_unmapped > 0) message(n_unmapped,"/",n_bed_rows," regions in ",to_lift," had parts which could not be liftOver'd:")
  system("cat /tmp/unmapped.tmp")

  out <- read.table(out_file, col.names=c("chr","start","end","name"))
  if(wrote_tmp_file) unlink(to_lift)
  unlink(c("/tmp/lifted.bed", "/tmp/unmapped.tmp"))
  out
}
