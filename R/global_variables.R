hg19_1kg_phase3_vcf_urls <- c(
  paste0("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr",1:22,".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"),
  "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"
)
hg38_1kg_30x_vcf_urls <- c(
  paste0("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr",1:22,".filtered.SNV_INDEL_SV_phased_panel.vcf.gz"),
  "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
)
hg38_gnomad_1kg_hgdp_vcf_urls <- c(
  paste0("https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr",c(1:22,"X","Y"),".vcf.bgz") # Only the Amazon links work, otherwise corrupts halfway through vcfs2Gds().
)

ucsc_human_chr_nm_map <- read.table(system.file("extdata", "ucsc_human_chr_nm_map.txt", package="refpanelutil"), header=T)
