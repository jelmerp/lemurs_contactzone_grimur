
## hzproj VCF
vcf_inds <- system("bcftools query -l results/geno/vcf/gatk/hzproj1.mac1.FS6.vcf.gz",
                   intern = TRUE)

## Allo
lookup_allo_file <- "metadata/hzlookup_allopatric.txt"
lookup_allo <- read.delim(lookup_file)
table(lookup_allo$pop)

## Sym
lookup_sym_file <- "results/manuscript/TableS1.txt"
lookup_sym <- read.delim(lookup_sym_file)
lookup_sym %>% filter(sp_radseq == "mmur") %>% count(poptype, filter_pass)
lookup_sym %>% filter(sp_radseq == "mgri") %>% count(poptype, filter_pass)

lookup_sym %>% filter(ID %in% vcf_inds) %>% count(pop)
