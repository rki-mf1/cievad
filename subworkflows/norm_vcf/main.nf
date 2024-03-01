// include modules
include { BCFOOTLS_INDEX } from '../../modules/bcftools/index/main.nf'
include { BCFOOTLS_NORM } from '../../modules/bcftools/norm/main.nf'
include { BCFOOTLS_SORT } from '../../modules/bcftools/sort/main.nf'


workflow NORM_VCF{
    take:
    vcffiles
    reference_genome
    
    main:
    BCFOOTLS_NORM(vcffiles,reference_genome) | BCFOOTLS_SORT | BCFOOTLS_INDEX

    emit:
    ch_normed_sorted_vcffiles = BCFOOTLS_SORT.out
}