// include modules - here, modules are single processes
//include { AMPLISIM } from './modules/amplisim/main.nf'
include { MASON_SIMULATOR } from './modules/mason/simulator/main.nf'
include { MASON_VARIATOR } from './modules/mason/variator/main.nf'
include { NANOSIM } from './modules/nanosim/main.nf'
//include { NORM_VCF } from './subworkflows/norm_vcf/main.nf'
include { SAMTOOLS_FAIDX } from './modules/samtools/faidx/main.nf'



workflow{
    // Input channels
    ch_ids      = Channel.of(1..params.n)
    ch_ref      = Channel.value("$baseDir/" + params.reference)
    ch_ref_idx  = SAMTOOLS_FAIDX(ch_ref)

    // Generate samples (haplotype consensus sequence + VCF)
    (ch_haplotypes,ch_vcf) = MASON_VARIATOR(ch_ids,ch_ref,ch_ref_idx)

    // Normalize, sort and index the VCF files
    //NORM_VCF(ch_vcf,ch_ref)

    ch_vcf
        .map { it -> tuple(it.toString().split('/')[-1].tokenize('_')[1].replaceFirst('hap', '').replaceFirst('\\.vcf', '').toInteger(), file(it)) }
        .set {ch_sample_vcf_map}

    // Generate reads
    if (params.read_type == 'ngs'){
        MASON_SIMULATOR(ch_sample_vcf_map,ch_ref,ch_ref_idx)
    }
    //else if (params.read_type == 'nano'){
    //    NANOSIM(ch_ids, ch_haplotypes)
    //}


}
