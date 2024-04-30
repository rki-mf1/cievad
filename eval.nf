// include modules - here, modules are single processes
include { SAMTOOLS_FAIDX } from './modules/samtools/faidx/main.nf'
include { HAPPY } from './modules/happy/main.nf'
include { SOMPY_SUMMARY } from './modules/misc/main.nf'


workflow{
    // ------------------
    // | Input channels |
    // ------------------
    ch_ref      = Channel.value("$baseDir/" + params.reference)
    ch_ref_idx  = SAMTOOLS_FAIDX(ch_ref)

    if (params.callsets_dir != "") {

        ch_callsets = Channel.fromPath(params.callsets_dir + "/" + "*.{vcf,vcf.gz}")
        ch_callsets
            .map { it -> tuple(it.toString().split('/')[-1].tokenize('_')[1].replaceFirst('.vcf', '').replaceFirst('.gz', '').toInteger(), file(it)) }
            .set {ch_callsets}
        ch_callsets.view()

        ch_truthsets = Channel.fromPath(params.outdir + "/" + "simulated_hap*.vcf")
        ch_truthsets
            .map { it -> tuple(it.toString().split('/')[-1].tokenize('_')[1].replaceFirst('hap', '').replaceFirst('.vcf', '').toInteger(), file(it)) }
            .set {ch_truthsets}
        // ch_truthsets.view()

        ch_truthsets.join(ch_callsets, by: 0)
            .set {ch_variantsets_map}
        // ch_variantsets_map.view()

    } else if (params.sample_sheet != "") { 

        ch_variantsets_map = Channel
            .fromPath(params.sample_sheet, checkIfExists: true)
            .splitCsv(header: true, sep: ",")
            .map {row -> [row["index"] as Integer, row["callset"], row["truthset"]]}
            .view()

    } else {

        exit 1, "ERROR: Either the sample_sheet or callsets_dir parameter has to be provided!\n"

    }

    // ------------------
    // | Main processes |
    // ------------------
    (ch_csv,ch_json) = HAPPY(ch_variantsets_map,ch_ref,ch_ref_idx)

    SOMPY_SUMMARY(ch_csv.collect())
}
