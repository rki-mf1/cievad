process SOMPY_SUMMARY {
    // Job label
    tag "sompy_summary"

    // Store results
    publishDir "${params.outdir}", mode: 'copy', pattern: "summary.sompy.stats.csv"

    // Engine settings
    //conda 'python3?'

    // Resources
    cpus 1

    // Process I/O
    input:
    path sompy_stats_csv

    output:
    path "summary.sompy.stats.csv"

    // Job script
    """
    cat ${sompy_stats_csv} \
        | sed 's/,type,total/idx,type,total/g' \
        | python ${projectDir}/aux/sompy_summary_of_summaries.py \
        > summary.sompy.stats.csv
    """


}