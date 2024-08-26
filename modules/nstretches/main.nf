process N_STRETCHES {

    // Job label
     tag "${ref}"

    // Store results
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.Nstretches.{bed,png}"

    // Engine settings
    conda 'conda-forge::python==3.12.1 conda-forge::matplotlib==3.9.2'

    // Resources
    cpus 1

    // Process I/O
    input:
    path ref

    output:
    path "${ref}.Nstretches.png",   emit: png
    path "${ref}.Nstretches.bed",   emit: bed

    // Job script
    script:
    """
    python ${projectDir}/aux/Nstretches.py \
        ${ref} \
        ${ref}.Nstretches.png \
        ${ref}.Nstretches.bed
    """

}
