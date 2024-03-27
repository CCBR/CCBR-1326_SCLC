
process STRIP_TAB {
    """
    Remove trailing tab character from file
    """
    tag { meta.id }
    container 'nciccbr/ccbr_ubuntu_base_20.04:v6'

    input:
    tuple val(meta), path(txt)

    output:
    tuple val(meta), path("*.fixed.*")

    script:
    output_txt = "${txt.baseName}.fixed.${txt.extension}"
    """
    sed 's/\\t\$//' ${txt} > ${output_txt}
    """
}
