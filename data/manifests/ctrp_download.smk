rule download_source:
    output:
        archive=RAW / "ctrpv2" / "{archive}",
    wildcard_constraints:
        archive="|".join(KNOWN_ARCHIVES).replace(".", r"\."),
    params:
        url=lambda wc: URL_BY_ARCHIVE[wc.archive],
    shell:
        """
        mkdir -p $(dirname {output.archive})
        curl --fail --location --retry 3 \
            --output {output.archive} \
            '{params.url}'
        """


rule download_gene_source:
    output:
        source=RAW_METADATA / "{gene_file}",
    wildcard_constraints:
        gene_file="|".join(DIRECT_GENE_FILES).replace(".", r"\."),
    params:
        url=lambda wc: GENE_SOURCE_BY_FILE[wc.gene_file]["url"],
    shell:
        """
        mkdir -p $(dirname {output.source})
        curl --fail --location --retry 3 \
            --output {output.source} \
            '{params.url}'
        """


rule download_human_gene_source:
    output:
        source=RAW_METADATA / "{gene_file}",
    wildcard_constraints:
        gene_file="|".join(HUMAN_GENE_FILES).replace(".", r"\."),
    params:
        url=lambda wc: GENE_SOURCE_BY_FILE[wc.gene_file]["url"],
    shell:
        """
        mkdir -p $(dirname {output.source})
        curl --fail --location --retry 3 '{params.url}' \
            | gzip -dc \
            | awk -F '\t' 'BEGIN {{OFS="\t"}} /^#/ || $1 == "9606"' \
            | gzip -c > {output.source}
        """
