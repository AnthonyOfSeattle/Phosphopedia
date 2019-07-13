rule dummy_retrieve:
    input:
    output:
        raw_files
    shell:
        """
        touch {output}
        """
