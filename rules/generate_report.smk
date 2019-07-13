rule calculate_stats:
    input:
    output:
        "qc/stats/version_0_stats.csv"
    shell:
        """
        touch {output}
        """
