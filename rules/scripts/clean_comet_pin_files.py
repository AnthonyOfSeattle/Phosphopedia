import re
import numpy as np
import pandas as pd

def process_input(file_name):
    with open(file_name, "r") as src:
        header = next(src).split()

        psm_list = []
        for line in src:
            line = line.split()
            psm_info = line[:len(header)] 
            psm_info[-1] = "\t".join(line[(len(header) - 1):])
            psm_list.append(psm_info)

        return pd.DataFrame(psm_list, index=range( len(psm_list) ), columns=header)

def clean_n_term(pep):
    def exchange_mods(match):
        mod = float(match.group(2)[1:-1])
        if match.group(4) is not None:
            mod += float(match.group(4)[1:-1])
        return "{}{}[{}]".format(match.group(1), match.group(3), mod)

    return re.sub("(-\.)(?:n)(\[[^A-Z]+\])([A-Z])(\[[^A-Z]+\])?", exchange_mods, pep)

if __name__ == "__main__":
    original_data = process_input(snakemake.input[0])
    original_data.Peptide = pd.Series(map(clean_n_term, original_data.Peptide))

    for col in snakemake.params[0]:
        original_data = original_data.drop(col, axis=1)

    original_data.to_csv(snakemake.output[0], sep="\t", index=False)
