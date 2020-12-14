from pysam import TabixFile
import pandas as pd
import re
import sys


class ExploreGnomad:
    def __init__(self, gnomad_file, frequency_table):
        self.gnomad = TabixFile(gnomad_file)
        self.frequencies = pd.read_csv(frequency_table, sep="\t", header=None)
        self.frequencies.columns = ["CHR:POS", "REF", "ALT", "AF"]
        self.frequencies[["CHR", "POS"]] = self.frequencies["CHR:POS"].str.split(":", expand=True)
        self.frequencies["POS"] = self.frequencies["POS"].astype(int)

    def search_position(self, chr, pos, ref, alt):
        query_lines = self.gnomad.fetch(chr, pos-1, pos)
        for variant in query_lines:
            variant_split = variant.split("\t")
            var_ref, var_alt = variant_split[3:5]
            if ref == var_ref and alt == var_alt:
                info_line = variant_split[-1]
                match = re.search(";AF_nfe=([0-9.e+\\-]+);", info_line)
                if match:
                    return match.group(1)
        return None

    def search_all(self, output_path):
        nfe_AF = [None] * len(self.frequencies)
        for i, row in self.frequencies.iterrows():
            if i % 1000 == 0:
                print(f"{round(100*i/len(self.frequencies))} % Done")
            nfe_AF[i] = self.search_position(
                row["CHR"], row["POS"], row["REF"], row["ALT"])

        self.frequencies["nfe_AF"] = nfe_AF
        self.frequencies.to_csv(
            output_path, sep="\t", index=None,
            columns=["CHR", "POS", "REF", "ALT", "AF", "nfe_AF"]
        )


if __name__ == '__main__':
    args = sys.argv[1:]
    gnomad = args[0]
    frequency_file = args[1]
    output_path = args[2]
    explorer = ExploreGnomad(gnomad, frequency_file)
    explorer.search_all(output_path)