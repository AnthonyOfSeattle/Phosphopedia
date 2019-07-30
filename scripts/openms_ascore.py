from __future__ import print_function
import re
from pyopenms import *
from numpy import isclose, nan

class AscoreAnalyzer:
    def __init__(self, spec_file, ident_file, hit_depth=1, **kwargs):
        """
        AscoreAnalyzer provides a wrapper to load Spectra from MzML files
        and peptide identifications from Comet from pepXML files. It then
        zips them together, and matches hits to spectra based on retention
        time and M/Z (guessing must occur since OpenMS does not retain scan
        information when loading identifications). This class then handles
        Ascore and allows printing TSVs.
        """
        # Load spectra
        exp = MSExperiment()
        MzMLFile().load(spec_file, exp)
        self.spectra = exp.getSpectra()
        # OpenMS does not retain scan indices, so we inject them here.
        [spec.setMetaValue("index", ind + 1) for ind, spec in enumerate(self.spectra)]
        self.spectra = [spec for spec in self.spectra if spec.getMSLevel() == 2]
        self.spectra.sort(key=lambda spec: spec.getRT())

        # Load identifications
        print(
              """
              #################################################
              Note:
              The following warnings are from OpenMs. We do not
              specify modifications and allow the XML parser to
              infer them on its own. Please read the warnings
              and make sure nothing appears off.
              ################################################# 
              """
             )
        self.peptide_records = []
        PepXMLFile().load(ident_file, [], self.peptide_records)
        self.peptide_records.sort(key=lambda pep: pep.getRT())

        # AScore Init
        self.ascore = AScore()
        ascore_params = self.ascore.getParameters()
        for key, val in kwargs.items():
            ascore_params.setValue(key, val)
        self.ascore.setParameters(ascore_params)

        # Misc
        self.results = []
        self.hit_depth=hit_depth
        self.n_spectra_with_match = 0

    def check_match_count(self):
        if self.n_spectra_with_match < len(self.peptide_records):
            print(
                "Only {} out of {} petide records matched with spectra".format(self.n_spectra_with_match,
                                                                               len(self.peptide_records))
                 )

    def generate_hits(self, record):
        nsupplied = 0
        was_supplied = {}
        for hit in record.getHits():
            stripped_sequence = re.sub("\[[^A-Z]+\]", "", hit.getSequence().toBracketString())
            if not was_supplied.get(stripped_sequence, 0):
                was_supplied[stripped_sequence] = 1
                nsupplied += 1
                yield hit
            
            if nsupplied == self.hit_depth:
                break

    def generate_pairs(self):
        spec_ind = 0
        record_ind = 0
        while spec_ind < len(self.spectra) and record_ind < len(self.peptide_records):
            record_rt = self.peptide_records[record_ind].getRT()
            spec_rt = self.spectra[spec_ind].getRT()
            rt_is_close = isclose(record_rt, spec_rt, 
                                  rtol=0., atol=.06)
            mz_is_close = isclose(self.peptide_records[record_ind].getMZ(), 
                                  self.spectra[spec_ind].getPrecursors()[0].getMZ(),
                                  rtol=10e-6, atol=.0)

            if rt_is_close and mz_is_close:
                for hit in self.generate_hits(self.peptide_records[record_ind]):
                    yield self.spectra[spec_ind], hit
                self.n_spectra_with_match += 1
                spec_ind += 1
                record_ind += 1

            elif record_rt > spec_rt:
                spec_ind += 1

            else:
                record_ind += 1

        self.check_match_count()

    def analyze(self):
        ind = 0
        for spectrum, hit in self.generate_pairs():
            nphospho = hit.getSequence().toString().count("Phospho")
            if nphospho > 0:
                ascore_hit = self.ascore.compute(hit, spectrum)
                scores = [str(ascore_hit.getMetaValue("AScore_{}".format(i))) for i in range(1, 3 + 1)]
                self.results.insert(0, (
                                        spectrum.getMetaValue("index"), 
                                        ascore_hit.getSequence().toBracketString(), 
                                        nphospho, ";".join(scores)
                                       )
                                   )
                ind += 1
                #if ind == 1000:
                #    break
        self.results.reverse()

    def to_tsv(self, filename):
        with open(filename, "w") as dest:
            dest.write("\t".join(["Scan", "Peptide", "NPhospho", "Ascores"]))
            dest.write("\n")
            for res in self.results:
                dest.write("\t".join([str(e) for e in res]))
                dest.write("\n")

if __name__ == "__main__":
    analyzer = AscoreAnalyzer("./data/09143.mzML", "./data/09143.pep.xml",
                              hit_depth=2, fragment_mass_tolerance=.05,
                              max_peptide_length=40)
    analyzer.analyze()
    analyzer.to_tsv("./output")
