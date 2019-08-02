from __future__ import print_function
import re
import sys
import time
import argparse
from pyopenms import *
from numpy import isclose, nan

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class Profiler:
    """
    Consumes Spectra and Peptides to determine speed of analysis.
    Will terminate when predetermined number of spectra analyzed.
    """
    def __init__(self, n_to_profile, dest_file):
        self.n_to_profile = n_to_profile
        self.n_analyzed = 0
        self.start_time = 0
        self.elapsed = []
        self.n_peaks = []
        self.n_phospho = []
        self.n_sty = []
        self.peptide_length = []
        self.sequences = []
        self.dest_file = dest_file
    
    def reset_time(self):
        self.start_time = time.time()

    def consume_spectrum(self, spectrum):
        mz = spectrum.get_peaks()[0]
        self.n_peaks.append(mz.shape[0])

    def consume_peptide(self, peptide):
        self.n_phospho.append(peptide.getSequence().toString().count("Phospho"))
        unmodified_peptide = re.sub("\[[^A-Z]+\]", "", peptide.getSequence().toBracketString())
        self.n_sty.append(len(re.findall("[STY]", unmodified_peptide)))
        self.peptide_length.append(len(unmodified_peptide))
        self.sequences.append(unmodified_peptide)

    def consume(self, spectrum, peptide):
        if self.n_to_profile > 0:
            self.elapsed.append(time.time() - self.start_time)
            self.consume_spectrum(spectrum)
            self.consume_peptide(peptide)
            self.n_analyzed += 1
            if self.n_analyzed % self.n_to_profile == 0:
                self.report_and_quit()

    def report_and_quit(self):
        with open(self.dest_file, "w") as dest:
            header = "\t".join(["Sequence",
                                "Length",
                                "N_STY",
                                "N_Phospho",
                                "N_Peaks",
                                "Elapsed_Time"])
            dest.write(header + "\n")
            for res in zip(self.sequences,
                           self.peptide_length,
                           self.n_sty,
                           self.n_phospho,
                           self.n_peaks,
                           self.elapsed):
                line = "\t".join([str(r) for r in res])
                dest.write(line + "\n")
        quit()


class AscoreAnalyzer:
    def __init__(self, spec_file, ident_file, out_file, hit_depth=1, n_to_profile=0, profile_file_name="profile.tsv", **kwargs):
        """
        AscoreAnalyzer provides a wrapper to load Spectra from MzML files
        and peptide identifications from Comet from pepXML files. It then
        zips them together, and matches hits to spectra based on retention
        time and M/Z (guessing must occur since OpenMS does not retain scan
        information when loading identifications). This class then handles
        Ascore and allows printing TSVs.
        """
        
        # Result properties
        self.out_file = out_file
        self.results = []

        # Spectra/Hit matching properties
        self.hit_depth=hit_depth
        self.n_spectra_with_match = 0

        # QC
        self.prof = Profiler(n_to_profile, profile_file_name)

        # Load spectra
        exp = MSExperiment()
        MzMLFile().load(spec_file, exp)
        self.spectra = exp.getSpectra()
        # OpenMS does not retain scan indices, so we inject them here.
        [spec.setMetaValue("index", ind + 1) for ind, spec in enumerate(self.spectra)]
        self.spectra = [spec for spec in self.spectra if spec.getMSLevel() == 2]
        self.spectra.sort(key=lambda spec: spec.getRT())

        # Load identifications
        eprint(
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
        self.protein_records = []
        PepXMLFile().load(ident_file, self.protein_records, self.peptide_records)
        self.peptide_records.sort(key=lambda pep: pep.getRT())
        print(len(self.peptide_records))
        # AScore Init
        self.ascore = AScore()
        ascore_params = self.ascore.getParameters()
        for key, val in kwargs.items():
            ascore_params.setValue(key, val)
        self.ascore.setParameters(ascore_params)

    def check_match_count(self):
        if self.n_spectra_with_match < len(self.peptide_records):
            eprint(
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

            # Grab the relevant statistics
            # What are the spectrum RTs? Are they Close?
            # Are the spectrum MZ values close?
            record_rt = self.peptide_records[record_ind].getRT()
            spec_rt = self.spectra[spec_ind].getRT()
            rt_is_close = isclose(record_rt, spec_rt, 
                                  rtol=0., atol=.06)
            mz_is_close = isclose(self.peptide_records[record_ind].getMZ(), 
                                  self.spectra[spec_ind].getPrecursors()[0].getMZ(),
                                  rtol=10e-6, atol=.0)

            # If everything matches up, return hits
            # Otherwise, increment the spectra or records to keep up
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
        for spectrum, hit in self.generate_pairs():
            nphospho = hit.getSequence().toString().count("Phospho")
            if nphospho > 0:
                self.prof.reset_time()
                ascore_hit = self.ascore.compute(hit, spectrum)
                scores = [str(ascore_hit.getMetaValue("AScore_{}".format(i))) for i in range(1, 3 + 1)]
                self.results.append((spectrum.getMetaValue("index"), 
                                     ascore_hit.getSequence().toBracketString(), 
                                     nphospho, ";".join(scores)))
                self.prof.consume(spectrum, hit)

    def to_tsv(self):
        with open(self.out_file, "w") as dest:
            dest.write("\t".join(["Scan", "Peptide", "NPhospho", "Ascores"]))
            dest.write("\n")
            for res in self.results:
                dest.write("\t".join([str(e) for e in res]))
                dest.write("\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Maps OpenMS' implementation of Ascore accross Comet identifications.",
        prog="OpenMS-Ascore Wrapper"
    )
    parser.add_argument("--hit_depth", default=1, type=int,
                        help="Number of unique peptide hits to analyze per spectrum. "
                             "May be less if not enough hits can be found.")
    parser.add_argument("--fragment_mass_tolerance", default=.05, type=float,
                        help="Mass tolerance for matching spectra peaks with theoretical "
                             "peaks. In Da.")
    parser.add_argument("--max_peptide_length", default=50, type=int,
                        help="Maximum length peptide hit to consider.")
    parser.add_argument("--n_to_profile", default=0, type=int,
                        help="Number of peptides to profile for analysis time. "
                             "Profiling only occurs if > 0.")
    parser.add_argument("--profile_file_name", default="profile.tsv", type=str,
                        help="Where to output profile data. Only relevant if n_to_profile > 0.")
    parser.add_argument("spec_file", type=str,
                        help="MS Spectra file supplied as MZML")
    parser.add_argument("ident_file", type=str,
                        help="Comet hits supplied as pepXML")
    parser.add_argument("out_file", type=str,
                        help="Destination for ascores")
    args = vars(parser.parse_args())
    
    # Algorithm script
    run_start = time.time()
    eprint("--> OpenMS-Ascore Wrapper started on: {}".format(time.ctime()))

    eprint("--> Loading MzML and pepXML files...")
    analyzer = AscoreAnalyzer(**args)
    eprint("--> MzML and pepXML parsing complete on: {}".format(time.ctime()))
    
    eprint("--> Calculating Ascores...")
    analyzer.analyze()
    eprint("--> Ascore calculations complete on: {}".format(time.ctime()))

    eprint("--> Writing output...")
    analyzer.to_tsv()
    eprint("--> OpenMS-Ascore Wrapper completed on: {}".format(time.ctime()))
    eprint("--> Total runtime: {:.1f} seconds".format(time.time() - run_start))
