from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, PrimaryKeyConstraint, ForeignKey, Integer, Numeric, String
from sqlalchemy.orm import relationship

PhosphopediaBase = declarative_base()

##########################
#                        #
# Sample metadata tables #
#                        #
##########################

class Dataset(PhosphopediaBase):

    __tablename__ = "dataset"

    accession = Column(String(9), primary_key=True)
    title = Column(String)

    def __repr__(self):
        title_print_len = min(len(self.title), 25)
        return "<Dataset(accession='{}', title='{}')>".format(
            self.accession, self.title[:title_print_len]
        )
    
class Sample(PhosphopediaBase):

    __tablename__ = "sample"

    id = Column(Integer, primary_key=True)
    parentDataset = Column(String(9), ForeignKey("dataset.accession"))
    sampleName = Column(String)
    fileName = Column(String)
    fileSize = Column(Integer)
    fileLocation = Column(String)

    def __repr__(self):
        return "<Sample(id={}, sampleName='{}', fileName='{}')>".format(
            self.id, self.sampleName, self.fileName
        )

Dataset.samples = relationship("Sample", order_by=Sample.id)

class Parameters(PhosphopediaBase):

    __tablename__ = "parameters"

    sampleId = Column(Integer, primary_key=True)
    ms1Analyzer = Column(String(10))
    ms2Analyzer = Column(String(10))

    def __repr__(self):
        return "<Parameter(sampleId={}, ms1Analyzer='{}', ms2Analyzer='{}')>".format(
            self.sampleId, self.ms1Analyzer, self.ms2Analyzer
        )

class Error(PhosphopediaBase):

    __tablename__ = "error"

    sampleId = Column(Integer, primary_key=True)
    errorCode = Column(String(10), primary_key=True)

    def __repr__(self):
        return "<Error(sampleId={}, errorCode='{}')>".format(
            self.sampleId, self.errorCode
        )

#########################
#                       #
# Identification tables #
#                       #
#########################

class IdentificationMixin:
    id = Column(Integer, primary_key=True)
    label = Column(String)
    base_sequence = Column(String, index=True)
    sequence = Column(String)
    score = Column(Numeric(asdecimal=False))
    qvalue = Column(Numeric(asdecimal=False))
    pep = Column(Numeric(asdecimal=False))

class PSM(PhosphopediaBase, IdentificationMixin):

    __tablename__ = "psms"

    sample_name = Column(String)
    scan_number = Column(Integer)

    def __repr__(self):
        return "<PSM(sample_name='{}', scan_number={}, label='{}', sequence='{}')>".format(
            self.sample_name, self.scan_number, self.label, self.sequence
        )

class Peptide(PhosphopediaBase, IdentificationMixin):

    __tablename__ = "peptides"

    def __repr__(self):
        return "<Peptide(label='{}', sequence='{}')>".format(
            self.label, self.sequence
        )

class ModificationMixin:
    position = Column(Integer)
    residue = Column(String)
    mass = Column(Numeric(asdecimal=False))
    localization_score = Column(Numeric(asdecimal=False))
    alternative_positions = Column(String)

    def __repr__(self):
        return "<Modification(position={}, residue={}, mass={}, localization_score={})>".format(
            self.position, self.residue, self.mass, self.localization_score
        )

class PSMModification(PhosphopediaBase, ModificationMixin):

    __tablename__ = "psm_modifications"
    __table_args__ = (
        PrimaryKeyConstraint('psm_id', 'position'),
        {},
    )
    psm_id = Column(Integer, ForeignKey("psms.id"))
    psm = relationship("PSM", back_populates="modifications")

PSM.modifications = relationship("PSMModification", order_by=PSMModification.position, back_populates="psm", lazy='joined')

class PeptideModification(PhosphopediaBase, ModificationMixin):

    __tablename__ = "peptide_modifications"
    __table_args__ = (
        PrimaryKeyConstraint('peptide_id', 'position'),
        {},
    )
    peptide_id = Column(Integer, ForeignKey("peptides.id"))
    peptide = relationship("Peptide", back_populates="modifications")

Peptide.modifications = relationship("PeptideModification", order_by=PeptideModification.position, back_populates="peptide", lazy='joined')
