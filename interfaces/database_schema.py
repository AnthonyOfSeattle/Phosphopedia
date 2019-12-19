from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, PrimaryKeyConstraint, ForeignKey, Integer, Numeric, String
from sqlalchemy.orm import relationship

PhosphopediaBase = declarative_base()

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
