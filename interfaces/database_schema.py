from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, ForeignKey, Integer, Numeric, String
from sqlalchemy.orm import relationship

PhosphopediaBase = declarative_base()

###########################################
#                                         #
# Storage of PSMs generated from pipeline #
#                                         #
###########################################

class PSM(PhosphopediaBase):

    __tablename__ = "psms"

    id = Column(Integer, primary_key=True)
    sample_name = Column(String)
    scan_number = Column(Integer)
    psm_label = Column(String)
    base_sequence = Column(String, index=True)
    psm_score = Column(Numeric(asdecimal=False))
    psm_qvalue = Column(Numeric(asdecimal=False))
    psm_pep = Column(Numeric(asdecimal=False))

    def __repr__(self):
        return "<PSM(sample_name='{}', scan_number={}, psm_label='{}', base_sequence='{}')>".format(
            self.sample_name, self.scan_number, self.psm_label, self.base_sequence
        )

class Modification(PhosphopediaBase):

    __tablename__ = "modifications"

    psm_id = Column(Integer, ForeignKey("psms.id"), primary_key=True)
    position = Column(Integer, primary_key=True)
    residue = Column(String)
    mass = Column(Numeric(asdecimal=False))
    localization_score = Column(Numeric(asdecimal=False))
    alternative_positions = Column(String)
    
    psm = relationship("PSM", back_populates="modifications")

    def __repr__(self):
        return "<Modification(position={}, residue={}, mass={}, localization_score={})>".format(
            self.position, self.residue, self.mass, self.localization_score
        )

PSM.modifications = relationship("Modification", order_by=Modification.position, back_populates="psm", lazy='joined')
