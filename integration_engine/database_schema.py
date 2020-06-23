from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, Numeric, String, BLOB

TempBase = declarative_base()

class PSM(TempBase):

    __tablename__ = "psms"

    hash_value = Column(Integer, index=True)

    id = Column(Integer, primary_key=True)
    sample_name = Column(String)
    scan_number = Column(Integer)
    label = Column(String)
    base_sequence = Column(String)
    score = Column(Numeric(asdecimal=False))
    qvalue = Column(Numeric(asdecimal=False))
    pep = Column(Numeric(asdecimal=False))
    modifications = Column(BLOB)

class PSMPeptide(TempBase):

    __tablename__ = "psms_peptides"

    psm_id = Column(Integer, primary_key=True)
    pep_id = Column(Integer, primary_key=True)

class Peptide(TempBase):

    __tablename__ = "peptides"

    psm_id = Column(Integer, primary_key = True)
    label = Column(String)
    base_sequence = Column(String)
    sequence = Column(String)
    score = Column(Numeric(asdecimal=False))
    qvalue = Column(Numeric(asdecimal=False))
    pep = Column(Numeric(asdecimal=False))
    modifications = Column(BLOB)

class PSMProtein(TempBase):

    __tablename__ = "psms_proteins"

    prot_id = Column(Integer, index=True)
    psm_id = Column(Integer, index=True)
    id = Column(Integer, primary_key=True)

class Protein(TempBase):

    __tablename__ = "proteins"

    id = Column(Integer, primary_key=True)
    accession = Column(String)
    label = Column(String)
    coverage = Column(Numeric(asdecimal=False))

class PTM(TempBase):

    __tablename__ = "ptms"

    prot_id = Column(Integer, primary_key=True)
    position = Column(Integer, primary_key=True)
    score = Column(Numeric(asdecimal=False))
    qvalue = Column(Numeric(asdecimal=False))
    pep = Column(Numeric(asdecimal=False))
