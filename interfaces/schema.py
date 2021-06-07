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
    title = Column(String(250))

    def __repr__(self):
        title_print_len = min(len(self.title), 25)
        return "<Dataset(accession='{}', title='{}')>".format(
            self.accession, self.title[:title_print_len]
        )
    
class Sample(PhosphopediaBase):

    __tablename__ = "sample"

    id = Column(Integer, primary_key=True)
    parentDataset = Column(String(9), 
                           ForeignKey("dataset.accession"))
    sampleName = Column(String(100))
    fileType = Column(String(4))
    fileName = Column(String(100))
    fileLocation = Column(String(250))

    def __repr__(self):
        return "<Sample(id={}, sampleName='{}', fileType='{}')>".format(
            self.id, self.sampleName, self.fileType
        )

Dataset.samples = relationship("Sample", order_by=Sample.id)

class Parameters(PhosphopediaBase):

    __tablename__ = "parameters"

    sampleId = Column(Integer, primary_key=True)
    ms1Analyzer = Column(String(10))
    ms2Analyzer = Column(String(10))

    def __repr__(self):
        return "<Parameters(sampleId={}, ms1Analyzer='{}', ms2Analyzer='{}')>".format(
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

class Flag(PhosphopediaBase):

    __tablename__ = "flag"

    sampleId = Column(Integer, primary_key=True)
    flagCode = Column(String(10), primary_key=True)

    def __repr__(self):
        return "<Flag(sampleId={}, flagCode='{}')>".format(
            self.sampleId, self.flagCode
        )

#########################
#                       #
# Identification tables #
#                       #
#########################

class PSM(PhosphopediaBase):

    __tablename__ = "psm"

    idPSM = Column(Integer, primary_key=True)
    idPeptide = Column(Integer)
    fdr = Column(Numeric(precision=10, 
                         scale=8,
                         asdecimal=False))
    sampleName = Column(String(75))
    scanNum = Column(Integer)
    pCharge = Column(Integer)
    pMz = Column(Numeric(precision=16,
                         scale=12,
                         asdecimal=False))
    
    def __repr__(self):
        return "<PSM(idPSM = {}, idPeptide = {}, sampleName='{}', scanNum={}".format(
            self.idPSM, self.idPeptide, self.sampleName, self.scanNum
        )

class Peptide(PhosphopediaBase):

    __tablename__ = "peptide"

    idPeptide = Column(Integer, primary_key=True)
    fdr = Column(Numeric(precision=10,   
                         scale=8,
                         asdecimal=False))
    sequence = Column(String(100))
    iRT = Column(Numeric(precision=16,
                         scale=14,
                         asdecimal=False))
    nRTExamples = Column(Integer)
    errorRT = Column(Numeric(precision=16,
                             scale=14,
                             asdecimal=False))
    z2 = Column(Integer)
    z3 = Column(Integer)
    z4 = Column(Integer)
    z5 = Column(Integer)
    z6 = Column(Integer)

    def __repr__(self):
        return "<Peptide(idPeptide = {}, sequence = '{}')>".format(
            self.idPeptide, self.sequence
        )

class PeptideSite(PhosphopediaBase):

    __tablename__ = "peptide_site"

    idPeptide = Column(Integer, primary_key=True)
    idSite  = Column(Integer, primary_key=True)
    ascore = Column(Numeric(precision=10,
                            scale=6,
                            asdecimal=False))

    def __repr__(self):
        return "<PeptideSite(idPeptide = {}, idSite = {}, ascore = {:.3f})>".format(
            self.idPeptide, self.idSite, self.ascore
        )

class Protein(PhosphopediaBase):

    __tablename__ = "protein"

    idProtein = Column(Integer, primary_key=True)
    accession = Column(String(10))
    reference = Column(String(15))
    description = Column(String(200))

    def __repr__(self):
        return "<Site(idProtein = {}, reference = '{}')>".format(
            self.idProtein, self.reference
        )

class Site(PhosphopediaBase):

    __tablename__ = "site"

    idSite = Column(Integer, primary_key=True)
    idProtein = Column(Integer)
    position = Column(Integer)
    residue = Column(String(1))
    fdr = Column(Numeric(precision=10,   
                         scale=8,
                         asdecimal=False))
    grade = Column(String(1))

    def __repr__(self):
        return "<Site(idSite = {}, idProtein = {}, site = {})>".format(
            self.idSite, self.idProtein, self.site
        )

class Pathway(PhosphopediaBase):

    __tablename__ = "pathway"

    idPathway = Column(Integer, primary_key=True)
    name = Column(String(25))

    def __repr__(self):
        return "<Pathway(idPathway = {}, name = '{}')>".format(
            self.idPathway, self.name
        )

class PathwaySite(PhosphopediaBase):

    __tablename__ = "pathway_site"

    idPathway = Column(Integer, primary_key=True)
    idSite = Column(Integer, primary_key=True)

    def __repr__(self):
        return "<PathwaySite(idPathway = {}, idSite = {})>".format(
            self.idPathway, self.idSite
        )
