##
# File:    MiscUtilsTests.py
# Author:  J. Westbrook
# Date:    25-Oct-2020
# Version: 0.001
#
# Updates:
##
"""
Test cases for sequence alignment miscellaneous utilities -
"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import time
import unittest

from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.seqalign.MiscUtils import MiscUtils


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class MiscUtilsTests(unittest.TestCase):
    """
    Test cases for sequence alignment miscellaneous utilities -
    """

    skipFlag = True

    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output", "CACHE")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        #
        configPath = os.path.join(self.__mockTopPath, "config", "dbload-setup-example.yml")
        configName = "site_info_remote_configuration"
        # configName = "site_info_configuration"
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName)
        #
        self.__pdbEntityTaxonPath = os.path.join(self.__workPath, "pdb_entity_taxon.tdd")
        self.__sabDabDumpPath = os.path.join(self.__dataPath, "TheraSAbDab_SeqStruc_OnlineDownload.csv")
        self.__sabDabFasta = os.path.join(self.__workPath, "fasta", "sabdab.fasta")
        #
        self.__proteinSequenceDataPath = os.path.join(self.__workPath, "pdb_protein_sequences.json")
        self.__proteinFastaPath = os.path.join(self.__workPath, "pdb_protein_sequences.fa")
        self.__unpRefSequenceDataPath = os.path.join(self.__workPath, "unp_reference_sequences.json")
        #
        self.__drugBankFastaPath = os.path.join(self.__dataPath, "drugbank_protein.fasta")
        self.__idMappingPath = os.path.join(self.__dataPath, "idmapping.dat.gz")
        self.__drugBankTargetTaxonPath = os.path.join(self.__workPath, "drugbank_taxon.tdd")
        self.__startTime = time.time()
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.debug("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testExportPDBProteinFasta(self):
        """Test case -  export protein sequence FASTA file"""
        try:
            _, fn = os.path.split(self.__proteinSequenceDataPath)
            pth = os.path.join(self.__dataPath, fn + ".gz")
            mU = MarshalUtil()
            pD = mU.doImport(pth, fmt="json")
            miscU = MiscUtils(cachePath=self.__workPath)
            miscU.exportFasta(pD, self.__proteinFastaPath, self.__pdbEntityTaxonPath)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMakeAntibodyFasta(self):
        """Test case -  create SabDab Antibody data dump to FASTA"""
        try:
            mU = MiscUtils()
            mU.convertSabDabDumpToFasta(inpPath=self.__sabDabDumpPath, outPath=self.__sabDabFasta)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    @unittest.skipIf(skipFlag, "Very long test")
    def testExportUniprotTaxonomyMapping(self):
        """Test case -  export the full uniprot taxonomy mapping"""
        try:
            mU = MiscUtils()
            outDirPath = os.path.join(self.__workPath, "uniprot_id_mapping_selected")
            taxMapFileName = "uniprot_taxonomy.pic"
            oD = mU.getUniprotXref(13, outDirPath, taxMapFileName, fmt="pickle", useCache=True)
            self.assertGreater(len(oD), 100000)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    @unittest.skipIf(skipFlag, "Requires supporting database")
    def testDrugBankTargetTaxonomy(self):
        """Test case -  export taxonomy for the DrugBank Targets """
        try:
            mU = MarshalUtil()
            dbSeqDict = mU.doImport(self.__drugBankFastaPath, fmt="fasta", commentStyle="default")
            upL = []
            for seqId in dbSeqDict:
                ff = seqId.split(" ")
                upL.append(ff[0].split("|")[1])
            #
            logger.info("Drugbank target count %d", len(upL))
            miscU = MiscUtils(cachePath=self.__workPath)
            mD = miscU.getUniprotTaxonomy(upL, self.__idMappingPath)
            #
            tL = []
            for rId, taxId in mD.items():
                lS = "%s\t%s" % (rId, taxId)
                tL.append(lS)
            ok = mU.doExport(self.__drugBankTargetTaxonPath, tL, fmt="list")
            #
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    @unittest.skipIf(skipFlag, "Requires supporting database")
    def testReferenceSequenceDataExport(self):
        """Test case -  export reference sequence data"""
        try:
            mU = MiscUtils()
            mU.exportRefSequenceTaxonomy(self.__cfgOb, self.__unpRefSequenceDataPath, fmt="json")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    @unittest.skipIf(skipFlag, "Requires supporting database")
    def testPDBProteinSequenceDataExport(self):
        """Test case -  export protein sequence data"""
        try:
            _, fn = os.path.split(self.__unpRefSequenceDataPath)
            pth = os.path.join(self.__dataPath, fn)
            mU = MarshalUtil()
            unpD = mU.doImport(pth, fmt="json")
            mU = MiscUtils()
            mU.exportProteinSequences(self.__cfgOb, unpD, self.__proteinSequenceDataPath, fmt="json")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    @unittest.skipIf(skipFlag, "Requires supporting database")
    def testTaxonomyExport(self):
        """Test case -  entity taxonomy export"""
        try:
            mU = MiscUtils()
            mU.exportEntityTaxonomy(self.__cfgOb, self.__pdbEntityTaxonPath)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteMiscTools():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(MiscUtilsTests("testTaxonomyExport"))
    suiteSelect.addTest(MiscUtilsTests("testMakeAntibodyFasta"))
    suiteSelect.addTest(MiscUtilsTests("testPDBProteinSequenceDataExport"))
    suiteSelect.addTest(MiscUtilsTests("testReferenceSequenceDataExport"))
    suiteSelect.addTest(MiscUtilsTests("testDrugBankTargetTaxonomy"))
    suiteSelect.addTest(MiscUtilsTests("testExportUniprotTaxonomyMapping"))
    suiteSelect.addTest(MiscUtilsTests("testExportPDBProteinFasta"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite = suiteMiscTools()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
