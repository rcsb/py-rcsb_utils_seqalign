##
# File:    MMseqsUtilsTests.py
# Author:  J. Westbrook
# Date:    25-Oct-2020
# Version: 0.001
#
# Updates:
##
"""
Test cases for MMSeqs2 sequence alignment utilities -
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
from rcsb.utils.seqalign.MMseqsUtils import MMseqsUtils


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class MMseqsUtilsTests(unittest.TestCase):
    """
    Test cases for MMSeqs2 sequence alignment utilities -
    """

    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output", "CACHE")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        #
        configPath = os.path.join(self.__mockTopPath, "config", "dbload-setup-example.yml")
        configName = "site_info_remote_configuration"
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName)
        #
        self.__seqDbTopPath = os.path.join(self.__workPath, "db")
        self.__timeOut = 3600
        #
        self.__drugBankFastaPath = os.path.join(self.__dataPath, "drugbank_protein.fasta")
        self.__drugBankTaxonPath = os.path.join(self.__dataPath, "drugbank_taxon.tdd")
        self.__drugBankSeqDbName = "drugbank"
        #
        self.__queryFastaPath = os.path.join(self.__dataPath, "sabdab.fasta")
        self.__sabDabFasta = os.path.join(self.__dataPath, "sabdab.fasta")
        self.__antibodySeqDbName = "sabdab"
        #
        # self.__pdbFastaPath = os.path.join(self.__mockTopPath, "cluster_data", "mmseqs_clusters_current", "pdb_seq_pr.fasta")
        self.__pdbFastaPath = os.path.join(self.__dataPath, "pdb_protein_sequences.fa.gz")
        self.__pdbEntityTaxonPath = os.path.join(self.__dataPath, "pdb_entity_taxon.tdd")
        self.__pdbSeqDbName = "pdbpr"
        #
        self.__antibodyEasySearchResultPath = os.path.join(self.__workPath, "antibodyFesultEasySearch.txt")
        self.__antibodySearchResultPath = os.path.join(self.__workPath, "antibodyResultSearch.txt")
        self.__antibodyMapResultPath = os.path.join(self.__workPath, "antibodyResultMapSearch.json")
        self.__antibodyHtmlMapResultPath = os.path.join(self.__workPath, "antibodyResultMapSearch.html")
        #
        self.__drugBankMapResultPath = os.path.join(self.__workPath, "drugBankResultMapSearch.json")
        #
        self.__startTime = time.time()
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.debug("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testACreateDatabaseAntibody(self):
        """Test case -  for SABDAB Antibody data, create sequence search database from fasta file"""
        try:
            mmS = MMseqsUtils(cachePath=self.__workPath)
            ok = mmS.createSearchDatabase(self.__sabDabFasta, self.__seqDbTopPath, self.__antibodySeqDbName, timeOut=self.__timeOut)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testACreateDatabasePDB(self):
        """Test case -  for PDB proteins, create sequence search database from fasta file"""
        try:
            mmS = MMseqsUtils(cachePath=self.__workPath)
            ok = mmS.createSearchDatabase(self.__pdbFastaPath, self.__seqDbTopPath, self.__pdbSeqDbName, timeOut=self.__timeOut)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testBCreateTaxonomyDatabasePDB(self):
        """Test case -  for PDB proteins, create taxonomy database from a search database"""
        try:
            mmS = MMseqsUtils(cachePath=self.__workPath)
            ok = mmS.createTaxonomySearchDatabase(self.__pdbEntityTaxonPath, self.__seqDbTopPath, self.__pdbSeqDbName, timeOut=self.__timeOut)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testACreateDatabaseDrugBank(self):
        """Test case -  for DrugBank targets, create sequence search database from fasta file"""
        try:
            mmS = MMseqsUtils(cachePath=self.__workPath)
            ok = mmS.createSearchDatabase(self.__drugBankFastaPath, self.__seqDbTopPath, self.__drugBankSeqDbName, timeOut=self.__timeOut)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testBCreateTaxonomyDatabaseDrugBank(self):
        """Test case -  for DrugBank targets create taxonomy database from a search database"""
        try:
            mmS = MMseqsUtils(cachePath=self.__workPath)
            ok = mmS.createTaxonomySearchDatabase(self.__drugBankTaxonPath, self.__seqDbTopPath, self.__drugBankSeqDbName, timeOut=self.__timeOut)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    #
    def testEasySearchDatabaseSabdabPDB(self):
        """Test case -  easy search workflow - SABDAB vs PDB """
        try:
            mmS = MMseqsUtils(cachePath=self.__workPath)
            ok = mmS.easySearchDatabase(
                self.__queryFastaPath,
                self.__seqDbTopPath,
                self.__pdbSeqDbName,
                self.__antibodyEasySearchResultPath,
                minSeqId=0.95,
                timeOut=self.__timeOut,
            )
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testSearchDatabaseSabdabPDB(self):
        """Test case -  search workflow  - SABDAB vs PDB"""
        try:
            mmS = MMseqsUtils(cachePath=self.__workPath)
            ok = mmS.searchDatabase(
                self.__queryFastaPath,
                self.__seqDbTopPath,
                self.__pdbSeqDbName,
                self.__antibodySearchResultPath,
                minSeqId=0.95,
                timeOut=self.__timeOut,
            )
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMapDatabaseSabdabPDB(self):
        """Test case -  map similar sequences workflow  - SABDAB vs PDB"""
        try:
            mmS = MMseqsUtils(cachePath=self.__workPath)
            ok = mmS.mapDatabaseFasta(
                self.__queryFastaPath,
                self.__seqDbTopPath,
                self.__pdbSeqDbName,
                self.__antibodyMapResultPath,
                minSeqId=0.95,
                timeOut=self.__timeOut,
            )
            self.assertTrue(ok)
            mL = mmS.getMatchResults(self.__antibodyMapResultPath, None, useTaxonomy=False, misMatchCutoff=10, sequenceIdentityCutoff=95.0)
            self.assertGreaterEqual(len(mL), 400)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMapDatabaseFormatOptSabdabPDB(self):
        """Test case -  map similar sequences workflow with format mode output  - SABDAB vs PDB"""
        try:
            mmS = MMseqsUtils(cachePath=self.__workPath)
            ok = mmS.mapDatabaseFasta(
                self.__queryFastaPath, self.__seqDbTopPath, self.__pdbSeqDbName, self.__antibodyHtmlMapResultPath, minSeqId=0.55, timeOut=self.__timeOut, formatMode=3
            )
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMapDatabaseDrugBankPDB(self):
        """Test case -  map similar sequences workflow with format mode output  - DrugBank vs PDB"""
        try:
            mmS = MMseqsUtils(cachePath=self.__workPath)
            ok = mmS.mapDatabase(self.__drugBankSeqDbName, self.__seqDbTopPath, self.__pdbSeqDbName, self.__drugBankMapResultPath, minSeqId=0.90, timeOut=self.__timeOut)
            self.assertTrue(ok)
            mL = mmS.getMatchResults(self.__drugBankMapResultPath, self.__drugBankTaxonPath)
            self.assertGreaterEqual(len(mL), 600)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteCreateDatabases():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(MMseqsUtilsTests("testACreateDatabaseAntibody"))
    suiteSelect.addTest(MMseqsUtilsTests("testACreateDatabasePDB"))
    suiteSelect.addTest(MMseqsUtilsTests("testBCreateTaxonomyDatabasePDB"))
    suiteSelect.addTest(MMseqsUtilsTests("testACreateDatabaseDrugBank"))
    suiteSelect.addTest(MMseqsUtilsTests("testBCreateTaxonomyDatabaseDrugBank"))
    #
    return suiteSelect


def suiteSearchModes():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(MMseqsUtilsTests("testEasySearchDatabaseSabdabPDB"))
    suiteSelect.addTest(MMseqsUtilsTests("testSearchDatabaseSabdabPDB"))
    suiteSelect.addTest(MMseqsUtilsTests("testMapDatabaseSabdabPDB"))
    suiteSelect.addTest(MMseqsUtilsTests("testMapDatabaseFormatOptSabdabPDB"))
    suiteSelect.addTest(MMseqsUtilsTests("testMapDatabaseDrugBankPDB"))
    #
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite = suiteCreateDatabases()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
    mySuite = suiteSearchModes()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
