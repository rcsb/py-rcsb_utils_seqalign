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
        self.__entityTaxonPath = os.path.join(self.__dataPath, "entity_taxon.tdd")
        # self.__queryFastaPath = os.path.join(self.__dataPath, "adalimumab.fasta")
        self.__queryFastaPath = os.path.join(self.__dataPath, "sabdab.fasta")

        self.__sabDabFasta = os.path.join(self.__dataPath, "sabdab.fasta")
        self.__antibodySeqDbName = "sabdab"
        self.__seqDbNamePDB = "sabdab"
        self.__seqDbTopPath = os.path.join(self.__workPath, "db")
        self.__timeOut = 3600
        self.__pdbFastaPath = os.path.join(self.__mockTopPath, "cluster_data", "mmseqs_clusters_current", "pdb_seq_pr.fasta")
        self.__pdbSeqDbName = "pdbpr"
        #
        self.__easySearchResultPath = os.path.join(self.__workPath, "resultEasySearch.txt")
        self.__searchResultPath = os.path.join(self.__workPath, "resultSearch.txt")
        self.__mapResultPath = os.path.join(self.__workPath, "resultMapSearch.json")
        self.__htmlMapResultPath = os.path.join(self.__workPath, "resultMapSearch.html")
        #
        self.__startTime = time.time()
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.debug("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testACreateDatabaseAntibody(self):
        """Test case -  create sequence search database from fasta file"""
        try:
            mmS = MMseqsUtils()
            ok = mmS.createSearchDatabase(self.__sabDabFasta, self.__seqDbTopPath, self.__antibodySeqDbName, timeOut=self.__timeOut)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testACreateDatabasePDB(self):
        """Test case -  create sequence search database from fasta file"""
        try:
            mmS = MMseqsUtils()
            ok = mmS.createSearchDatabase(self.__pdbFastaPath, self.__seqDbTopPath, self.__pdbSeqDbName, timeOut=self.__timeOut)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testBCreateTaxonomyDatabasePDB(self):
        """Test case -  create taxonomy database from a search database"""
        try:
            mmS = MMseqsUtils()
            ok = mmS.createTaxonomySearchDatabase(self.__entityTaxonPath, self.__seqDbTopPath, self.__pdbSeqDbName, timeOut=self.__timeOut)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testEasySearchDatabase(self):
        """Test case -  easy search workflow"""
        try:
            mmS = MMseqsUtils()
            ok = mmS.easySearchDatabase(
                self.__queryFastaPath,
                self.__seqDbTopPath,
                self.__pdbSeqDbName,
                self.__easySearchResultPath,
                minSeqId=0.95,
                timeOut=self.__timeOut,
            )
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testSearchDatabase(self):
        """Test case -  search workflow"""
        try:
            mmS = MMseqsUtils()
            ok = mmS.searchDatabase(
                self.__queryFastaPath,
                self.__seqDbTopPath,
                self.__pdbSeqDbName,
                self.__searchResultPath,
                minSeqId=0.95,
                timeOut=self.__timeOut,
            )
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMapDatabase(self):
        """Test case -  map similar sequences workflow"""
        try:
            mmS = MMseqsUtils()
            ok = mmS.mapDatabase(
                self.__queryFastaPath,
                self.__seqDbTopPath,
                self.__pdbSeqDbName,
                self.__mapResultPath,
                minSeqId=0.95,
                timeOut=self.__timeOut,
            )
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMapDatabaseFormatOpt(self):
        """Test case -  map similar sequences workflow with format mode output"""
        try:
            mmS = MMseqsUtils()
            ok = mmS.mapDatabase(self.__queryFastaPath, self.__seqDbTopPath, self.__pdbSeqDbName, self.__htmlMapResultPath, minSeqId=0.55, timeOut=self.__timeOut, formatMode=3)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteCreateDatabases():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(MMseqsUtilsTests("testACreateDatabaseAntibody"))
    suiteSelect.addTest(MMseqsUtilsTests("testACreateDatabasePDB"))
    suiteSelect.addTest(MMseqsUtilsTests("testBCreateTaxonomyDatabasePDB"))
    #
    return suiteSelect


def suiteSearchModes():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(MMseqsUtilsTests("testEasySearchDatabase"))
    suiteSelect.addTest(MMseqsUtilsTests("testSearchDatabase"))
    suiteSelect.addTest(MMseqsUtilsTests("testMapDatabase"))
    #
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite = suiteCreateDatabases()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
    mySuite = suiteSearchModes()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
