##
# File:    MMseqsUtilsTests.py
# Author:  J. Westbrook
# Date:    25-Oct-2020
# Version: 0.001
#
# Updates:
#  6-Dec-2020 jdw Refactor tests with all target data
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
import platform
import resource
import time
import unittest

from rcsb.utils.io.MarshalUtil import MarshalUtil
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
        self.__seqDbTopPath = os.path.join(self.__workPath, "db")
        self.__timeOut = 3600
        self.__identityCutoff = 0.99
        #
        self.__seqDataTupL = [
            (os.path.join(self.__dataPath, "pdb-protein-entity.fa.gz"), os.path.join(self.__dataPath, "pdb-protein-entity-taxon.tdd.gz"), "pdbpr"),
            (os.path.join(self.__dataPath, "drugbank-targets.fa"), os.path.join(self.__dataPath, "drugbank-targets-taxon.tdd"), "drugbank"),
            (os.path.join(self.__dataPath, "card-targets.fa"), os.path.join(self.__dataPath, "card-targets-taxon.tdd"), "card"),
            (os.path.join(self.__dataPath, "chembl-targets.fa"), os.path.join(self.__dataPath, "chembl-targets-taxon.tdd"), "chembl"),
            (os.path.join(self.__dataPath, "pharos-targets.fa"), os.path.join(self.__dataPath, "pharos-targets-taxon.tdd"), "pharos"),
            (os.path.join(self.__dataPath, "sabdab-targets.fa"), None, "sabdab"),
        ]
        # ---

        self.__startTime = time.time()
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)\n", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testACreateDatabases(self):
        """Test case -  for PDB proteins, create sequence search database from fasta file"""
        try:
            mmS = MMseqsUtils(cachePath=self.__workPath)
            for (fastaPath, taxonPath, dbName) in self.__seqDataTupL:
                ok = mmS.createSearchDatabase(fastaPath, self.__seqDbTopPath, dbName, timeOut=self.__timeOut)
                self.assertTrue(ok)
                if taxonPath:
                    ok = mmS.createTaxonomySearchDatabase(taxonPath, self.__seqDbTopPath, dbName, timeOut=self.__timeOut)
                    self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    #
    def testEasySearch(self):
        """Test case -  easy search workflow - all vs PDB """
        try:
            resultDirPath = os.path.join(self.__workPath, "easy-search-results")
            mU = MarshalUtil()
            mU.mkdir(resultDirPath)
            #
            mmS = MMseqsUtils(cachePath=self.__workPath)
            qTup = self.__seqDataTupL[0]
            for (fastaPath, _, dbName) in self.__seqDataTupL[1:]:
                resultPath = os.path.join(self.__workPath, "easy-search-results", dbName + "-results.txt")
                ok = mmS.easySearchDatabase(
                    fastaPath,
                    self.__seqDbTopPath,
                    qTup[2],
                    resultPath,
                    minSeqId=self.__identityCutoff,
                    timeOut=self.__timeOut,
                )
                self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testSearchDatabase(self):
        """Test case -  search workflow  - all vs PDB"""
        try:
            resultDirPath = os.path.join(self.__workPath, "search-results")
            mU = MarshalUtil()
            mU.mkdir(resultDirPath)
            #
            mmS = MMseqsUtils(cachePath=self.__workPath)
            qTup = self.__seqDataTupL[0]
            for (fastaPath, _, dbName) in self.__seqDataTupL[1:]:
                resultPath = os.path.join(self.__workPath, "search-results", dbName + "-results.txt")
                ok = mmS.searchDatabase(
                    fastaPath,
                    self.__seqDbTopPath,
                    qTup[2],
                    resultPath,
                    minSeqId=self.__identityCutoff,
                    timeOut=self.__timeOut,
                )
                self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMapFastaDatabase(self):
        """Test case -  map similar sequences (fasta) workflow  - all vs PDB"""
        try:
            resultDirPath = os.path.join(self.__workPath, "map-results")
            mU = MarshalUtil()
            mU.mkdir(resultDirPath)
            #
            mmS = MMseqsUtils(cachePath=self.__workPath)
            qTup = self.__seqDataTupL[0]
            for (fastaPath, _, dbName) in self.__seqDataTupL[1:]:
                resultPath = os.path.join(self.__workPath, "map-results", dbName + "-results.txt")
                ok = mmS.mapDatabaseFasta(
                    fastaPath,
                    self.__seqDbTopPath,
                    qTup[2],
                    resultPath,
                    minSeqId=self.__identityCutoff,
                    timeOut=self.__timeOut,
                )
                self.assertTrue(ok)
                mL = mmS.getMatchResults(resultPath, None, useTaxonomy=False, misMatchCutoff=10, sequenceIdentityCutoff=95.0)
                self.assertGreaterEqual(len(mL), 400)
                #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMapFastaDatabaseFormatHtml(self):
        """Test case -  map similar sequences (fasta) workflow (HTML FORMAT) - all vs PDB"""
        try:
            resultDirPath = os.path.join(self.__workPath, "map-results-html")
            mU = MarshalUtil()
            mU.mkdir(resultDirPath)
            #
            mmS = MMseqsUtils(cachePath=self.__workPath)
            qTup = self.__seqDataTupL[0]
            for (fastaPath, _, dbName) in self.__seqDataTupL[1:]:
                resultPath = os.path.join(self.__workPath, "map-results-html", dbName + "-results.txt")
                ok = mmS.mapDatabaseFasta(fastaPath, self.__seqDbTopPath, qTup[2], resultPath, minSeqId=0.55, timeOut=self.__timeOut, formatMode=3)
                self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMapDatabaseDatabase(self):
        """Test case -  map similar sequences (database) workflow  - all vs PDB"""
        try:
            resultDirPath = os.path.join(self.__workPath, "map-results-db")
            mU = MarshalUtil()
            mU.mkdir(resultDirPath)
            #
            mmS = MMseqsUtils(cachePath=self.__workPath)
            qTup = self.__seqDataTupL[0]
            for (_, taxonPath, dbName) in self.__seqDataTupL[1:]:
                resultPath = os.path.join(self.__workPath, "map-results-db", dbName + "-results.txt")
                mmS = MMseqsUtils(cachePath=self.__workPath)
                ok = mmS.mapDatabase(dbName, self.__seqDbTopPath, qTup[2], resultPath, minSeqId=0.90, timeOut=self.__timeOut)
                self.assertTrue(ok)
                mL = mmS.getMatchResults(resultPath, taxonPath)
                self.assertGreaterEqual(len(mL), 600)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteCreateDatabases():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(MMseqsUtilsTests("testACreateDatabases"))
    #
    return suiteSelect


def suiteSearchModes():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(MMseqsUtilsTests("testEasySearch"))
    suiteSelect.addTest(MMseqsUtilsTests("testSearchDatabase"))
    suiteSelect.addTest(MMseqsUtilsTests("testMapFastaDatabase"))
    suiteSelect.addTest(MMseqsUtilsTests("testMapFastaDatabaseFormatHtml"))
    suiteSelect.addTest(MMseqsUtilsTests("testMapDatabaseDatabase"))
    #
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite = suiteCreateDatabases()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
    mySuite = suiteSearchModes()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
