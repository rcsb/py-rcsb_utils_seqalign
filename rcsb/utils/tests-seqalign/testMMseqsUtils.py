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

    skipFull = True

    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output", "CACHE")
        self.__dataPath = os.path.join(HERE, "test-data")
        #
        self.__seqDbTopPath = os.path.join(self.__workPath, "db")
        self.__timeOut = 3600
        self.__identityCutoff = 0.80
        self.__sensitivity = 1.0
        #
        if MMseqsUtilsTests.skipFull:
            self.__seqDataTupL = [
                (os.path.join(self.__dataPath, "pdb-protein-entity.fa.gz"), os.path.join(self.__dataPath, "pdb-protein-entity-taxon.tdd"), "pdbpr", 0),
                (os.path.join(self.__dataPath, "drugbank-targets.fa"), os.path.join(self.__dataPath, "drugbank-targets-taxon.tdd"), "drugbank", 2000),
                # (os.path.join(self.__dataPath, "imgt-reference.fa"), os.path.join(self.__dataPath, "imgt-reference-taxon.tdd"), "imgt", 200),
            ]
        else:
            self.__seqDataTupL = [
                (os.path.join(self.__dataPath, "pdb-protein-entity.fa.gz"), os.path.join(self.__dataPath, "pdb-protein-entity-taxon.tdd"), "pdbpr", 0),
                (os.path.join(self.__dataPath, "drugbank-targets.fa"), os.path.join(self.__dataPath, "drugbank-targets-taxon.tdd"), "drugbank", 2000),
                (os.path.join(self.__dataPath, "card-targets.fa"), os.path.join(self.__dataPath, "card-targets-taxon.tdd"), "card", 500),
                (os.path.join(self.__dataPath, "chembl-targets.fa"), os.path.join(self.__dataPath, "chembl-targets-taxon.tdd"), "chembl", 2000),
                (os.path.join(self.__dataPath, "pharos-targets.fa"), os.path.join(self.__dataPath, "pharos-targets-taxon.tdd"), "pharos", 2100),
                (os.path.join(self.__dataPath, "sabdab-targets.fa"), None, "sabdab", 150),
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
            for (fastaPath, taxonPath, dbName, _) in self.__seqDataTupL:
                ok = mmS.createSearchDatabase(fastaPath, self.__seqDbTopPath, dbName, timeOut=self.__timeOut, verbose=True)
                self.assertTrue(ok)
                if taxonPath:
                    ok = mmS.createTaxonomySearchDatabase(taxonPath, self.__seqDbTopPath, dbName, timeOut=self.__timeOut)
                    self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    #
    def testEasySearch(self):
        """Test case -  easy search workflow - all vs PDB"""
        try:
            resultDirPath = os.path.join(self.__workPath, "easy-search-results")
            mU = MarshalUtil()
            mU.mkdir(resultDirPath)
            #
            mmS = MMseqsUtils(cachePath=self.__workPath)
            qTup = self.__seqDataTupL[0]
            for (fastaPath, _, dbName, _) in self.__seqDataTupL[1:]:
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
            for (fastaPath, _, dbName, _) in self.__seqDataTupL[1:]:
                resultPath = os.path.join(self.__workPath, "search-results", dbName + "-results.json")
                ok = mmS.searchDatabaseFasta(
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
            for (fastaPath, _, dbName, minMatch) in self.__seqDataTupL[1:]:
                resultPath = os.path.join(self.__workPath, "map-results", dbName + "-results.json")
                ok = mmS.mapDatabaseFasta(
                    fastaPath,
                    self.__seqDbTopPath,
                    qTup[2],
                    resultPath,
                    minSeqId=self.__identityCutoff,
                    timeOut=self.__timeOut,
                )
                self.assertTrue(ok)
                mL = mmS.getMatchResults(resultPath, None, useTaxonomy=False, misMatchCutoff=-1, sequenceIdentityCutoff=self.__identityCutoff)
                logger.info("Search result %r (%d)", dbName, len(mL))
                self.assertGreaterEqual(len(mL), minMatch)
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
            for (fastaPath, _, dbName, _) in self.__seqDataTupL[1:]:
                resultPath = os.path.join(self.__workPath, "map-results-html", dbName + "-results.html")
                ok = mmS.mapDatabaseFasta(fastaPath, self.__seqDbTopPath, qTup[2], resultPath, minSeqId=self.__identityCutoff, timeOut=self.__timeOut, formatMode=3)
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
            for (_, taxonPath, dbName, minMatch) in self.__seqDataTupL[1:]:
                resultPath = os.path.join(self.__workPath, "map-results-db", dbName + "-results.json")
                mmS = MMseqsUtils(cachePath=self.__workPath)
                ok = mmS.mapDatabase(dbName, self.__seqDbTopPath, qTup[2], resultPath, minSeqId=self.__identityCutoff, sensitivity=self.__sensitivity, timeOut=self.__timeOut)
                self.assertTrue(ok)
                #
                if taxonPath:
                    mL = mmS.getMatchResults(resultPath, taxonPath, useTaxonomy=True, misMatchCutoff=-1, sequenceIdentityCutoff=self.__identityCutoff)
                else:
                    mL = mmS.getMatchResults(resultPath, None, useTaxonomy=False, misMatchCutoff=-1, sequenceIdentityCutoff=self.__identityCutoff)

                logger.info("Search result %r (%d)", dbName, len(mL))
                self.assertGreaterEqual(len(mL), minMatch)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testSearchDatabaseDatabase(self):
        """Test case -  map similar sequences (database) workflow  - all vs PDB"""
        try:
            resultDirPath = os.path.join(self.__workPath, "search-results-db")
            mU = MarshalUtil()
            mU.mkdir(resultDirPath)
            #
            mmS = MMseqsUtils(cachePath=self.__workPath)
            qTup = self.__seqDataTupL[0]
            for (_, taxonPath, dbName, minMatch) in self.__seqDataTupL[1:]:
                resultPath = os.path.join(self.__workPath, "search-results-db", dbName + "-results.json")
                mmS = MMseqsUtils(cachePath=self.__workPath)
                ok = mmS.searchDatabase(dbName, self.__seqDbTopPath, qTup[2], resultPath, minSeqId=self.__identityCutoff, sensitivity=self.__sensitivity, timeOut=self.__timeOut)
                self.assertTrue(ok)
                #
                if taxonPath:
                    mD = mmS.getMatchResults(resultPath, taxonPath, useTaxonomy=True, misMatchCutoff=-1, sequenceIdentityCutoff=self.__identityCutoff)
                else:
                    mD = mmS.getMatchResults(resultPath, None, useTaxonomy=False, misMatchCutoff=-1, sequenceIdentityCutoff=self.__identityCutoff)

                logger.info("Search result %r (%d)", dbName, len(mD))
                self.assertGreaterEqual(len(mD), minMatch)
                resultPath = os.path.join(self.__workPath, "search-results-db", dbName + "-filtered-results.json")
                ok = mU.doExport(resultPath, mD, fmt="json", indent=3)
                self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testAlignedRegions(self):
        """Test case:  reconstruct aligned regions from cigar string."""
        sD = {
            "queryStart": 1,
            "queryEnd": 842,
            "targetStart": 1,
            "targetEnd": 825,
            "cigar": "8M4I19M1I354M4I2M2D194M10I246M",
        }
        mmS = MMseqsUtils(cachePath=self.__workPath)
        alR = mmS.getAlignedRegions(sD["cigar"], sD["queryStart"], sD["targetStart"])
        self.assertEqual(alR[-1]["queryEnd"], sD["queryEnd"])
        self.assertEqual(alR[-1]["targetEnd"], sD["targetEnd"])
        logger.debug("Aligned regions %r", alR)


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
