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
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName)
        #
        self.__entityTaxonPath = os.path.join(self.__workPath, "entity_taxon.tdd")
        self.__sabDabDumpPath = os.path.join(self.__dataPath, "TheraSAbDab_SeqStruc_OnlineDownload.csv")
        self.__sabDabFasta = os.path.join(self.__workPath, "fasta", "sabdab.fasta")
        #
        self.__startTime = time.time()
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.debug("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    @unittest.skipIf(skipFlag, "Requires supporting database")
    def testTaxonomyExport(self):
        """Test case -  entity taxonomy export"""
        try:
            mU = MiscUtils()
            mU.exportEntityTaxonomy(self.__cfgOb, self.__entityTaxonPath)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMakeFasta(self):
        """Test case -  create SabDab Antibody data dump to FASTA"""
        try:
            mU = MiscUtils()
            mU.convertSabDabDumpToFasta(inpPath=self.__sabDabDumpPath, outPath=self.__sabDabFasta)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testFetchNcbiTaxonomyDump(self):
        """Test case -  fetch NCBI taxonomy database dump file bundle"""
        try:
            mU = MiscUtils()
            mU.fetchNcbiTaxonomyDump(os.path.join(self.__workPath, "ncbi-taxonomy"))
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteMiscTools():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(MiscUtilsTests("testTaxonomyExport"))
    suiteSelect.addTest(MiscUtilsTests("testMakeFasta"))
    suiteSelect.addTest(MiscUtilsTests("testFetchNcbiTaxonomyDump"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite = suiteMiscTools()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
