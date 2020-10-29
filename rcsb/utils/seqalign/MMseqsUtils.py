##
# File:    MMseqsUtils.py
# Author:  J. Westbrook
# Date:    26-Oct-2020
#
# Updates:
#
##
"""

./mmseqs/bin/mmseqs easy-search adalimumab.fasta  db/db_pdb_entity entityResult.txt tmp
--min-seq-id 0.75
--format-output "query,target,pident,evalue,qlen,tlen,alnlen,taxid,taxname"
# ---
./mmseqs/bin/mmseqs createdb ./FASTA/pdb_seq_pr.fasta  db/db_pdb_entity
mkdir ncbi-taxdump && cd ncbi-taxdump
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xzvf taxdump.tar.gz
cd ..
./mmseqs/bin/mmseqs createtaxdb db/db_pdb_entity tmp --ncbi-tax-dump ncbi-taxdump --tax-mapping-file FASTA/entity_taxon.tdd

./mmseqs/bin/mmseqs easy-search antibody-seq.fasta db/db_pdb_entity antibodyAlign.txt tmp --min-seq-id 0.75 --format-output "query,target,pident,evalue,qlen,tlen,alnlen,taxid,taxname"

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import uuid

from rcsb.utils.io.ExecUtils import ExecUtils
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.seqalign.MiscUtils import MiscUtils

logger = logging.getLogger("__name__")


class MMseqsUtils(object):
    def __init__(self, **kwargs):
        self.__mmseqs2BinPath = kwargs.get("mmseqsBinPath", os.path.join("/usr", "local", "bin", "mmseqs"))
        self.__mU = MarshalUtil()
        self.__reportColsSimple = "query,target,pident,evalue,qlen,tlen,alnlen,taxid,taxname"
        self.__reportColsDefault = "query,target,taxid,taxname,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,raw,bits,qlen,tlen,qaln,taln,cigar"
        self.__reportCols = self.__reportColsDefault
        self.__keyMap = {
            "query": ("query", "str"),
            "target": ("target", "str"),
            "taxid": ("taxId", "int"),
            "taxname": ("taxName", "str"),
            "pident": ("sequenceIdentity", "float"),
            "alnlen": ("alignLen", "int"),
            "mismatch": ("mismatch", "int"),
            "gapopen": ("gapOpen", "int"),
            "qstart": ("queryStart", "int"),
            "qend": ("queryEnd", "int"),
            "tstart": ("targetStart", "int"),
            "tend": ("targetEnd", "int"),
            "evalue": ("eValue", "float"),
            "raw": ("rawScore", "float"),
            "bits": ("bitScore", "float"),
            "qlen": ("queryLen", "int"),
            "tlen": ("targetLen", "int"),
            "qaln": ("queryAlign", "str"),
            "taln": ("targetAlign", "str"),
            "cigar": ("cigar", "str"),
        }

    def createSearchDatabase(self, fastaPath, seqDbTopPath, seqDbName, **kwargs):
        """Create sequence search database from a FASTA file

        Args:
            fastaPath (str): input FASTA file path
            seqDbTopPath (str):  top path to search sequence database directories
            seqDbName (str): name of the sequence search database
            timeOut (int, optional): time out for the process excecution. Defaults to 3600 secs.
        Returns:
            (bool): True for success or False otherwise
        """
        ok = False
        try:
            timeOut = kwargs.get("timeOut", 3600)
            tmpDir = os.path.join(seqDbTopPath, "tmp")
            self.__mU.mkdir(seqDbTopPath)
            dbDir = os.path.join(seqDbTopPath, seqDbName)
            self.__mU.mkdir(dbDir)
            dbPath = os.path.join(dbDir, seqDbName)
            dbLogPath = os.path.join(seqDbTopPath, seqDbName + ".log")
            ok1 = self.__createSearchDatabase(fastaPath, dbPath, dbLogPath, timeOut=timeOut)
            ok2 = self.__createSearchIndex(dbPath, tmpDir, dbLogPath, timeOut=timeOut)
            ok = ok1 & ok2
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __createSearchDatabase(self, fastaPath, dbPath, outPath, timeOut=3600):
        """Create search database for the input fasta file"""
        ok = False
        try:
            exU = ExecUtils()

            ok = exU.run(self.__mmseqs2BinPath, execArgList=["createdb", fastaPath, dbPath], outPath=outPath, outAppend=True, timeOut=timeOut)
            logger.info("status %r", ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __createSearchIndex(self, dbPath, tmpDir, outPath, timeOut=3600):
        """Create search index for the input search database file"""
        ok = False
        try:
            exU = ExecUtils()
            ok = exU.run(self.__mmseqs2BinPath, execArgList=["createindex", dbPath, tmpDir], outPath=outPath, outAppend=True, timeOut=timeOut)
            logger.info("status %r", ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def createTaxonomySearchDatabase(self, taxonomyMappingPath, seqDbTopPath, seqDbName, timeOut=100):
        """Create taxonomy search database for existing search database

        Args:
            taxonomyMappingPath (str): taxonomy mapping file path
            seqDbTopPath (str):  top path to search sequence database directories
            seqDbName (str): name of the sequence search database
            timeOut (int, optional): execution process time out. Defaults to 100 secs.

        Returns:
            (bool): True for success or False otherwise
        """
        ok = False
        try:
            tmpDir = os.path.join(seqDbTopPath, "tmp")
            self.__mU.mkdir(tmpDir)
            dbDir = os.path.join(seqDbTopPath, seqDbName)
            dbPath = os.path.join(dbDir, seqDbName)
            dbLogPath = os.path.join(seqDbTopPath, seqDbName + ".log")
            #
            ncbiTaxonomyDumpPath = os.path.join(seqDbTopPath, "ncbi-taxonomy")
            #
            if not self.__mU.exists(ncbiTaxonomyDumpPath):
                ok1 = self.__getNcbiTaxonomyDatabaseDump(ncbiTaxonomyDumpPath)
                if not ok1:
                    logger.error("Fetching NCBI taxonomy database dump failing")
                    return ok1
            #
            ok = self.__createTaxonomySearchDatabase(taxonomyMappingPath, ncbiTaxonomyDumpPath, dbPath, tmpDir, dbLogPath, timeOut=timeOut)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __createTaxonomySearchDatabase(self, taxonomyMappingPath, ncbiTaxonomyDumpPath, dbPath, tmpDir, outPath, timeOut=100):
        """Create taxonomy search database for the input search database

        Args:
            taxonomyMappingPath (str): taxonomy mapping file path
            ncbiTaxonomyDumpPath (str): NCBI taxonomy database dump path
            dbPath (str): sequence search database path
            tmpDir (str): temporary directory
            outPath (str, optional): output log path. Defaults to "createTaxonomySearchDb.log".
            timeOut (int, optional): execution process time out. Defaults to 100 secs.

        Returns:
            (bool): True for success or False otherwise
        """
        ok = False
        try:
            exU = ExecUtils()
            ok = exU.run(
                self.__mmseqs2BinPath,
                execArgList=["createtaxdb", dbPath, tmpDir, "--ncbi-tax-dump", ncbiTaxonomyDumpPath, "--tax-mapping-file", taxonomyMappingPath],
                outPath=outPath,
                outAppend=True,
                timeOut=timeOut,
            )
            logger.info("status is %r", ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __getNcbiTaxonomyDatabaseDump(self, ncbiTaxonomyDumpPath):
        try:
            mU = MiscUtils()
            ok = mU.fetchNcbiTaxonomyDump(ncbiTaxonomyDumpPath)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def easySearchDatabase(self, fastaPath, seqDbTopPath, seqDbName, resultPath, **kwargs):
        """Search sequence database with the input FASTA file

        Args:
            fastaPath (str): query FASTA file path
            seqDbTopPath (str):  top path to search sequence database directories
            seqDbName (str): name of the sequence search database
            resultPath (str): search results path
            minSeqId (float): minimun sequence identity
            timeOut (int, optional): time out for the process excecution. Defaults to 100 secs.
            sensitivity (float, optional): sentivity for prefilter search (1-8) (default = 1)
            eValCutoff  (int, optional): e-Value cuttof (default= 100)
            formatMode (int, optional): 0: BLAST 1: SAM 2: BLAST+ 3: HTML (default: None)
            formatOutput (str, optional): output column selection (default: "query,target,pident,evalue,qlen,tlen,alnlen,taxid,taxname")

        Returns:
            (bool): True for success or False otherwise
        """
        ok = False
        #
        try:
            tmpDir = os.path.join(seqDbTopPath, "tmp")
            self.__mU.mkdir(tmpDir)
            dbDir = os.path.join(seqDbTopPath, seqDbName)
            dbPath = os.path.join(dbDir, seqDbName)
            dbLogPath = os.path.join(seqDbTopPath, seqDbName + ".log")
            ok = self.__easySearchDatabase(fastaPath, dbPath, resultPath, tmpDir, dbLogPath, **kwargs)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def searchDatabase(self, fastaPath, seqDbTopPath, seqDbName, resultPath, **kwargs):
        """Search sequence database with the input FASTA file

        Args:
            fastaPath (str): query FASTA file path
            seqDbTopPath (str):  top path to search sequence database directories
            seqDbName (str): name of the sequence search database
            resultPath (str): search results path
            minSeqId (float): minimun sequence identity
            timeOut (int, optional): time out for the process excecution. Defaults to 100 secs.
            sensitivity (float, optional): sentivity for prefilter search (1-8) (default = 1)
            eValCutoff  (int, optional): e-Value cuttof (default= 100)
            formatMode (int, optional): 0: BLAST 1: SAM 2: BLAST+ 3: HTML (default: None)
            formatOutput (str, optional): output column selection (default: "query,target,pident,evalue,qlen,tlen,alnlen,taxid,taxname")

        Returns:
            (bool): True for success or False otherwise
        """
        ok = False
        _ = resultPath
        try:
            qId = "query" + str(uuid.uuid1())
            ok = self.createSearchDatabase(fastaPath, seqDbTopPath, qId, **kwargs)
        except Exception as e:
            logger.exception("Failing indexing query with %s", str(e))
            return ok
        #
        ok = False
        try:
            tmpDir = os.path.join(seqDbTopPath, "tmp")
            self.__mU.mkdir(tmpDir)
            resultDirPath = os.path.join(seqDbTopPath, "results")
            self.__mU.mkdir(resultDirPath)
            resultDbPath = os.path.join(resultDirPath, qId)
            #
            targetDbPath = os.path.join(seqDbTopPath, seqDbName, seqDbName)
            queryDbPath = os.path.join(seqDbTopPath, qId, qId)
            dbLogPath = os.path.join(seqDbTopPath, seqDbName + ".log")
            ok1 = self.__searchDatabase(queryDbPath, targetDbPath, resultDbPath, tmpDir, dbLogPath, **kwargs)
            ok2 = self.__formatSearchResults(queryDbPath, targetDbPath, resultDbPath, resultPath, dbLogPath, **kwargs)
            ok = ok1 & ok2
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __searchDatabase(self, queryDbPath, targetDbPath, resultDbPath, tmpDir, outPath, **kwargs):
        """Search database for the input FASTA file"""
        ok = False
        try:
            minSeqId = kwargs.get("minSeqId", 0.95)
            timeOut = kwargs.get("timeOut", 100)
            sensitivity = kwargs.get("sensitivity", 1)
            eValCutoff = kwargs.get("eValCutoff", 100)

            exU = ExecUtils()
            ok = exU.run(
                self.__mmseqs2BinPath,
                execArgList=["search", queryDbPath, targetDbPath, resultDbPath, tmpDir, "--min-seq-id", str(minSeqId), "-a", "true", "-e", str(eValCutoff), "-s", str(sensitivity)],
                outPath=outPath,
                outAppend=True,
                timeOut=timeOut,
            )
            logger.info("status is %r", ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __easySearchDatabase(self, fastaPath, dbPath, resultPath, tmpDir, outPath, **kwargs):
        """Search database for the input FASTA file"""
        ok = False
        try:
            minSeqId = kwargs.get("minSeqId", 0.95)
            timeOut = kwargs.get("timeOut", 3600)
            sensitivity = kwargs.get("sensitivity", 1)
            eValCutoff = kwargs.get("eValCutoff", 100)
            formatMode = kwargs.get("formatMode", None)
            formatOutput = kwargs.get("formatOutput", "query,target,pident,evalue,qlen,tlen,alnlen,taxid,taxname")
            timeOut = kwargs.get("timeOut", 3600)
            if formatMode:
                fOptL = ["format_mode", str(formatMode)]
            else:
                fOptL = ["--format-output", formatOutput]
            exU = ExecUtils()
            ok = exU.run(
                self.__mmseqs2BinPath,
                execArgList=["easy-search", fastaPath, dbPath, resultPath, tmpDir, "--min-seq-id", str(minSeqId), "-s", str(sensitivity), "-a", "true", "-e", str(eValCutoff)]
                + fOptL,
                outPath=outPath,
                outAppend=True,
                timeOut=timeOut,
            )
            logger.info("status is %r", ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def mapDatabase(self, fastaPath, seqDbTopPath, seqDbName, resultPath, **kwargs):
        """Map similar sequences between input FASTA file and input target database.

        Args:
            fastaPath (str): query FASTA file path
            seqDbTopPath (str):  top path to search sequence database directories
            seqDbName (str): name of the sequence search database
            resultPath (str): search results path
            minSeqId (float, optional): minimun sequence identity (default=0.95)
            timeOut (int, optional): time out for the process excecution. Defaults to 3600 secs.
            sensitivity (float, optional): sentivity for prefilter search (1-8) (default = 1)
            eValCutoff  (int, optional): e-Value cuttof (default= 100)
            formatMode (int, optional): 0: BLAST 1: SAM 2: BLAST+ 3: HTML (default: None)
            formatOutput (str, optional): output column selection (default: "query,target,pident,evalue,qlen,tlen,alnlen,taxid,taxname")

        Returns:
            (bool): True for success or False otherwise
        """
        ok = False
        try:
            qId = "query" + str(uuid.uuid1())
            ok = self.createSearchDatabase(fastaPath, seqDbTopPath, qId, **kwargs)
        except Exception as e:
            logger.exception("Failing indexing query with %s", str(e))
            return ok
        #
        ok = False
        try:
            tmpDir = os.path.join(seqDbTopPath, "tmp")
            self.__mU.mkdir(tmpDir)
            resultDirPath = os.path.join(seqDbTopPath, "results")
            self.__mU.mkdir(resultDirPath)
            resultDbPath = os.path.join(resultDirPath, qId)
            #
            targetDbPath = os.path.join(seqDbTopPath, seqDbName, seqDbName)
            queryDbPath = os.path.join(seqDbTopPath, qId, qId)
            dbLogPath = os.path.join(seqDbTopPath, seqDbName + ".log")
            #
            ok1 = self.__mapDatabase(queryDbPath, targetDbPath, resultDbPath, tmpDir, dbLogPath, **kwargs)
            ok2 = self.__formatSearchResults(queryDbPath, targetDbPath, resultDbPath, resultPath, dbLogPath, **kwargs)
            ok = ok1 & ok2
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __mapDatabase(self, queryDbPath, targetDbPath, resultDbPath, tmpDir, outPath, **kwargs):
        """Map close matching sequences between query and target databases"""
        ok = False
        try:
            minSeqId = kwargs.get("minSeqId", 0.95)
            timeOut = kwargs.get("timeOut", 3600)
            sensitivity = kwargs.get("sensitivity", 1)
            eValCutoff = kwargs.get("eValCutoff", 100)
            #
            exU = ExecUtils()
            ok = exU.run(
                self.__mmseqs2BinPath,
                execArgList=["map", queryDbPath, targetDbPath, resultDbPath, tmpDir, "--min-seq-id", str(minSeqId), "-s", str(sensitivity), "-a", "true", "-e", str(eValCutoff)],
                outPath=outPath,
                outAppend=True,
                timeOut=timeOut,
            )
            logger.info("status is %r", ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __formatSearchResults(self, queryDbPath, targetDbPath, resultDbPath, resultPath, logPath, **kwargs):
        """Format search results"""
        ok = False
        try:
            formatMode = kwargs.get("formatMode", None)
            formatOutput = kwargs.get("formatOutput", self.__reportColsDefault)
            self.__reportCols = formatOutput
            timeOut = kwargs.get("timeOut", 3600)
            if formatMode:
                fOptL = ["--format-mode", str(formatMode)]
            else:
                fOptL = ["--format-output", formatOutput]
            #
            if resultPath.endswith(".csv"):
                fmtOut = "csv"
                tmpPath = resultPath[:-4] + ".txt"
            elif resultPath.endswith(".json"):
                fmtOut = "json"
                tmpPath = resultPath[:-5] + ".txt"
            else:
                fmtOut = None
                tmpPath = resultPath
            #
            exU = ExecUtils()
            ok = exU.run(
                self.__mmseqs2BinPath,
                execArgList=["convertalis", queryDbPath, targetDbPath, resultDbPath, tmpPath] + fOptL,
                outPath=logPath,
                outAppend=True,
                timeOut=timeOut,
            )
            if fmtOut:
                self.__parseAlignmentFile(tmpPath, resultPath, fmtOut)
            logger.info("status is %r", ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __parseAlignmentFile(self, inpFileName, outFileName, fmtOut="json"):
        ok = False
        try:
            headerL = self.__reportCols.split(",")
            rowL = self.__mU.doImport(inpFileName, fmt="list", rowFormat="list")
            rowDictL = []
            for row in rowL:
                rD = {}
                rL = row.split("\t")
                for ii, header in enumerate(headerL):
                    ky, dType = self.__keyMap[header]
                    if dType == "str":
                        rD[ky] = str(rL[ii]).strip()
                    elif dType == "int":
                        rD[ky] = int(rL[ii])
                    elif dType == "float":
                        rD[ky] = float(rL[ii])
                rowDictL.append(rD)
            if fmtOut == "csv":
                hL = [self.__keyMap[nm][0] for nm in headerL]
                ok = self.__mU.doExport(outFileName, rowDictL, fmt=fmtOut, fieldNames=hL)
            else:
                ok = self.__mU.doExport(outFileName, rowDictL, fmt=fmtOut)
            #
        except Exception as e:
            logger.exception("Failing for %r (%r) with %s", inpFileName, fmtOut, str(e))
        return ok
