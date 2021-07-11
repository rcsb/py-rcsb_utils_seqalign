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
import re
import uuid

from rcsb.utils.io.ExecUtils import ExecUtils
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.taxonomy.TaxonomyProvider import TaxonomyProvider

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
            "taxid": ("targetTaxId", "int"),
            "taxname": ("targetTaxName", "str"),
            "pident": ("sequenceIdentity", "float"),
            "alnlen": ("alignLen", "int"),
            "mismatch": ("mismatchCount", "int"),
            "gapopen": ("gapOpenCount", "int"),
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
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__taxDirPath = os.path.join(self.__cachePath, "NCBI")
        self.__taxU = None

    def createSearchDatabase(self, fastaPath, seqDbTopPath, seqDbName, **kwargs):
        """Create sequence search database from a FASTA file

        Args:
            fastaPath (str): input FASTA file path
            seqDbTopPath (str):  top path to search sequence database directories
            seqDbName (str): name of the sequence search database
            timeOut (int, optional): time out for the process execution. Defaults to 3600 secs.
        Returns:
            (bool): True for success or False otherwise
        """
        ok = False
        try:
            logger.info("Creating sequence search database for %r", seqDbName)
            timeOut = kwargs.get("timeOut", 3600)
            verbose = kwargs.get("verbose", False)
            tmpDir = os.path.join(seqDbTopPath, "tmp")
            self.__mU.mkdir(seqDbTopPath)
            dbDir = os.path.join(seqDbTopPath, seqDbName)
            self.__mU.mkdir(dbDir)
            dbPath = os.path.join(dbDir, seqDbName)
            dbLogPath = os.path.join(seqDbTopPath, seqDbName + ".log")
            ok1 = self.__createSearchDatabase(fastaPath, dbPath, dbLogPath, timeOut=timeOut)
            ok2 = self.__createSearchIndex(dbPath, tmpDir, dbLogPath, timeOut=timeOut)
            if verbose:
                logger.info("create db %r status %r", seqDbName, ok1 & ok2)
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
            logger.debug("create db %r status %r", dbPath, ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __createSearchIndex(self, dbPath, tmpDir, outPath, timeOut=3600):
        """Create search index for the input search database file"""
        ok = False
        try:
            exU = ExecUtils()
            ok = exU.run(self.__mmseqs2BinPath, execArgList=["createindex", dbPath, tmpDir], outPath=outPath, outAppend=True, timeOut=timeOut)
            logger.debug("Create index %r status %r", dbPath, ok)
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
            logger.info("Creating sequence search database for with taxonomy enabled for %r", seqDbName)
            tmpDir = os.path.join(seqDbTopPath, "tmp")
            self.__mU.mkdir(tmpDir)
            dbDir = os.path.join(seqDbTopPath, seqDbName)
            dbPath = os.path.join(dbDir, seqDbName)
            dbLogPath = os.path.join(seqDbTopPath, seqDbName + ".log")
            #
            if not self.__mU.exists(self.__taxDirPath):
                ok1 = self.__getNcbiTaxonomyDatabaseDump(self.__taxDirPath)
                if not ok1:
                    logger.error("Fetching NCBI taxonomy database dump failing")
                    return ok1
            #
            ok = self.__createTaxonomySearchDatabase(taxonomyMappingPath, dbPath, tmpDir, dbLogPath, timeOut=timeOut)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __createTaxonomySearchDatabase(self, taxonomyMappingPath, dbPath, tmpDir, outPath, timeOut=100):
        """Create taxonomy search database for the input search database

        Args:
            taxonomyMappingPath (str): taxonomy mapping file path
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
                execArgList=["createtaxdb", dbPath, tmpDir, "--ncbi-tax-dump", self.__taxDirPath, "--tax-mapping-file", taxonomyMappingPath],
                outPath=outPath,
                outAppend=True,
                timeOut=timeOut,
            )
            logger.info("create tax db %r status is %r", dbPath, ok)
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
            timeOut (int, optional): time out for the process execution. Defaults to 100 secs.
            sensitivity (float, optional): sensitivity for prefilter search (1-8) (default = 1)
            eValCutoff  (int, optional): e-Value cuttoff (default= 100)
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

    def searchDatabaseFasta(self, fastaPath, seqDbTopPath, seqDbName, resultPath, **kwargs):
        """Search sequence database with the input FASTA file

        Args:
            fastaPath (str): query FASTA file path
            seqDbTopPath (str):  top path to search sequence database directories
            seqDbName (str): name of the sequence search database
            resultPath (str): search results path
            minSeqId (float): minimun sequence identity
            timeOut (int, optional): time out for the process execution. Defaults to 100 secs.
            sensitivity (float, optional): sensitivity for prefilter search (1-8) (default = 1)
            eValCutoff  (int, optional): e-Value cuttoff (default= 100)
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

    def searchDatabase(self, queryDbName, seqDbTopPath, seqDbName, resultPath, **kwargs):
        """Search sequences from input query database with sequence in the target search database.

        Args:
            queryDbName (str): name of the query sequence database
            seqDbTopPath (str):  top path to search sequence database directories
            seqDbName (str): name of the sequence search database
            resultPath (str): search results path
            minSeqId (float): minimun sequence identity
            timeOut (int, optional): time out for the process execution. Defaults to 100 secs.
            sensitivity (float, optional): sensitivity for prefilter search (1-8) (default = 1)
            eValCutoff  (int, optional): e-Value cuttoff (default= 100)
            formatMode (int, optional): 0: BLAST 1: SAM 2: BLAST+ 3: HTML (default: None)
            formatOutput (str, optional): output column selection (default: "query,target,pident,evalue,qlen,tlen,alnlen,taxid,taxname")

        Returns:
            (bool): True for success or False otherwise
        """
        ok = False
        try:
            tmpDir = os.path.join(seqDbTopPath, "tmp")
            self.__mU.mkdir(tmpDir)
            resultDirPath = os.path.join(seqDbTopPath, "results")
            self.__mU.mkdir(resultDirPath)
            #
            # resultDbPath = os.path.join(resultDirPath, queryDbName)
            qId = "query-" + str(uuid.uuid1())
            resultDbPath = os.path.join(resultDirPath, qId)
            targetDbPath = os.path.join(seqDbTopPath, seqDbName, seqDbName)
            queryDbPath = os.path.join(seqDbTopPath, queryDbName, queryDbName)

            targetDbPath = os.path.join(seqDbTopPath, seqDbName, seqDbName)
            dbLogPath = os.path.join(seqDbTopPath, seqDbName + ".log")
            ok1 = self.__searchDatabase(queryDbPath, targetDbPath, resultDbPath, tmpDir, dbLogPath, **kwargs)
            if not ok1:
                logger.info("search %s with %s returning %r", queryDbName, seqDbName, ok1)
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
            sensitivity = kwargs.get("sensitivity", 1.0)
            eValCutoff = kwargs.get("eValCutoff", 100)
            appendMode = kwargs.get("appendMode", False)

            exU = ExecUtils()
            ok = exU.run(
                self.__mmseqs2BinPath,
                # execArgList=["search", queryDbPath, targetDbPath, resultDbPath, tmpDir, "--min-seq-id", str(minSeqId), "-a", "true", "-e", str(eValCutoff), "-s", str(sensitivity)],
                execArgList=[
                    "search",
                    queryDbPath,
                    targetDbPath,
                    resultDbPath,
                    tmpDir,
                    "--min-seq-id",
                    str(minSeqId),
                    "-a",
                    "true",
                    "-e",
                    str(eValCutoff),
                    "-s",
                    str(sensitivity),
                    # "--cov-mode",
                    # "1",
                    # "-c",
                    # "0.70",
                ],
                outPath=outPath,
                outAppend=appendMode,
                timeOut=timeOut,
            )
            logger.debug("status is %r", ok)
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
                execArgList=["easy-search", fastaPath, dbPath, resultPath, tmpDir, "--min-seq-id", str(minSeqId), "-s", str(sensitivity), "-a", "true", "-e", str(eValCutoff)] + fOptL,
                outPath=outPath,
                outAppend=True,
                timeOut=timeOut,
            )
            logger.debug("Easy-search %r status is %r", dbPath, ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def mapDatabaseFasta(self, fastaPath, seqDbTopPath, seqDbName, resultPath, **kwargs):
        """Map similar sequences between input FASTA file and input target database.

        Args:
            fastaPath (str): query FASTA file path
            seqDbTopPath (str):  top path to search sequence database directories
            seqDbName (str): name of the sequence search database
            resultPath (str): search results path
            minSeqId (float, optional): minimun sequence identity (default=0.95)
            timeOut (int, optional): time out for the process execution. Defaults to 3600 secs.
            sensitivity (float, optional): sensitivity for prefilter search (1-8) (default = 1)
            eValCutoff  (int, optional): e-Value cuttoff (default= 100)
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

    def mapDatabase(self, queryDbName, seqDbTopPath, seqDbName, resultPath, **kwargs):
        """Map similar sequences between input query and target database.

        Args:
            queryDbName (str): name of the query sequence database
            seqDbTopPath (str):  top path to search sequence database directories
            seqDbName (str): name of the sequence search database
            resultPath (str): search results path
            minSeqId (float, optional): minimun sequence identity (default=0.95)
            timeOut (int, optional): time out for the process execution. Defaults to 3600 secs.
            sensitivity (float, optional): sensitivity for prefilter search (1-8) (default = 1)
            eValCutoff  (int, optional): e-Value cuttoff (default= 100)
            formatMode (int, optional): 0: BLAST 1: SAM 2: BLAST+ 3: HTML (default: None)
            formatOutput (str, optional): output column selection (default: "query,target,pident,evalue,qlen,tlen,alnlen,taxid,taxname")

        Returns:
            (bool): True for success or False otherwise
        """
        ok = False
        try:
            tmpDir = os.path.join(seqDbTopPath, "tmp")
            self.__mU.mkdir(tmpDir)
            resultDirPath = os.path.join(seqDbTopPath, "results")
            self.__mU.mkdir(resultDirPath)
            resultDbPath = os.path.join(resultDirPath, queryDbName)
            #
            targetDbPath = os.path.join(seqDbTopPath, seqDbName, seqDbName)
            queryDbPath = os.path.join(seqDbTopPath, queryDbName, queryDbName)
            #
            dbLogPath = os.path.join(seqDbTopPath, queryDbName + ".log")
            #
            iP = os.path.join(resultDirPath, queryDbName + ".dbtype")
            self.__mU.remove(iP)
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
            logger.debug("status is %r", ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __formatSearchResults(self, queryDbPath, targetDbPath, resultDbPath, resultPath, logPath, **kwargs):
        """Format search results"""
        ok = False
        try:
            formatMode = kwargs.get("formatMode", None)
            formatOutput = kwargs.get("formatOutput", self.__reportColsDefault)
            formatOutput = formatOutput if formatOutput else self.__reportColsDefault
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
            logger.info("convertalis status is %r", ok)
            if fmtOut:
                self.__parseAlignmentFile(tmpPath, resultPath, fmtOut)

        except Exception as e:
            logger.exception("Failing for resultPath %r with %s", resultPath, str(e))
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
                ok = self.__mU.doExport(outFileName, rowDictL, fmt=fmtOut, indent=3)
            #
        except Exception as e:
            logger.exception("Failing for %r (%r) with %s", inpFileName, fmtOut, str(e))
        return ok

    def __getNcbiTaxonomyDatabaseDump(self, ncbiTaxonomyDumpPath):
        ok = False
        try:
            self.__taxU = TaxonomyProvider(taxDirPath=ncbiTaxonomyDumpPath, useCache=True)
            ok = self.__taxU.testCache()
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def getMatchResults(self, jsonSearchResultPath, queryTaxonomyMappingPath=None, useTaxonomy=False, misMatchCutoff=-1, sequenceIdentityCutoff=95.0, useBitScore=False):
        """Get matches applying input criteria on taxonomy, mismatches and sequence identity.

        Args:
            jsonSearchResultPath (str): input search results path (format=json)
            queryTaxonomyMappingPath (str, optional): input taxon mapping file (format=tdd). Defaults to None.
            useTaxonomy (bool, optional): filter results on matching taxonomy. Defaults to False.
            misMatchCutoff (int, optional): filter output on residue mismatches. Defaults to -1. (not used)
            sequenceIdentityCutoff (float, optional): minimum sequence identity. Defaults to 95.0.
            useBitScore (bool, optional): filter results comparing bit scores in the query comment. Defaults to False.

        Returns:
            dict: dictionary of lists of dictionaries containing match results
        """

        rL = self.__mU.doImport(jsonSearchResultPath, fmt="json")
        queryTaxonD = None
        if useTaxonomy and queryTaxonomyMappingPath:
            rowL = self.__mU.doImport(queryTaxonomyMappingPath, fmt="tdd", rowFormat="list")
            queryTaxonD = {row[0]: int(row[1]) for row in rowL if isinstance(row[1], int)}
        mD = self.__getMatchResults(rL, queryTaxonD, useTaxonomy=useTaxonomy, misMatchCutoff=misMatchCutoff, sequenceIdentityCutoff=sequenceIdentityCutoff, useBitScore=useBitScore)
        return mD

    def __getMatchResults(self, searchDictL, queryTaxonD, useTaxonomy=False, misMatchCutoff=-1, sequenceIdentityCutoff=95.0, useBitScore=False):
        """Get matches applying input criteria on taxonomy, mismatches and sequence identity.

        Args:
            searchDictL (list): list of candidate match dictionaries
            queryTaxonD (dict): {qSeqId: taxId, ... }
            useTaxonomy (bool, optional): filter results on matching taxonomy. Defaults to False.
            misMatchCutoff (int, optional): filter output on residue mismatches. Defaults to -1. (not used)
            sequenceIdentityCutoff (float, optional): minimum sequence identity. Defaults to 95.0.
            useBitScore (bool, optional): filter results comparing bit scores in the query comment. Defaults to False.

        Returns:
            dict: {query: [
                  {
                    "query": "drugbank_target|P54289",
                    "target": "6JP5_3|1|1|1073|1073|9986",
                    "targetTaxId": 9986,
                    "targetTaxName": "Oryctolagus cuniculus",
                    "sequenceIdentity": 97.7,
                    "alignLen": 1063,
                    "mismatchCount": 24,
                    "gapOpenCount": 0,
                    "queryStart": 9,
                    "queryEnd": 1071,
                    "targetStart": 11,
                    "targetEnd": 1073,
                    "eValue": 0.0,
                    "rawScore": 6439.0,
                    "bitScore": 2546.0,
                    "queryLen": 1103,
                    "targetLen": 1073,
                    "queryAlign": "LTLTLFQSLLIGPSSEEPFPSAVTIKSWVDKM...",
                    "targetAlign": "LTLWQAWLILIGPSSEEPFPSAVTIKSWVDKM...",
                    "cigar": "1063M"
                    },..,],...
                    }
        """
        if useTaxonomy and not self.__taxU:
            logger.info("useTaxonomy flag (%r)", useTaxonomy)
            self.__getNcbiTaxonomyDatabaseDump(self.__taxDirPath)
        mD = {}
        logger.info("Starting search result with (%r) records", len(searchDictL))
        for sD in searchDictL:
            if sD["sequenceIdentity"] < sequenceIdentityCutoff:
                continue
            if misMatchCutoff >= 0 and sD["mismatch"] > misMatchCutoff:
                continue
            # --
            query = sD["query"]
            if useBitScore:
                qL = query.split("|")
                bs = float(qL[qL.index("bitScore") - 1]) if "bitScore" in qL else None
                if "bitScore" in sD and bs and float(sD["bitScore"]) < bs:
                    # logger.debug("skipping %r < %r %r", sD["bitScore"], bs, query)
                    continue
            # --
            alR = []
            if "cigar" in sD and "queryStart" in sD and "targetStart" in sD:
                alR = self.getAlignedRegions(sD["cigar"], sD["queryStart"], sD["targetStart"])
            sD["alignedRegions"] = alR
            # --
            if useTaxonomy:
                targetTaxId = int(sD["targetTaxId"]) if "targetTaxId" in sD else None
                queryTaxId = int(queryTaxonD[query]) if queryTaxonD and query in queryTaxonD else None
                if not queryTaxId:
                    # Try extracting the query taxId from the query identifier
                    qL = query.split("|")
                    qD = {ky: val for val, ky in zip(qL, qL[1:])}
                    queryTaxId = int(qD["taxId"]) if "taxId" in qD else None
                    if targetTaxId and queryTaxId:
                        status, lcaTaxId, rank = self.__taxU.compareTaxons(queryTaxId, targetTaxId)
                        sD["taxonomyMatchStatus"] = status
                        sD["lcaTaxId"] = lcaTaxId
                        sD["lcaTaxName"] = self.__taxU.getScientificName(lcaTaxId)
                        sD["lcaRank"] = rank
                        sD["queryTaxId"] = queryTaxId
                        sD["queryRank"] = self.__taxU.getRank(queryTaxId)
                        sD["queryTaxName"] = self.__taxU.getScientificName(queryTaxId)
                        sD["targetRank"] = self.__taxU.getRank(targetTaxId)

            #
            mD.setdefault(query, []).append(sD)
        #
        logger.info("Query match count %d", len(mD))
        return mD

    def __parseCigar(self, cigarString):
        """Parse the subset of CIGAR (Compact Idiosyncratic Gapped Alignment Report) codes used by mmseqs2
           string into tuple of count an action code.

        Args:
            cigarString (str): CIGAR string

        Yields:
            tuples: (count, CIGAR code)
        """
        pattern = re.compile(r"\d+[MID]{1}")
        elementL = pattern.findall(cigarString)
        for element in elementL:
            #      count        ,    CIGAR code
            yield int(element[:-1]), element[-1]

    def getAlignedRegions(self, cigarString, queryStart, targetStart):
        """Return the list of aligned region boundaries for the input query an target sequences
           based on the alignment described by the input CIGAR (Compact Idiosyncratic Gapped Alignment Report) strings.

           Alignment cigar. Each position contains either M (match), D (deletion, gap in query), or I (Insertion, gap in target)

        Args:
            cigarString (str): CIGAR string
            queryStart (int): alignment starting index in query sequence
            targetStart (int): alignment starting index in target sequence

        Returns:
            list: tuples of (count, CIGAR code)
        """
        #
        qStart = queryStart
        tStart = targetStart
        #
        alR = []
        for count, cType in self.__parseCigar(cigarString):
            logger.debug("count %d cType %s", count, cType)
            #
            if cType in ["M", "X", "="]:
                alR.append({"queryBegin": qStart, "queryEnd": qStart + count - 1, "targetBegin": tStart, "targetEnd": tStart + count - 1})
                qStart += count
                tStart += count
            elif cType in ["D"]:
                tStart += count
            elif cType in ["I"]:
                qStart += count
            else:
                logger.error("Unrecognized CIGAR code %r", cType)
        logger.debug("Aligned regions %r", alR)
        return alR
