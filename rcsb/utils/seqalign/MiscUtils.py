##
# File:    MiscUtils.py
# Author:  J. Westbrook
# Date:    26-Oct-2020
#
# Updates:
#
##
"""[summary]
"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os

from rcsb.exdb.utils.ObjectExtractor import ObjectExtractor
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger("__name__")


class MiscUtils(object):
    def __init__(self, **kwargs):
        _ = kwargs

    def fetchNcbiTaxonomyDump(self, outDirPath):
        """
        mkdir ncbi-taxdump && cd ncbi-taxdump
        wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
        tar xzvf taxdump.tar.gz
        cd ..
        """
        targetUrl = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
        #
        try:
            fileU = FileUtil()
            logger.info("Fetch taxonomy data from %r in %r", targetUrl, outDirPath)
            tarPath = os.path.join(outDirPath, fileU.getFileName(targetUrl))
            ok1 = fileU.get(targetUrl, tarPath)
            ok2 = fileU.unbundleTarfile(tarPath, dirPath=outDirPath)
            ok = ok1 & ok2
            logger.info("status %r %r", ok1, ok2)
        except Exception as e:
            logger.exception("Failing for %r with %s", targetUrl, str(e))
            ok = False
        return ok

    def convertSabDabDumpToFasta(self, inpPath="TheraSAbDab_SeqStruc_OnlineDownload.csv", outPath="antibody-seq.fasta"):
        try:
            mU = MarshalUtil()
            rDL = mU.doImport(inpPath, fmt="csv", rowFormat="dict")
            logger.info("rD keys %r", list(rDL[0].keys()))
            seqObj = {}
            for rD in rDL:
                tS = "|".join(rD[ky].strip() for ky in ["Therapeutic", "Format", "CH1 Isotype", "VD LC", "Highest_Clin_Trial (Jan '20)", "Est. Status"])
                tS = "".join(tS.split())
                hSeq = rD["Heavy Sequence"] if rD["Heavy Sequence"] != "na" else None
                lSeq = rD["Light Sequence"] if rD["Light Sequence"] != "na" else None
                if hSeq:
                    seqObj[tS + "|heavy"] = {"sequence": hSeq.strip()}
                if lSeq:
                    seqObj[tS + "|light"] = {"sequence": lSeq.strip()}
            mU.doExport(outPath, seqObj, fmt="fasta")
        except Exception as e:
            logger.exception("Failing for %r with %s", inpPath, str(e))
        return

    def exportEntityTaxonomy(self, cfgOb, entityTaxonPath):
        """Extract unique entity source and host taxonomies"""
        tL = []
        try:
            obEx = ObjectExtractor(
                cfgOb,
                databaseName="pdbx_core",
                collectionName="pdbx_core_polymer_entity",
                useCache=False,
                keyAttribute="entity",
                uniqueAttributes=["rcsb_id"],
                selectionQuery=None,
                selectionList=["rcsb_id", "rcsb_entity_source_organism.ncbi_taxonomy_id", "rcsb_entity_host_organism.ncbi_taxonomy_id"],
            )
            eCount = obEx.getCount()
            logger.info("Polymer entity count is %d", eCount)
            objD = obEx.getObjects()
            sD = {}
            hD = {}
            for rId, eD in objD.items():
                try:
                    for tD in eD["rcsb_entity_source_organism"]:
                        sD.setdefault(rId, []).append(str(tD["ncbi_taxonomy_id"]))
                except Exception:
                    pass
                try:
                    for tD in eD["rcsb_entity_host_organism"]:
                        hD.setdefault(rId, []).append(str(tD["ncbi_taxonomy_id"]))
                except Exception:
                    pass
            for rId, taxIdL in sD.items():
                tS = "|".join(sorted(set(taxIdL)))
                if tS:
                    lS = "%s\t%s" % (rId, "|".join(sorted(set(taxIdL))))
                tL.append(lS)
            mU = MarshalUtil()
            mU.doExport(entityTaxonPath, tL, fmt="list")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
