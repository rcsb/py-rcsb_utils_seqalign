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
from rcsb.utils.io.IoUtil import IoUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.taxonomy.TaxonomyProvider import TaxonomyProvider

logger = logging.getLogger("__name__")


def getRangeOverlap(entityBeg, entityEnd, refBeg, refEnd):
    r1 = range(entityBeg, entityEnd)
    r2 = range(refBeg, refEnd)
    if r1.start == r1.stop or r2.start == r2.stop:
        return set()
    if not ((r1.start < r2.stop and r1.stop > r2.start) or (r1.stop > r2.start and r2.stop > r1.start)):
        return set()
    return set(range(max(r1.start, r2.start), min(r1.stop, r2.stop) + 1))


def firstCommonElement(l1, l2):
    try:
        for ii, v in enumerate(l1, 1):
            if v in l2:
                return v, ii
    except Exception:
        pass
    return None, None


class MiscUtils(object):
    def __init__(self, **kwargs):
        self.__cachePath = kwargs.get("cachePath", ".")
        _ = kwargs

    def getUniprotTaxonomy(self, unpIdList, idMapPath="idmapping.dat.gz"):
        """Fetch the taxonomy assignments for the input UniProt identifier list

        ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz

        ~4300 seconds ~1.2 hrs to process a single pass. macbook pro
        """
        oD = {}
        try:
            targetUrl = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz"
            if not os.access(idMapPath):
                outDirPath = os.path.join(self.__cachePath, "uniprot-idmapping")
                fileU = FileUtil()
                logger.info("Fetch idmapping data from %r in %r", targetUrl, outDirPath)
                idPath = os.path.join(outDirPath, fileU.getFileName(targetUrl))
                ok = fileU.get(targetUrl, idPath)
                if not ok:
                    logger.error("Failed to downlowd %r", targetUrl)
                    return oD
            # ---
            ioU = IoUtil()
            unpIdD = {k: True for k in unpIdList}
            iCount = 0
            iMatch = 0
            for row in ioU.deserializeCsvIter(idMapPath, delimiter="\t", rowFormat="list", encodingErrors="ignore"):
                if row[0] in unpIdD and row[1] == "NCBI_TaxID":
                    oD[row[0]] = int(row[2])
                    iMatch += 1
                    if iMatch % 500 == 0:
                        logger.info("Matched %d", iMatch)
                #
                iCount += 1
                if iCount % 50000000 == 0:
                    logger.info("Processing %d", iCount)
                #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return oD

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

    def exportRefSequenceTaxonomy(self, cfgOb, filePath, fmt="json"):
        """Export reference protein sequence taxonomy data """
        ok = False
        uD = None
        try:
            obEx = ObjectExtractor(
                cfgOb,
                databaseName="uniprot_exdb",
                collectionName="reference_entry",
                useCache=False,
                keyAttribute="uniprot",
                uniqueAttributes=["rcsb_id"],
                selectionQuery={},
                selectionList=[
                    "source_scientific",
                    "taxonomy_id",
                    "rcsb_id",
                    "gene",
                    "names",
                    "sequence",
                ],
            )
            #
            eCount = obEx.getCount()
            logger.info("Reference entry count is %d", eCount)
            objD = obEx.getObjects()
            rD = {}
            for rId, uD in objD.items():
                taxId = uD["taxonomy_id"]
                sn = uD["source_scientific"]
                sequence = uD["sequence"]
                if "gene" in uD:
                    for tD in uD["gene"]:
                        if tD["type"] == "primary":
                            gn = tD["name"]
                            break
                for tD in uD["names"]:
                    if tD["nameType"] == "recommendedName":
                        pn = tD["name"]
                        break
                rD[rId] = {"accession": rId, "taxId": taxId, "scientific_name": sn, "gene": gn, "name": pn, "sequence": sequence}
            mU = MarshalUtil()
            ok = mU.doExport(filePath, rD, fmt=fmt, indent=3)
        except Exception as e:
            logger.exception("Failing uD %r with %s", uD, str(e))
        #
        return ok

    def exportEntityTaxonomy(self, cfgOb, entityTaxonPath):
        """Export unique entity source and host taxonomies"""
        tL = []
        ok = False
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
            ok = mU.doExport(entityTaxonPath, tL, fmt="list")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def exportProteinSequences(self, cfgOb, unpD, filePath, fmt="json"):
        """Export protein sequence and taxonomy data """
        ok = False
        missingSrcD = {}
        try:
            obEx = ObjectExtractor(
                cfgOb,
                databaseName="pdbx_core",
                collectionName="pdbx_core_polymer_entity",
                useCache=False,
                keyAttribute="entity",
                uniqueAttributes=["rcsb_id"],
                selectionQuery={"entity_poly.rcsb_entity_polymer_type": "Protein"},
                selectionList=[
                    "rcsb_id",
                    "rcsb_entity_source_organism",
                    "rcsb_polymer_entity.rcsb_source_part_count",
                    "rcsb_polymer_entity.rcsb_source_taxonomy_count",
                    "rcsb_polymer_entity.src_method",
                    "entity_poly",
                    "rcsb_polymer_entity_align",
                ],
            )
            #
            eCount = obEx.getCount()
            logger.info("Polymer entity count is %d", eCount)
            objD = obEx.getObjects()
            rD = {}
            for rId, eD in objD.items():

                try:
                    pD = eD["entity_poly"]
                    seqS = pD["pdbx_seq_one_letter_code_can"]
                    seqLen = len(seqS)
                except Exception:
                    logger.warning("%s no one-letter-code sequence", rId)
                #
                if seqLen < 20:
                    continue
                #
                srcMethod = None
                try:
                    pD = eD["rcsb_polymer_entity"]
                    srcMethod = pD["src_method"]
                except Exception:
                    pass
                #
                if "rcsb_entity_source_organism" not in eD:
                    logger.debug("%s No source information (%r) skipping (seqLen %d)", rId, srcMethod, seqLen)
                    continue
                try:
                    sL = []
                    for tD in eD["rcsb_entity_source_organism"]:
                        srcName = tD["scientific_name"] if "scientific_name" in tD else None
                        if "beg_seq_num" in tD and "end_seq_num" in tD:
                            begSeqNum = tD["beg_seq_num"]
                            endSeqNum = tD["end_seq_num"] if tD["end_seq_num"] <= seqLen else seqLen
                        else:
                            begSeqNum = 1
                            endSeqNum = seqLen
                        srcId = tD["pdbx_src_id"]
                        srcType = tD["source_type"]
                        taxId = tD["ncbi_taxonomy_id"] if "ncbi_taxonomy_id" in tD else -1
                        if srcName and taxId == -1:
                            missingSrcD.setdefault(srcName, []).append(rId)
                        orgName = tD["ncbi_scientific_name"] if "ncbi_scientific_name" in tD else ""
                        sL.append({"srcId": srcId, "taxId": taxId, "orgName": orgName, "entitySeqBeg": begSeqNum, "entitySeqEnd": endSeqNum})
                    if len(sL) == 1:
                        sL[0]["entitySeqBeg"] = 1
                        sL[0]["entitySeqEnd"] = seqLen

                except Exception as e:
                    logger.exception("Failing for (%r) tD %r with %s", rId, tD, str(e))
                #
                try:
                    pD = eD["rcsb_polymer_entity"]
                    partCount = pD["rcsb_source_part_count"]
                except Exception:
                    logger.warning("%s no source part count", rId)
                    partCount = 1
                try:
                    pD = eD["rcsb_polymer_entity"]
                    taxCount = pD["rcsb_source_taxonomy_count"]
                except Exception:
                    if srcType == "synthetic":
                        taxCount = 0
                    else:
                        logger.warning("%s (srcName %r) no source taxonomy count type %r", rId, srcName, srcType)
                        if srcName:
                            taxCount = 1
                        else:
                            taxCount = 0
                #
                uDL = []
                try:
                    for tD in eD["rcsb_polymer_entity_align"]:
                        uD = {}
                        if tD["reference_database_name"] in ["UniProt", "GenBank", "PIR", "EMBL", "NORINE", "PRF"]:
                            uD["refDbId"] = tD["reference_database_accession"]
                            uD["refDbName"] = tD["reference_database_name"]
                            uD["provSource"] = tD["provenance_source"]
                            if tD["reference_database_accession"] in unpD:
                                uD.update(unpD[tD["reference_database_accession"]])
                            aL = []
                            for qD in tD["aligned_regions"]:
                                if qD["entity_beg_seq_id"] + qD["length"] - 1 > seqLen:
                                    qD["length"] = seqLen - qD["entity_beg_seq_id"] + 1
                                srcId = self.__getSourcePart(rId, sL, qD["entity_beg_seq_id"], qD["length"])

                                aL.append({"srcId": srcId, "entitySeqBeg": qD["entity_beg_seq_id"], "refSeqBeg": qD["ref_beg_seq_id"], "length": qD["length"]})
                            uD["alignList"] = aL
                            uDL.append(uD)
                        else:
                            logger.info("%s reference database %s", rId, tD["reference_database_name"])

                except Exception:
                    pass
                rD[rId] = {"alignmentL": uDL, "sourceOrgL": sL, "partCount": partCount, "taxCount": taxCount, "sequence": seqS, "seqLen": seqLen}
            # ----
            mU = MarshalUtil()
            ok = mU.doExport(filePath, rD, fmt=fmt, indent=3)
            pth, _ = os.path.split(filePath)
            mU = MarshalUtil()
            ok = mU.doExport(os.path.join(pth, "missingSrcNames.json"), missingSrcD, fmt="json")
        except Exception as e:
            logger.exception("Failing %r (%r) with %s", filePath, fmt, str(e))
        return ok

    def __getSourcePart(self, entityId, sourceOrgL, entityBeg, seqLen):
        """Return the source part containing the input entity range -

        Args:
            sourceOrgL (list): list of source dictionaries
            entityBeg (int):  begining entity sequence position (matched region)
            seqLen (int):  length sequence range (matched region)

        Returns:
            (int): corresponding source part id or None
        """
        entityEnd = entityBeg + seqLen - 1
        for sD in sourceOrgL:
            srcId = sD["srcId"]
            if sD["entitySeqBeg"] <= entityBeg and sD["entitySeqEnd"] >= entityEnd:
                return srcId
        #

        # logger.error("%r (%d) Inconsistent range for beg %r end %r sourceOrgL %r", entityId, len(sourceOrgL), entityBeg, entityEnd, sourceOrgL)
        if len(sourceOrgL) == 1:
            logger.error("%r (%d) Inconsistent range for beg %r end %r sourceOrgL %r", entityId, len(sourceOrgL), entityBeg, entityEnd, sourceOrgL)
            return 1
        else:
            ovTupL = []
            for sD in sourceOrgL:
                srcId = sD["srcId"]
                logger.debug("%r %r beg %r end %r beg %r end %r", entityId, srcId, sD["entitySeqBeg"], sD["entitySeqEnd"], entityBeg, entityEnd)
                oVS = getRangeOverlap(sD["entitySeqBeg"], sD["entitySeqEnd"], entityBeg, entityEnd)
                ovTupL.append((srcId, len(oVS)))
            rL = sorted(ovTupL, key=lambda x: x[1], reverse=True)
            logger.debug("ovTupL %r", rL)
            #
            return rL[0][0]

    def exportFasta(self, proteinSeqD, fastaPath, taxonPath, fmt="fasta"):
        """[summary]

        Args:
            proteinSeqD (dict): protein sequence and taxonomy data dictionary
            filePath (str): FASTA output file path

        Returns:
            bool: True for success or False otherwise

        Example:
            "5H7D_1": {
                    "alignmentL": [
                        {
                            "refDbId": "P42588",
                            "refDbName": "UniProt",
                            "provSource": "PDB",
                            "accession": "P42588",
                            "taxId": 83333,
                            "scientific_name": "Escherichia coli (strain K12)",
                            "gene": "patA",
                            "name": "PATase",
                            "alignList": [
                            {
                                "srcId": "1",
                                "entitySeqBeg": 5,
                                "refSeqBeg": 7,
                                "length": 447
                            }
                            ]
                        },
                        {
                            "refDbId": "P38507",
                            "refDbName": "UniProt",
                            "provSource": "PDB",
                            "accession": "P38507",
                            "taxId": 1280,
                            "scientific_name": "Staphylococcus aureus",
                            "gene": "spa",
                            "name": "IgG-binding protein A",
                            "alignList": [
                            {
                                "srcId": "2",
                                "entitySeqBeg": 452,
                                entitySeqBeg"220,
                                "length": 48
                            }
                            ]
                        }
                    ],
                    "sourceOrgL": [
                        {
                            "srcId": "1",
                            "taxId": 83333,
                            "orgName": "Escherichia coli K-12",
                            "entitySeqBeg": 1,
                            "entitySeqEnd": 451
                        },
                        {
                            "srcId": "2",
                            "taxId": 1280,
                            "orgName": "Staphylococcus aureus",
                            "entitySeqBeg": 452,
                            "entitySeqEnd": 499
                        }
                    ],
                    "partCount": 2,
                    "taxCount": 2,
                    "sequence": "GSHMSASALACSAHALNLIEKRTLDHEEMKALNREVIEYFKEHVNPGF...",
                    "seqLen": 499
                },
        >1ABC_#|prt|<taxid>|beg|end|refdb|refId|refTaxId|refbeg|refend|ref_gn|ref_nm
        """
        ok = False
        # tU = TaxonomyProvider(taxDirPath=self.__cachePath)
        # taxId = 9606
        # tL = tU.getLineage(taxId)
        # logger.info("Lineage for %d (%d): %r", taxId, len(tL), tL)
        try:
            taxonL = []
            seqDict = {}
            for eId, eD in proteinSeqD.items():
                # partCount = eD["partCount"]
                # taxCount = eD["taxCount"]
                # alignL = eD["alignmentL"] if "alignmentL" in eD else []
                seq = eD["sequence"]
                for sD in eD["sourceOrgL"]:
                    srcId = sD["srcId"]
                    taxId = sD["taxId"]
                    seqBeg = int(sD["entitySeqBeg"])
                    seqEnd = int(sD["entitySeqEnd"])
                    seqLen = 1 + (seqEnd - seqBeg)
                    # orgName = sD["orgName"]
                    seqIdL = [eId, srcId, seqBeg, seqEnd, seqLen, taxId]
                    seqId = "|".join([str(t) for t in seqIdL])
                    seqDict[seqId] = {"sequence": seq[seqBeg - 1 : seqEnd]}
                    taxonL.append("%s\t%s" % (seqId, taxId))
                # ----
            mU = MarshalUtil()
            ok = mU.doExport(fastaPath, seqDict, fmt=fmt)
            ok = mU.doExport(taxonPath, taxonL, fmt="list")
        except Exception as e:
            logger.exception("Failing %r with %s", fastaPath, str(e))
        return ok

    def __exportFasta(self, proteinSeqD, filePath, fmt="fasta"):
        """[summary]

        Args:
            proteinSeqD (dict): protein sequence and taxonomy data dictionary
            filePath (str): FASTA output file path

        Returns:
            bool: True for success or False otherwise

        Example:
            "5H7D_1": {
                    "alignmentL": [
                        {
                            "refDbId": "P42588",
                            "refDbName": "UniProt",
                            "provSource": "PDB",
                            "accession": "P42588",
                            "taxId": 83333,
                            "scientific_name": "Escherichia coli (strain K12)",
                            "gene": "patA",
                            "name": "PATase",
                            "alignList": [
                            {
                                "srcId": "1",
                                "entitySeqBeg": 5,
                                "refSeqBeg": 7,
                                "length": 447
                            }
                            ]
                        },
                        {
                            "refDbId": "P38507",
                            "refDbName": "UniProt",
                            "provSource": "PDB",
                            "accession": "P38507",
                            "taxId": 1280,
                            "scientific_name": "Staphylococcus aureus",
                            "gene": "spa",
                            "name": "IgG-binding protein A",
                            "alignList": [
                            {
                                "srcId": "2",
                                "entitySeqBeg": 452,
                                entitySeqBeg"220,
                                "length": 48
                            }
                            ]
                        }
                    ],
                    "sourceOrgL": [
                        {
                            "srcId": "1",
                            "taxId": 83333,
                            "orgName": "Escherichia coli K-12",
                            "entitySeqBeg": 1,
                            "entitySeqEnd": 451
                        },
                        {
                            "srcId": "2",
                            "taxId": 1280,
                            "orgName": "Staphylococcus aureus",
                            "entitySeqBeg": 452,
                            "entitySeqEnd": 499
                        }
                    ],
                    "partCount": 2,
                    "taxCount": 2,
                    "sequence": "GSHMSASALACSAHALNLIEKRTLDHEEMKALNREVIEYFKEHVNPGF...",
                    "seqLen": 499
                },
        >1ABC_#|prt|<taxid>|beg|end|refdb|refId|refTaxId|refbeg|refend|ref_gn|ref_nm
        """
        ok = False
        taxId = 9606
        tU = TaxonomyProvider(taxDirPath=self.__cachePath)
        tL = tU.getLineage(taxId)
        logger.info("Lineage for %d (%d): %r", taxId, len(tL), tL)
        try:
            seqDict = {}
            for eId, eD in proteinSeqD.items():
                partCount = eD["partCount"]
                taxCount = eD["taxCount"]
                alignL = eD["alignmentL"] if "alignmentL" in eD else []
                # ----
                if partCount == 1 and taxCount == 1 and len(alignL) == 1:
                    tD = eD["sourceOrgL"][0]
                    seq = eD["sequence"]
                    srcId = tD["srcId"]
                    taxId = tD["taxId"]
                    orgName = tD["orgName"]
                    seqBeg = tD["entitySeqBeg"]
                    seqEnd = tD["entitySeqEnd"]
                    seqLen = 1 + (seqEnd - seqBeg)
                    #
                    aD = eD["alignmentL"][0]
                    refDbId = aD["refDbId"]
                    refDbName = aD["refDbName"]
                    provSource = aD["provSource"]
                    refTaxId = aD["taxId"] if "taxId" in aD else -1
                    refOrgName = aD["scientific_name"] if "scientific_name" in aD else ""
                    refGeneName = aD["gene"] if "gene" in aD else ""
                    refProteinName = aD["name"] if "name" in aD else ""
                    #
                    #  >1ABC_#|prt|<taxid>|beg|end|refdb|refId|refTaxId|refbeg|refend|ref_gn|ref_nm

                    seqIdL = [eId, srcId, taxId, orgName, seqBeg, seqEnd, seqLen, refDbName, refDbId, refTaxId, refOrgName, refGeneName, refProteinName]
                    seqIdL = [str(v) for v in seqIdL]
                    seqId = "|".join(seqIdL)
                    seqDict[seqId] = {"sequence": seq}
                    if taxId != refTaxId and refTaxId != -1:
                        tL = tU.getLineage(refTaxId)
                        rtL = tU.getLineage(taxId)
                        if not ((taxId in tL) or (refTaxId in rtL)):
                            fcTaxId, nSteps = firstCommonElement(list(reversed(tL)), list(reversed(rtL)))
                            if fcTaxId:
                                sn = tU.getScientificName(fcTaxId)
                                logger.error("taxonomy mismatch common ancestor %s (%d in %d) (%r) for %r", sn, fcTaxId, nSteps, provSource, seqId)
                            else:
                                logger.error("taxonomy mismatch (%r) for %r", provSource, seqId)

                elif partCount == 1 and taxCount == 1 and len(alignL) > 1:
                    logger.info("%s Multiple reference sequences (%d) for with partCount %d", eId, len(alignL), partCount)

                elif partCount == 1 and taxCount == 1 and len(alignL) < 1:
                    # No reference sequence assignments
                    #
                    logger.debug("%s No reference sequences (%d) with partCount %d taxCount %d", eId, len(alignL), partCount, taxCount)
                elif partCount == taxCount:
                    pass

                elif partCount < taxCount:
                    logger.warning("%s partCount %d taxCount %d", eId, partCount, taxCount)
                # -
            # ----
            mU = MarshalUtil()
            ok = mU.doExport(filePath, seqDict, fmt=fmt)
        except Exception as e:
            logger.exception("Failing %r with %s", filePath, str(e))
        return ok
