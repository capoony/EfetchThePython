from optparse import OptionParser, OptionGroup
from Bio import Entrez
import sys
import json
from collections import defaultdict as d

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, 'Requires BioPython: "pip install biopython"')

#########################################################   CODE   #########################################################################

parser.add_option("--Term", dest="Term",
                  help="Search Term in Quotes, e.g. \"Torpediniformes[Organism], COI\"")
parser.add_option("--Output", dest="OUT", help="Prefix of output files")
parser.add_option("--email", dest="email",
                  help="NCBI email ID; see here https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/")
parser.add_option("--api_key", dest="api",
                  help="api key for NCBI efetch, see here: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/")
parser.add_option("--FASTA", dest="FA",
                  help="optionally produce FASTA file", action="store_true")


(options, args) = parser.parse_args()
parser.add_option_group(group)

Entrez.email = options.email
Entrez.api_key = options.api

handle = Entrez.esearch(db="nuccore",
                        term=options.Term,
                        idtype="acc")
record = Entrez.read(handle)
handle.close()
count = record['Count']
# print(count)
handle = Entrez.esearch(db="nuccore",
                        term=options.Term,
                        idtype="acc",
                        retmax=count)
records = Entrez.read(handle)
handle.close()
# print(len(records["IdList"]))

Out1 = open(options.OUT+".tsv", "wt")
if options.FA:
    Out2 = open(options.OUT+".fasta", "wt")
Out1.write("GenBank_ID\tAccession\tTaxName\tGenBankTitle\tJournal\tSubmissionTitle\tConsortium\tLat_Lon\tCollCountry\tTissue\tIdentifier\tCollector\tSeq\n")
C = 1
for id in records["IdList"]:

    handle2a = Entrez.efetch(db="nuccore",
                             id=id,
                             rettype="gb",
                             retmode="xml")

    TMPa = Entrez.read(handle2a)[0]
    for x, y in TMPa.items():
        QualDict = d(list)
        if x == "GBSeq_feature-table":
            for j in y:
                for x1, y1 in j.items():
                    if x1 == "GBFeature_quals":
                        for k in y1:
                            for x2, y2 in k.items():
                                if x2 == "GBQualifier_name":
                                    ID = y2
                                elif x2 == "GBQualifier_value":
                                    if ID == "PCR_primers":
                                        QualDict[ID].append(y2)
                                    else:
                                        QualDict[ID] = y2

        elif x == "GBSeq_references":
            RefDict = d(lambda: d(str))
            for j in y:
                for x1, y1 in j.items():
                    if x1 == "GBReference_reference":
                        ID = y1
                        continue
                    RefDict[ID][x1.split("_")[1]] = y1

    handle2a.close()

    handle2b = Entrez.esummary(db="nuccore",
                               id=id,
                               retmode="json",
                               version="2.0")
    TMPb = json.load(handle2b)
    handle2b.close()
    for UID in TMPb["result"]["uids"]:
        TMP2 = TMPb["result"][UID]
        SUB = dict(zip(TMP2["subtype"].split("|"),
                       TMP2["subname"].split("|")))

        SUB = d(str, SUB)
        for k, v in RefDict.items():
            Out1.write("\t".join([str(x) for x in [TMP2["gi"],
                       TMPa["GBSeq_accession-version"],
                       TMPa["GBSeq_organism"],
                       TMP2["title"].replace(" ", "_"),
                       RefDict[k]["journal"],
                       RefDict[k]["title"],
                       RefDict[k]["consortium"],
                       SUB["lat_lon"],
                       SUB["country"],
                       SUB["tissue_type"],
                       SUB["identified_by"],
                       SUB["collected_by"],
                       TMPa["GBSeq_sequence"]]])+"\n")
        if options.FA:
            Out2.write(">"+TMPa["GBSeq_organism"].replace(" ", "_")+" " +
                       TMPa["GBSeq_accession-version"]+"\n")
            Out2.write(TMPa["GBSeq_sequence"]+"\n")
    C += 1
