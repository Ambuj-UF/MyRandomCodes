
from Bio import SeqIO
import os
from Bio.Blast.Applications import NcbitblastnCommandline
import re

from Bio.Seq import Seq

from Bio.Alphabet import generic_dna


def merge_dict(d1, d2):
    if d1 == {}:
        return d2
    else:
        newDict = dict()
        for key, val in d1.items():
            newDict[key] = (val)
            for inkey, inval in d2.items():
                if key == inkey:
                    for childKey, childVal in inval.items():
                        newDict[key][childKey] = (childVal)
            
                if inkey not in d1.keys():
                    newDict[inkey] = (inval)

    return newDict


def tblastnWrapper(rec, seqObj):
    retData = dict()
    store_blasted = dict()
    

    with open("infile.fas", "w") as fp:
        SeqIO.write([rec], fp, "fasta")

    os.system("makeblastdb -in infile.fas -out blast/HIV1_O_blast -dbtype nucl -hash_index")
    tblastn_cline = NcbitblastnCommandline(cmd='tblastn',
                                                   db="blast/HIV1_O_blast",
                                                   num_threads=6,
                                                   outfmt=5,
                                                   seg="no",
                                                   out=("Result.xml")
                                                   )
                

    stdout, stderr = tblastn_cline(stdin = seqObj)
                                                                   
    """Parse Blast output"""
    with open("Result.xml",'r') as xml:
        cFlag = False; eFlag = False; e_val = list(); negFlag = False; fcount = 0
        hitFlag = False
        for line in xml:
            if re.search('No hits found', line) == None:
                """Check if the sequence belong to the group"""
                cFlag = True
                                                                                                   
            if re.search('<Hit_id>', line) != None:
                hitFlag = True
                                                                                                           
            if hitFlag == True:
                if re.search('<Hit_def>', line) != None:
                    description = line.strip().strip("\n")[9:-10]
                    retData["description"] = description
                if re.search('<Hsp_query-from>', line) != None:
                    retData["qstart"] = int(line.strip().rstrip().strip("<Hsp_query-from>").strip('</'))
                if re.search('<Hsp_query-to>', line) != None:
                    retData["qstop"] = int(line.strip().rstrip().strip("<Hsp_query-to>").strip('</'))
                if re.search('<Hsp_hit-from>', line) != None:
                    retData["start"] = int(line.strip().rstrip().strip("<Hsp_hit-from>").strip('</'))
                if re.search('<Hsp_hit-to>', line) != None:
                    retData["stop"] = int(line.strip().rstrip().strip("<Hsp_hit-to>").strip('</'))
                if re.search('<Hsp_evalue>', line) != None:
                    retData["eval"] = float(line.strip().rstrip().strip("<Hsp_evalue>").strip('</'))
                if re.search('<Hsp_query-frame>', line) != None:
                    retData["frame"] = int(line.strip().rstrip().strip("<Hsp_query-frame>").strip('</'))
                    if retData["frame"] < 0:
                        retData["start"], retData[rec.id][id]["stop"] = reverser(retData[rec.id][id]["start"], retData[rec.id][id]["stop"])
                if re.search('<Hsp_hseq>', line) != None:
                    retData["seq"] = line.strip().rstrip().strip("<Hsp_hseq>").strip('</')
                    hitFlag = False
                            
        store_blasted = merge_dict(store_blasted, retData)
    
    retData = store_blasted
    
    return retData


def run_blast(filename, seqObj, idObj):
    print seqObj
    record = list(SeqIO.parse(open(filename, "rU"), "fasta"))
    newrec = list()
    store = list()
    for i, rec in enumerate(record):
        seqObjrec = "".join([x for x in str(rec.seq) if x != "-"])
        rec.seq = Seq(seqObjrec, generic_dna)
        output = tblastnWrapper(rec, seqObj)
        if output.keys() == []:
            continue

        record[i].seq = rec.seq[int(output["start"])-1: int(output["stop"])]
        aminoSeq = record[i].seq.translate()
        print seqObj
        print aminoSeq
        print len(seqObj), len(aminoSeq)
        if len(seqObj) == len(aminoSeq) and "*" not in aminoSeq:
            newrec.append(record[i])

    try:
        os.mkdir("Sequences/"+idObj)
    except:
        pass

    with open("Sequences/"+ idObj + "/seq.phy", "w") as fp:
        SeqIO.write(newrec, fp, "phylip-relaxed")
    
    with open("Sequences/"+ idObj + "/seq.fas", "w") as fp:
        SeqIO.write(newrec, fp, "fasta")

    with open("Sequences/"+ idObj + "/seq.nex", "w") as fp:
        SeqIO.write(newrec, fp, "nexus")

    print idObj, "Length = ", len(newrec)



allRecord = list(SeqIO.parse(open("zika_sequences.txt", "rU"), "fasta"))
for allRec in allRecord:
    if allRec.id == "AMM39804.1_All":
        continue
    run_blast("genome_1.fas", str(allRec.seq), str(allRec.id))

