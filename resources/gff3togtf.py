#!/usr/bin/env python

"""
Name: gff3togtf.py
Created by: Dr. Tovah Markowitz, NCBR
Date: 11/22/21
Purpose:
    The script converts a NCBI Nucleotide GFF3 into a GTF 
    that will work with RSEM and qualimap. In the gtf, each 
    gene will have at least one transcript and at least one
    exon. Transcripts will be defined by ncRNA or mRNA when 
    available or by CDS otherwise. Also, adds the necessary 
    information for column 9 to match up the multiple rows 
    for each gene and other required functionality.
"""

from __future__ import print_function, division
import re
import argparse


def readGff3(inputName):
    f = open(inputName,'r')
    inputData = f.readlines()
    f.close()
    inputData = [ row.strip().split('\t') for row in inputData ]
    header = [ row[0] for row in inputData if row[0].startswith("#") ]
    inputData = [ row for row in inputData if not row[0].startswith("#") ]
    inputData = [ row for row in inputData if len(row) == 9 ]
    return(header, inputData)


def makeGTF(header, inputData):
    genes = [ row for row in inputData if row[2] == "gene" ]
    RNAs = [ row for row in inputData if (row[2] == "mRNA") | (row[2] == "ncRNA") ]
    CDS = [ row for row in inputData if row[2] == "CDS" ]
    exons = [ row for row in inputData if row[2] == "exon" ]
    outAll = header
    for i in range(len(genes)):
        column9 = genes[i][8].split(';')
        geneID = column9[0].split("-",1)[1]
        geneName = [ pos.split("=")[1] for pos in column9 if pos.startswith("Name=") ][0]
        geneBiotype = [ pos.split("=")[1] for pos in column9 if pos.startswith("gene_biotype") ][0]
        outGene = genes[i][0:8]
        outGene.append('gene_id "' + geneID + '"; gene_name "' + geneName + '"; gene_biotype "' + 
                       geneBiotype + '";')
        outAll.append( "\t".join(outGene) )
        geneRNAs = [ RNA for RNA in RNAs if re.search(geneID + ";", RNA[8]) ]
        if len(geneRNAs) != 0:
            geneTranscript = geneRNAs
        else:
            geneTranscript = [ CDSrow for CDSrow in CDS if re.search(geneID + ";", CDSrow[8]) ]
        geneExons = [ exon for exon in exons if re.search(geneID + ";", exon[8]) ]
        if len(geneTranscript) != 0:
            for row in geneTranscript:
                transcriptRow = row[0:8]
                outTranscript = transcriptRow
                outTranscript[2] = "transcript"
                tmp = row[8].split(";")
                transcriptID = tmp[0].split("-",1)[1]
                if row[2] == "CDS":
                    transcriptName = [ pos.split("=")[1] for pos in tmp if pos.startswith("Name=") ][0]
                else:
                    transcriptName = [ pos.split("=")[1] for pos in tmp if pos.startswith("gene=") ][0]
                outTranscript.append('gene_id "' + geneID + '"; gene_name "' + geneName +
                             '"; gene_biotype "' + geneBiotype + '"; transcript_id "' +
                             transcriptID + '"; transcript_name "' + transcriptName + '";')
                outAll.append( "\t".join(outTranscript) )
        else:
            transcriptRow = genes[i][0:8]
            outTranscript = transcriptRow
            outTranscript[2] = "transcript"
            outTranscript.append('gene_id "' + geneID + '"; gene_name "' + geneName +
                             '"; gene_biotype "' + geneBiotype + '"; transcript_id "' +
                             geneID + '"; transcript_name "' + geneName + '";')
            outAll.append( "\t".join(outTranscript) )
        if len(geneExons) != 0:
            for row in geneExons:
                outExon = row[0:8]
                tmp = row[8].split(";")
                transcriptID2 = [ pos.split("=")[1] for pos in tmp if pos.startswith("Parent=") ][0].split("-",1)[1]
                outExon.append('gene_id "' + geneID + '"; transcript_id "' +
                             transcriptID2 + '"; gene_biotype "' + geneBiotype + '";')
                outAll.append( "\t".join(outExon) )
        else:
            outExon = transcriptRow[0:8]
            outExon[2] = "exon"
            if len(geneTranscript) == 0:
                # Cannot find transcript ID,
                # set to gene id
                transcriptID = geneID
            outExon.append('gene_id "' + geneID + '"; transcript_id "' + transcriptID + '"; gene_biotype "' + geneBiotype + '";')
            outAll.append( "\t".join(outExon) )
    return(outAll)


def writeGTF(GTF, outputName):
    f = open(outputName, 'w')
    f.write( "\n".join(GTF) )
    f.close()


if __name__ == "__main__":

    descriptionText = """
    The script converts a NCBI Nucleotide GFF3 into a GTF
    that will work with RSEM and qualimap. In the gtf, each 
    gene will have at least one transcript and at least one 
    exon. Transcripts will be defined by ncRNA or mRNA when 
    available or by CDS otherwise. Also, adds the necessary 
    information for column 9 to match up the multiple rows 
    for each gene and other required functionality.
    """

    parser = argparse.ArgumentParser(description = descriptionText, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", action="store", required="true", dest="inputName", help="Input GFF3 file name.")
    parser.add_argument("-o", action="store", required="true", dest="outputName", help="Output GTF file name.")
    args = parser.parse_args()

    inputName = args.inputName
    outputName = args.outputName

    (header, inputData) = readGff3(inputName)
    GTF = makeGTF(header, inputData)
    writeGTF(GTF, outputName)

