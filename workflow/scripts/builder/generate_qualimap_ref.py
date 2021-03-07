import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import Bio.SeqUtils as SeqUtils
import HTSeq
import numpy as np


def idsContainGiven(givenId, transcriptIds):
    for tId in transcriptIds:
        if givenId.find(tId) != -1:
            return True

    return False

if __name__ == "__main__":
    
    descriptionText = "The script extracts features from a GTF file and a FASTA file into Qualimap annotation format. Note: exons have to be sorted according to exon number! This important for correct reverse transcribed cDNA sequences extraction."

    parser = argparse.ArgumentParser(description = descriptionText,formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-g", action="store", required="true", dest="gtfFile",
        help="Input file with list of genes in GTF format")
    parser.add_argument("-f", action="store", required="true", dest="fastaFile",
        help="Input genome sequence. ")
    parser.add_argument("-o", action="store", dest="outFile", default="annotations.txt",
        help="Output file. Default is annotations.txt")

    parser.add_argument("--filter", action="store", dest="filterStr", default="",
                        help="Comma-separted list of entries to filter from GTF file \
                        based on given attribute id")

    parser.add_argument("--ignore-strange-chrom", action="store_true", default=False,
        dest="ignoreStrangeChromosomes", help="All chromosomes except numbered and X,Y,MT are ignored ")
    
    args = parser.parse_args()
    
    print args
    
    gtfFileName = args.gtfFile
    fastaFileName = args.fastaFile
    outFileName = args.outFile
    attr_id = "gene_id"

    # parse GTF file

    gtf_file = HTSeq.GFF_Reader( gtfFileName )

    features = {}

    filtered_transcripts = args.filterStr.split(",")
    if filtered_transcripts:
        print "Filtering for: ", filtered_transcripts

    for feature in gtf_file:
        if feature.type == 'exon':
            geneName = feature.attr[ attr_id ]
            #print transcriptName
            if geneName in features:
                features[geneName].append(feature)
            else:
                features[geneName] = [feature]

    # load & save sequences 

    seqData = SeqIO.to_dict(SeqIO.parse(fastaFileName, "fasta"))
    
    outFile = open(outFileName, "w")

    header = "\"%s\"\t\"%s\"\t\"%s\"\n" % ("biotypes","length","gc")
    outFile.write(header)


    for geneId in features:
       
        exons = features[geneId]

        print "Processing %s" % geneId

        if len(exons) == 0:
            continue
        try:
        	biotype = exons[0].attr["gene_type"]
	except KeyError:
		biotype = exons[0].attr["gene_biotype"]
        length = 0
        transcripts = {}

        for exon in exons:
            transcriptId = exon.attr["transcript_id"]
            
            tSeq = transcripts.get(transcriptId, Seq(""))

            iv = exon.iv
            seqName = iv.chrom
            if seqName in seqData:
                #print "Exon (%s,%d,%d) " % (iv.chrom,iv.start,iv.end)
                buf = seqData[ iv.chrom ].seq[ iv.start  : iv.end ]
                if iv.strand == '-':
                    buf = buf.reverse_complement()
                tSeq += buf
            else:
                print "Can not locate sequence  %s in %s, skipping region..." % (seqName, fastaFileName)
            transcripts[transcriptId] = tSeq

        gc_array = []
        lengths = []

        for tSeq in transcripts.values():
            lengths.append( len(tSeq) )
            gc_array.append ( SeqUtils.GC(tSeq) )

        #gene_length = len(seq_rec) / len(transcripts)
        #gc = SeqUtils.GC(seq_rec)
        
        #print gc_array, lengths
        gene_length = np.mean(lengths)
        gene_gc = np.mean(gc_array)
        
        #print len(seq_rec), len(transcripts),gc

        line = "\"%s\"\t\"%s\"\t%d\t%.2f\n" % (geneId, biotype, gene_length, gene_gc)

        outFile.write ( line )
                   
        #outFile.flush()
        #sys.stdin.readline()

    outFile.close()

