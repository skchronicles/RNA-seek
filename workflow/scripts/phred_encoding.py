#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import print_function, division
import sys, gzip

# USAGE
# sys.argv[1] = sample_name.R1.fastq.gz

# EXAMPLE
# $ python phred_encoding.py input.R1.fastq.gz

# ABOUT
# Older FastQ files may have quality scores that are encoded with Phred-64.
# This can cause problems with BBmerge, and can cause it to error out. As so,
# this script is used to infer the encoding type to pass that information to
# the qin option of bbmerge:

# @J00170:88:ANYVJBBXX:8:1101:1600:1244 1:N:0:ACTTGA
# GGGAAGTTGAAAGCTTCCAGTGCTCCCTGTCAATTCTAGTCCCTCCAGTCT
# +
# AAAFFJJFJJJJJJFJJJJJJJJJJFJAJJJJJFJJJJJFFJJAJJJJ7JJ <- Determine if Phred-33 encoding or Phred-64

# Current FastQ Format Encodings
# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................   Sanger-33
# ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................   Solexa-64
# ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................   Illumina-64
# .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ.....................   Illumina-64
# LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................   Illumina-33
# NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN...........................................   Nanopore-33
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE......................................   ElemBio-33
# PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP   PacBio-33
# **************************....................................................................   Phred-33
# ........................................................*****************.....................   Phred-64
# !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
# |                         |    |        |              |               |                     |
#  3                        59   64       73             88             104                   126
# 0........................26...31.......40                                                        Sanger-33
#                          -5....0........9.............................40                         Solexa-64
#                                0........9.............................40                         Illumina-64
#                                   3.....9..............................41                        Illumina-64
# 0.2......................26...31........41                                                       Illumina-33
# 0..................20........30........40........50                                              Nanopore-33
# 0..................20........30........40........50...55                                         ElemBio-33
# 0..................20........30........40........50..........................................93  PacBio-33


def usage(message = '', exitcode = 0):
    """Displays help and usage information. If provided invalid usage
    returns non-zero exit-code. Additional message can be displayed with
    the 'message' parameter.
    """
    print('Usage: python {} sampleName.R1.fastq.gz'.format(sys.argv[0]))
    if message:
        print(message)
    sys.exit(exitcode)


def reader(fname):
    """Returns correct file object handler or reader for gzipped
    or non-gzipped FastQ files based on the file extension. Assumes
    gzipped files endwith the '.gz' extension.
    """
    if fname.endswith('.gz'):
        # Opens up file with gzip handler
        return gzip.open
    else:
        # Opens up file normal, uncompressed handler
        return open


def decoded(qscore):
    """Returns Phred ASCII encoding type of FastQ quality scores.
    Older FastQ files may use Phred 64 encoding.
    See the Encoding section here for more information:
    https://en.wikipedia.org/wiki/FASTQ_format
    """
    encoding = ''
    # Unique set of characters across both Phred encoding types
    encodings = { # Pred-33 Encoding characters
                '!': '33', '#': '33', '"': '33', '%': '33', '$': '33', "'": '33',
                '&': '33', ')': '33', '(': '33', '+': '33', '*': '33', '-': '33',
                ',': '33', '/': '33', '.': '33', '1': '33', '0': '33', '3': '33',
                '2': '33', '5': '33', '4': '33', '7': '33', '6': '33', '9': '33',
                '8': '33', ';': '33', ':': '33', '=': '33', '<': '33', '?': '33',
                '>': '33',
                 # Pred-64 Encoding characters, decreased range
                 # due to overlap in ElemBio AVITI Phred+33 scores
                'Y': '64', 'Z': '64', '[': '64', ']': '64', '\\': '64', '_': '64', 
                '^': '64', 'a': '64', '`': '64', 'c': '64', 'b': '64', 'e': '64', 
                'd': '64', 'g': '64', 'f': '64', 'i': '64', 'h': '64'
                }

    for char in qscore:
        try:
            encoding = encodings[char]
            break
        except KeyError:
            pass

    return encoding


if __name__ == '__main__':

    # Check Arguments
    if '-h' in sys.argv or '--help' in sys.argv or '-help' in sys.argv:
        usage(exitcode = 0)
    elif len(sys.argv) != 2:
        usage(message = 'Error: failed to provide all required positional arguments!', exitcode = 1)

    # Get file name
    filename = sys.argv[1]

    # Set handler for gzipped or uncompressed file
    handle = reader(filename)
    # Default encoding if not found
    encoding = '33'

    # Open in 'rt' mode to maintain compatibility across python2 and python3
    # python3 default mode is 'rb' and will return a byte string representation
    with handle(filename, 'rt') as fastq:
        i = 0
        for line in fastq:
            line = line.strip()
            if i%4 == 3: # Quality scores
                encoded = decoded(line)
                if encoded:
                    # Found Phred ASCII encoding type (33 vs. 64)
                    encoding = encoded
                    break # Stop Iteration
            i+=1

    # Print encoding to standard output
    print(encoding)
