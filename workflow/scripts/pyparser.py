#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, re
import pandas as pd

# Configuration for defining valid files, cleaning sample names, parse fields, rename fields
# Add new files to parse and define their specifications below
config = {
    ".warning": ["\033[93m", "\033[00m"], ".error": ["\033[91m", "\033[00m"],

    ".rnaseq": {
        ".default": {
            ".output_preference": [
                "Sample", "Encoding", "total_read_pairs", "trimmed_read_pairs",
                "avg_sequence_length", "sequence_length", "gc_content", "percent_duplication",
                "percent_aligned","inner_distance_maxima", "median_insert_size", "mean_insert_size", "mean_mapping_quality",
                "mean_coverage", "avg_aligned_read_length", "pct_mrna_bases", "pct_coding_bases",
                "pct_intronic_bases", "pct_utr_bases", "pct_intergenic_bases", "median_cv_coverage",
                "median_5prime_to_3prime_bias", "median_5prime_bias", "median_3prime_bias",
                "rRNA_percent_aligned", "uni_vec_percent_aligned", "percent_antisense_strand",
                "percent_sense_strand", "median_tin", "flowcell_lanes"
            ]
        }
    },

	"multiqc_cutadapt.txt": {
        "delimeter": "\t",
		"clean_sample_name": ["\.R1$", "\.R2$"],
		"parse_column": ["Sample", "pairs_processed", "r_processed"],
		"rename_field": {
			"pairs_processed": "total_read_pairs",
            "r_processed": "total_read_pairs"
		},
        "typecast": {
            "total_read_pairs": int
        }
	},

	"multiqc_fastqc.txt": {
        "delimeter": "\t",
		"clean_sample_name": ["^QC \\| ", "^rawQC \\| ", "\.trim$", "\.R1$", "\.R2$"],
        "collapse": True,
		"parse_column": ["Sample", "Encoding", "Total Sequences", "Sequence length", "%GC", "avg_sequence_length"],
		"rename_field": {
			"Total Sequences": "trimmed_read_pairs",
			"Sequence length": "sequence_length",
            "%GC": "gc_content",
		},
        "typecast": {
            "trimmed_read_pairs": int,
            "avg_sequence_length": float
        }
	},

	"multiqc_fastq_screen.txt": {
        "delimeter": "\t",
		"clean_sample_name": ["^FQscreen \\| ", "^FQscreen2 \\| ", "_screen$", "\.trim$", "\.R1$", "\.R2$"],
		"parse_column": ["Sample", "Uni_Vec percentage", "rRNA percentage", "Human percentage", "Mouse percentage", "Bacteria percentage", "Fungi percentage", "Virus percentage"],
		"rename_field": {
			"Uni_Vec percentage": "uni_vec_percent_aligned",
			"rRNA percentage": "rRNA_percent_aligned",
            "Human percentage": "human_percent_aligned",
            "Mouse percentage": "mouse_percent_aligned",
            "Bacteria percentage": "bacteria_percent_aligned",
            "Fungi percentage": "fungi_percent_aligned",
            "Virus percentage": "virus_percent_aligned"
		},
        "typecast": {
            "uni_vec_percent_aligned": float,
            "rRNA_percent_aligned": float,
            "human_percent_aligned": float,
            "mouse_percent_aligned": float,
            "bacteria_percent_aligned": float,
            "fungi_percent_aligned": float,
            "virus_percent_aligned": float
        }
	},

	"multiqc_picard_dups.txt": {
        "delimeter": "\t",
		"clean_sample_name": ["\.p2$"],
		"parse_column": ["Sample", "PERCENT_DUPLICATION"],
		"rename_field": {
			"PERCENT_DUPLICATION": "percent_duplication"
		},
        "typecast": {
            "percent_duplication": float
        },
        "scaling_factor": {
            "percent_duplication": 100.0
        },
	},

	"multiqc_picard_RnaSeqMetrics.txt": {
        "delimeter": "\t",
		"clean_sample_name": ["\.p2$"],
		"parse_column": ["Sample", "PCT_CODING_BASES", "PCT_MRNA_BASES", "MEDIAN_CV_COVERAGE", "PCT_INTRONIC_BASES", "MEDIAN_3PRIME_BIAS", "MEDIAN_5PRIME_BIAS", "MEDIAN_5PRIME_TO_3PRIME_BIAS", "PCT_INTERGENIC_BASES", "PCT_UTR_BASES"],
		"rename_field": {
			"PCT_CODING_BASES": "pct_coding_bases",
			"PCT_MRNA_BASES": "pct_mrna_bases",
			"MEDIAN_CV_COVERAGE": "median_cv_coverage",
			"PCT_INTRONIC_BASES": "pct_intronic_bases",
			"MEDIAN_3PRIME_BIAS": "median_3prime_bias",
			"MEDIAN_5PRIME_BIAS": "median_5prime_bias",
			"MEDIAN_5PRIME_TO_3PRIME_BIAS": "median_5prime_to_3prime_bias",
			"PCT_INTERGENIC_BASES": "pct_intergenic_bases",
			"PCT_UTR_BASES": "pct_utr_bases"
		},
        "typecast": {
            "pct_coding_bases": float,
            "pct_mrna_bases": float,
            "median_cv_coverage": float,
            "pct_intronic_bases": float,
            "median_3prime_bias": float,
            "median_5prime_bias": float,
            "median_5prime_to_3prime_bias": float,
            "pct_intergenic_bases": float,
            "pct_utr_bases": float
        }
	},

	"multiqc_rseqc_infer_experiment.txt": {
        "delimeter": "\t",
		"clean_sample_name": ["^RSeQC \\| ", "\.strand\.info$","\.info\.strand$", "^output\.", "\.p2$"],
		"parse_column": ["Sample", "pe_sense", "se_sense", "pe_antisense", "se_antisense"],
		"rename_field": {
			"pe_sense": "percent_sense_strand",
            "se_sense": "percent_sense_strand",
			"pe_antisense": "percent_antisense_strand",
            "se_antisense": "percent_antisense_strand",
		},
        "typecast": {
            "percent_sense_strand": float,
            "percent_antisense_strand": float
        },
        "scaling_factor": {
            "percent_sense_strand": 100.0,
            "percent_antisense_strand": 100.0
        },
	},

    "rseqc_inner_distances.txt": {
        "delimeter": "\t",
        "clean_sample_name": ["\.inner_distance_freq\.txt$"],
        "parse_column": ["Sample", "Inner_Dist_Maxima"],
        "rename_field": {
            "Inner_Dist_Maxima": "inner_distance_maxima"
        },
        "typecast": {
            "inner_distance_maxima": float
        },
    },

	"rseqc_median_tin.txt": {
        "delimeter": "\t",
		"clean_sample_name": ["\.star_rg_added\.sorted\.dmark\.bam$"],
		"parse_column": ["Sample", "median_tin"],
        "typecast": {
            "median_tin": float
        }
	},

	"sample_group.txt": {
        "delimeter": "\t",
		"clean_sample_name": [],
		"parse_column": ["Sample", "TissueType"],
	},

	"fastq_flowcell_lanes.txt": {
        "delimeter": "\t",
		"clean_sample_name": [],
		"parse_column": ["Sample", "flowcell_lanes"],
	},

	"multiqc_star.txt": {
        "delimeter": "\t",
		"clean_sample_name": ["\.p2$"],
		"parse_column": ["Sample", "uniquely_mapped_percent", "avg_input_read_length"],
		"rename_field": {
			"uniquely_mapped_percent": "percent_aligned",
			"avg_input_read_length": "avg_aligned_read_length"
		},
        "typecast": {
            "percent_aligned": float,
            "avg_aligned_read_length": int
        }
	},

	"multiqc_qualimap_bamqc_genome_results.txt": {
        "delimeter": "\t",
		"clean_sample_name": ["\.p2$"],
		"parse_column": ["Sample", "mean_insert_size", "median_insert_size", "mean_mapping_quality", "mean_coverage"],
		"rename_field": {},
        "typecast": {
            "mean_insert_size": float,
            "median_insert_size": float,
            "mean_mapping_quality": float,
            "mean_coverage": float
        }
	}
}


def help():
        return """
pyparser.py - a config based file parser.

USAGE:
    python pyparser.py <file_1> <file_2> <file_3> <file_N-1> <output_dir> [-h]

    Positional Arguments:
        [1...N-1]       Type [File]: An output file from MultiQC, or a list of
                                   output files generated from MultiQC. Each provided
                                   file is parsed, and information is aggregated
                                   across all samples into a single tab-seperated
                                   ouput file: multiqc_matrix.tsv

            Currently Supported MultiQC Files:
            multiqc_cutadapt.txt, multiqc_star.txt, multiqc_picard_dups.txt,
            multiqc_fastq_screen.txt, multiqc_picard_RnaSeqMetrics.txt,
            multiqc_fastqc.txt, multiqc_rseqc_infer_experiment.txt,
            multiqc_qualimap_bamqc_genome_results.txt

        [N]             Type [PATH]: An existing PATH on the filesystem to write
                                  the output file (i.e. multiqc_matrix.tsv).

    Optional Arguments:
        [-h, --help]  Displays usage and help information for the script.

    Example:
        # Creates QC table: multiqc_matrix.tsv in users working directory
        $ python pyparser.py multiqc_cutadapt.txt multiqc_fastqc.txt multiqc_fastq_screen.txt $PWD
        # Supports globbing
        $ python pyparser.py /path/to/MultiQC/ouput/folder/*.txt $PWD

    Requirements:
        multiqc == 1.9
        python >= 2.7
          + pandas
"""


def args(argslist):
    """Parses command-line args from "sys.argv". Returns a list of filenames to parse."""
    # Input list of filenames to parse
    *files, odir = argslist[1:]

    # Check for optional args
    if '-h' in files or '--help' in files:
        print(help())
        sys.exit(0)
    # Check to see if user provided input files to parse
    elif not files:
        print("\n{}Error: Failed to provide input files to parse!{}".format(*config['.error']), file=sys.stderr)
        print(help())
        sys.exit(1)

    return files, odir


def isvalid(file):
    """Checks if a file is a recognized/supported file to parse.
    Comparse file's name against list of supported files in config specification"""

    # Supported files are keys in config
    supported = config.keys()

    # Remove absolute or relateive PATH
    if os.path.basename(file) not in supported:
        cstart, cend = config['.warning']
        print("{}Warning:{} {} is a not supported file to parse... Skipping over file!".format(cstart, cend, file))
        return False

    return True


def exists(file):
    """Checks to see if file exists or is accessible.
    Avoids problem with inconsistencies across python2.7 and >= python3.4 and
    works in both major versions of python"""

    try:
        fh = open(file)
        fh.close()
    # File cannot be opened for reading (may not exist) or permissions problem
    except IOError:
        cstart, cend = config['.warning']
        print("{}Warning:{} Cannot open {}... File may not exist... Skipping over file!".format(cstart, cend, file))
        return False

    return True


def column_indexes(line, filename, verbose=True):
    """Parses header of file to find fields of interest defined in config[filename][parse_column]
    Returns a list of integers corresponding to an index of a column to parse."""

    indices = []
    found = []
    header = line
    # Remove file's PATH before cross-referencing config
    filename = os.path.basename(filename)
    fields2parse = config[filename]["parse_column"] # Attributes or columns of interest

    # Get index of column to parse
    for i in range(0, len(header), 1):
        if header[i] in fields2parse:
            indices.append(i)
            found.append(header[i])

    if verbose:
        # Warning that an expected field could not be found
        fields_not_found = set(fields2parse) - set(found)
        for field in fields_not_found:
            cstart, cend = config['.warning']
            print("{}Warning:{} Cannot find expected field '{}' in {}... skipping over parsing that field!".format(cstart, cend, field, filename))

    return indices


def clean(linelist, sample_name_index, filename):
    """Cleans sample name from suffixes defined in config[filename]['clean_sample_name'] and
    renames fields defined in config[filename]['rename_field']. Returns a list of cleaned fields."""

    samplename = linelist[sample_name_index]

    # Remove file's PATH before cross-referencing config
    filename = os.path.basename(filename)
    for suffix in config[filename]['clean_sample_name']:
            regex = '{}'.format(suffix)
            samplename = re.sub(regex, '', samplename)

    # Update linelist with new sample name
    linelist[sample_name_index] = samplename

    return linelist


def rename(header, filename):
    """Renames fields defined in config[filename]['rename_field']. Returns a list of re-named columns."""
    renamed = []
    filename = os.path.basename(filename)
    for field in header:
        try:
            newname = config[filename]['rename_field'][field]
            renamed.append(newname)
        # Field is not in config, keep old name
        except KeyError:
            renamed.append(field)
    return renamed


def cast_typed(value, column, filename, decimals=3):
    """Cast types data in a row/column according to specification in config[filename]["typecast"][column_name].
    Converts string to either an integer or float (rounded to three decimal places).
    """
    filename = os.path.basename(filename)
    try:
        # Python witch-craft, functions are first-class objects and can be used accordingly
        # Storing function object into caster variable for typecasting as int() or float()
        caster = config[filename]["typecast"][column]
        value = caster(value)        # typecast to spec defined in config
        if type(value) is float:
            value = round(value, decimals)
    except ValueError:
        # Must convert to float before converting to integer
        # cannot pass a string representation of a float into int()
        if value: # case for when row/column is empty string
            value = caster(float(value))
    except KeyError:
        # No type is defined in config, pass
        pass
    return value


def scaled(value, column, filename):
    """Scales data in a row/column according to specification in config[filename]["scaling_factor"][column_name].
    User-defined optional multiplier
    """
    filename = os.path.basename(filename)
    try:
        # Get the scaling factor
        scaling_unit = config[filename]["scaling_factor"][column]  # KeyError if DNE
        value = value * scaling_unit  # TypeError if string
        value = round(value, 3)
    except TypeError:
        # Did not typecast value using the config
        # Remove warning by typecasting value as float or int
        cstart, cend = config['.warning']
        print("{}Warning:{} Attribute {} in {} is NOT defined in config... defaulting to float".format(cstart, cend, column, filename))
        if value: # case for when row/column is empty string
            value = float(value) * scaling_unit
            value = round(value, 3)
    except KeyError:
        pass
    return value


def populate_table(parsed_header, parsed_line, file, data_dict):
    """Appends parsed sample metadata to a nested dictionary where
    dictionary['Sample_Name']['QC_Attribute'] = QC_Metadata.
    Returns an updated dictionary containing new information for N-th line."""

    sample_index = parsed_header.index('Sample')
    sample_name = parsed_line[sample_index]

    # Add sample name to dictionary, if does not exist [key1]
    if parsed_line[sample_index] not in data_dict:
        data_dict[sample_name] = {}

    for i in range(0, len(parsed_line), 1):
        # Skip over sample name (already first key)
        if parsed_line[i]: # check if empty string
            metadata = cast_typed(parsed_line[i], parsed_header[i], file)
            metadata = scaled(metadata, parsed_header[i], file)
            data_dict[sample_name][parsed_header[i]] = metadata

    return data_dict


def parsed(file, delimeter='\t'):
    """Parses columns of file according to specification in config[filename]['parse_column'].
    Column names are renamed according to specification in config[filename]['rename_field'].
    Sample names are cleaned to removed any prefixes or suffixes specified in config[filename]['clean_sample_name'].
    Yields a tuple consisting of the parsed header and N-th parsed line of the file.
    """

    #print('\nBeginning to parse {}'.format(file))
    with open(file, 'r') as fh:
        # Parse header
        header = next(fh).strip().split(delimeter) # Get file header
        indexes = column_indexes(header, file)     # Indexes of columns to parse
        header = [header[i] for i in indexes]  # Parse each column of interest
        header = rename(header, file)   # Rename columns

        # Parse QC metadata from file
        sample_index = header.index('Sample')
        for line in fh:
            #linelist = line.strip().split(delimiter)
            linelist = line.rstrip('\n').split(delimeter)
            parsed_line = [linelist[i] for i in indexes]
            parsed_line = clean(parsed_line, sample_index, file) # remove extensions from sample name

            yield header, parsed_line

def main():

    # Minor Todo(s):
    #       1. Get rid of pandas dependency (add transpose function and loop through dict to print table)
    #       2. Add more advanced argument parsing, make path to config an arg

    # Check for usage and optional arguements, get list of files to parse and output directory
    ifiles, outdir = args(sys.argv)

    # Check if files are supported, see config specification, and if file is readable
    ifiles = [file for file in ifiles if isvalid(file) and exists(file)]

    # Parse each file and add to the QC metadata dicitionary
    QC = {}
    for file in ifiles:
        for header, line in parsed(file):
            QC = populate_table(header, line, file, QC)

    df = pd.DataFrame(QC).transpose()

    # Get default output peference
    try:
        output_preference = config['.rnaseq']['.default']['.output_preference']
        df = df.reindex(columns = output_preference)
    except KeyError:
        # Output peference is not defined in config
        pass

    # Write to file
    df.to_csv(os.path.join(outdir, 'multiqc_matrix.tsv'), index = False, sep='\t')


if __name__ == '__main__':

    main()
