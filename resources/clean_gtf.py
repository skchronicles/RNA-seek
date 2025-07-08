#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Author: Skyler Kuhn

# Standard Library
from __future__ import print_function
from datetime import datetime
from textwrap import dedent
import argparse, gzip, os, re, sys


# Constants
_VERSION = '1.0.0'
_NAME    = os.path.basename(sys.argv[0])
_HELP    = dedent("""
@Usage:
    $ ./{0} [-h] [--version] \\
            [--gene-name-attribute GENE_NAME] \\
            --input-gtf INPUT_GTF_FILE \\
            --output-gtf OUTPUT_GTF_FILE
@About:
    This script cleans a GTF file by adding
    extra required fields in the 9th column
    of a GTF file. The 9th column contains
    additional metadata in the form of key,
    value pairs. This script will also remove
    any nested quotes in the 9th column that
    may cause issues with other programs or
    tools that parse GTF files (e.g. STAR,
    RSEM, and other custom scripts to build
    reference files).

    It is recommended to run this program with
    a version of Python greater than '3.8'.
    This ensures that the order of the key,
    value pairs in the 9th column is retained.

@Required:
    -i, --input-gtf INPUT_GTF_FILE
        Input GTF file to be cleaned. It is
        recommended to use a GTF file that
        has been pre-processed by AGAT. AGAT
        can be used to convert a GFF file to
        GTF file format, and it can also fix
        most common issues with GTF files. If
        you have a GTF file from a provider
        outside of GENCODE, we highly recommend
        running it through AGAT first. Please
        see the example below for how you can
        pre-process a GFF/GTF with AGAT before
        running this script. If you downloaded
        the GTF file from GENCODE, you do not
        need to pre-process it with AGAT or
        this script. It is already in the
        correct format.
    -o, --output-gtf OUTPUT_GTF_FILE
        Output GTF file with cleaned metadata.
        This file will contain any required
        key, value pairs in the 9th column of
        the GTF file and it will also have any
        nested quotes removed from the 9th
        column. This output file can be used
        as an input file to the build sub-
        command of the RNA-seek pipeline.
@Options:
    -g, --gene-name-attribute GENE_NAME
        This option allows you to specify
        the attribute to use for the gene
        name in the GTF file. By default,
        this script will use the 'gene_name'
        attribute in the 9th column of the
        GTF file. If the 'gene_name' attribute
        is not present in the GTF file, it
        will use the 'gene_id' attribute as
        the gene name. If you want to use a
        different attribute for the gene name,
        you can specify it with this option.
        This maybe needed for annotations
        coming from NCBI or GenBank. Some
        annotations from NCBI have 'gene'
        as the attribute for the gene name.
        In that case, you can use set this
        option to 'gene'.
          • Default: 'gene_name'
    -h, --help
        Shows help message and exits.
    -v, --version
        Prints the version and exits.

@Example:
    # Steps for converting messy a GFF/GTF file
    # into a properly formatted GTF file
    # 1. Pull the docker image from their registry
    #    and create SIF
    module load singularity
    SINGULARITY_CACHEDIR=${{PWD}}/.${{USER}} \\
    singularity pull -F \\
    docker://quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0

    # 2. Run AGAT todo the heavy lifting and correct
    #    any major issues. AGAT can take a GTF or GFF
    #    file as input, but we recommend to providing
    #    a GTF file to AGAT if you have access to it.
    singularity exec -B $PWD \\
    agat_1.4.2--pl5321hdfd78af_0.sif agat_convert_sp_gff2gtf.pl \\
        --gtf_version 3 \\
        --gff input.gtf \\
        -o converted.gtf

    # 3. Finally run this script to add any required
    #    metadata to the 9th column of the GTF file
    #    and to remove any nested quotes.
    module load python/3.9
    ./{0} \\
        --gene-name-attribute gene_name \\
        --input-gtf converted.gtf \\
        --output-gtf clean.gtf
""".format(_NAME))


# Helper functions
def err(*message, **kwargs):
    """Prints any provided args to standard error.
    kwargs can be provided to modify print functions
    behavior.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    print(*message, file=sys.stderr, **kwargs)


def fatal(*message, **kwargs):
    """Prints any provided args to standard error
    and exits with an exit code of 1.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    err(*message, **kwargs)
    sys.exit(1)


def timestamp(format="%Y-%m-%d %H:%M:%S"):
    """Returns a formatted timestamp string
    for the current time.
    @param format <str>:
        Format string for the timestamp, default:
        "%Y-%m-%d %H:%M:%S" which is equivalent to
        "2023-10-01 12:00:00" for example.
    @return <str>:
        Formatted timestamp string, i.e. "2023-10-01 12:00:00"
    """
    return datetime.now().strftime(format)


def log(*message):
    """Logs a message to standard output with a timestamp.
    @param message <any>:
        Values printed to log
    """
    print("[{0}] {1}".format(
        timestamp(),
        " ".join([str(m) for m in message]))
    )


def check_permissions(parser, path, *args, **kwargs):
    """Checks permissions using os.access() to see the
    user is authorized to access a file/directory. Checks
    for existence, read, write and execute via args:
        • os.F_OK, tests existence
        • os.R_OK, tests read
        • os.W_OK, tests write
        • os.X_OK, tests exec
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param path <str>:
        Name of path to check
    @param args <any>:
        Positional args to pass to os.access()
    @param kwargs <any>:
        Named kwargs to pass to os.access()
    @return path <str>:
        Returns absolute path if it exists and the
        checked permssions are setup are correct.
    """
    if not os.path.exists(path):
        parser.error(
            "Path '{}' does not exists! Failed to provide vaild input.".format(path)
        )
    if not os.access(path, *args, **kwargs):
        parser.error(
            "Path '{}' exists, but cannot read path due to permissions!".format(path)
        )
    return os.path.abspath(path)


def parse_cli_arguments():
    """Parses command line arguments and returns
    an argparse.parse_args object.
    @return <argparse.parse_args()>:
        Parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        add_help=False,
        description=_HELP,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage = argparse.SUPPRESS,
    )
    # Input GTF file (from AGAT)
    parser.add_argument(
        '-i', '--input-gtf',
        type = lambda file: \
            check_permissions(parser, file, os.R_OK),
        required=True,
        help=argparse.SUPPRESS
    )
    # Cleaned output GTF file
    parser.add_argument(
        '-o', '--output-gtf',
        type=str,
        required=True,
        help=argparse.SUPPRESS
    )
    # Gene name attribute,
    # default: gene_name
    parser.add_argument(
        '-g', '--gene-name-attribute',
        type=str,
        required=False,
        default='gene_name',
        help=argparse.SUPPRESS
    )
    # Get version information
    parser.add_argument(
        '-v', '--version',
        action='version',
        help = argparse.SUPPRESS,
        version='%(prog)s {0}'.format(_VERSION)
    )
    # Add custom help message
    parser.add_argument(
        '-h', '--help',
        action='help',
        help=argparse.SUPPRESS
    )
    return parser.parse_args()


def replace_nested_quotes(line, find_char = '"', replace_char = ''):
    """
    Replaces nested quotes within a string. This function assumes
    the quote character in the 9th column is a double quote or
    <"> character. This is the correct character to use based on
    the GTF specification. This function is used to fix issues
    with nested quotes in the 9th column of a GTF file. This will
    cause problems with other programs or tools that parse
    GTF files downstream of this script.
    @param line <str>:
        The line from the GTF file to process. This should be the
        value in the key, value pairs in the 9th column.
    @param find_char <str>:
        The character to find and replace in the string. Default is
        double quote <"> character. This is the correct character to
        use based on the GTF specification.
    @param replace_char <str>:
        The character to replace the found character with. Default is
        an empty string <''>. This will remove the nested quotes.
    @return fixed <str>:
        Returns the fixed string with nested quotes removed.
    """
    # Normal:
    #   protein_id "XP_040355194.1";
    # Bad:
    #   transl_except "(pos:956284..956286)" "(pos:956290..956292)";
    # Fixed:
    #   transl_except "(pos:956284..956286) (pos:956290..956292)";
    quote_count = 0
    inside_quotes = False
    fixed = ''

    for i in range(len(line)):
        curr_char = line[i]
        # Scan for next character to determine
        # if reached end of quotation.
        try: next_char = line[i+1]
        except IndexError: next_char = ''

        if curr_char  == '"':
            # Entered the border or ending of
            # a quote, increase the counter and
            # check where we are in the string
            quote_count += 1
            if quote_count == 1:
                inside_quotes = True

        if next_char  == ';':
            # Reached end border of quote,
            # reset boolean flag and counters
            inside_quotes = False
            quote_count = 0

        if inside_quotes:
            # Fix evil mistakes of the past,
            # replace reserved delimeter with
            # another character, let's use a
            # url encoding of the character
            if curr_char  == find_char and quote_count > 1:
                curr_char  = replace_char

        # Add the existing/converted character
        fixed += curr_char

    return fixed



def url_escape_inside_quotes(line, delimiter=';', url_encoding = '%3B'):
    """Replaces any reserved delimiter in the 9th column with a URL
    encoding of the character. This is needed to fix issues with any
    reserved character in the 9th column of a GTF file (specifically
    the value in a key, value pair). If this issue is not fixed, it
    will cause problems with other programs or tools that parse GTF
    files downstream of this script. For more information, please see
    the following issue for description and context:
    https://github.com/NBISweden/AGAT/issues/250
    NOTE: This function assumes the quote character in the 9th column
    is a double quote or <"> character. This is the correct character
    to use based on the specification for GTF files.
    @param line <str>:
        The line from the GTF file to process. This should be the
        value in the key, value pairs in the 9th column.
    @param delimiter <str>:
        The character to replace with a URL encoding. Default is
        semicolon <;> character.
    @param url_encoding <str>:
        The URL encoding of the character to replace. Default is
        '%3B' which is the URL encoding of the semicolon <;> character.
    @return fixed <str>:
        Returns the fixed string with reserved delimiters replaced
        with the URL encoding of the character.
    """
    quote_count = 0
    inside_quotes = False
    fixed = ''
    for c in line:
        if c == '"':
            # Entered the border or ending of
            # a quote, increase the counter and
            # check where we are in the string
            quote_count += 1
            inside_quotes = True

            if quote_count > 1:
                # Reached end border of quote,
                # reset boolean flag and counters
                inside_quotes = False
                quote_count = 0

        if inside_quotes:
            # Fix evil mistakes of the past,
            # replace reserved delimeter with
            # another character, let's use a
            # url encoding of the character
            if c == delimiter:
                c = url_encoding

        # Add the existing/converted character
        fixed += c

    return fixed


def stripped(v):
    """Cleans string to remove quotes. This is used to strip a string
    of any quotes that may be present in the value of a key, value pair
    in the 9th column of a GTF file.
    @param v <str>:
        The value to strip of quotes. This should be the value in a
        key, value pair in the 9th column.
    @return v <str>:
        Returns the stripped string with any quotes removed. This is
        used to ensure that the value is clean and does not contain any
        quotes that may cause issues with other programs or tools that
        parse GTF files downstream of this script.
    """
    return v.strip('"').strip("'")


def lookup(mykey, dictionary):
    """ Tries to lookup value in dictionary using an exact match if the
    key. Returns empty string if not found.
    @param mykey <str>:
        The key to lookup in the dictionary. This should be the key in a
        key, value pair in the 9th column of a GTF file.
    @param dictionary <dict>:
        The dictionary to lookup the key in. This should be the parsed
        metadata from the 9th column of a GTF file.
    @return v <str>:
        Returns the value associated with the key in the dictionary. If
        the key is not found, it returns an empty string. This is used
        to ensure that the value is clean and does not contain any quotes
        that may cause issues with other programs or tools that parse GTF
        files downstream of this script.
    """
    v = ''
    if mykey in dictionary:
        v = dictionary[mykey]
        v = stripped(v)
    return v


def contains(pattern, dictionary):
    """Flexible lookup that searches for a pattern instead of a key.
    Returns empty string if pattern is not found in dictionary.
    @param pattern <str>:
        The pattern to search for in the dictionary keys.
        This should be a substring of the key you are looking
        for in the 9th column of a GTF file.
    @param dictionary <dict>:
        The dictionary to search for the pattern in. This should
        be the parsed metadata from the 9th column of a GTF
        file.
    @return v <str>:
        Returns the value associated with the first key that
        contains the pattern in the dictionary. If no key
        contains the pattern, it returns an empty string.
    """
    v = ''
    kys = dictionary.keys()
    for k in kys:
        if pattern in k:
            v = dictionary[k]
            break
    return v


def parse(linelist):
    """Parses key, value pairs in 9th column and returns an index, i.e.
    dictionary, of all the fields.
    @param linelist <list>:
        List of strings representing a line from a GTF file. This is a
        string that has been split on the tab char. The 9th column should
        contain key, value pairs in the form of 'key "value";'.
    @return tags <dict>:
        Returns a dictionary with keys and values parsed from the 9th column
        of the GTF file. This is used to ensure that all required metadata is
        present in the 9th column of the GTF file and to remove any nested
        quotes that may cause issues with other programs or tools that parse
        GTF files downstream of this script.
    """
    tags = {}  # store key, value pairs in 9th column
    metadata = re.split('; ', replace_nested_quotes(url_escape_inside_quotes(linelist[8].rstrip(';'))))
    for field in metadata:
        k,v = field.split(' ', 1)
        tags[k] = v.strip('"').strip("'")
    return tags


def default(v, d):
    """Returns d when v is empty (i.e an empty string, list, etc) or false.
    @param v <any>:
        The value to check if it is empty or false. This can be any type
        of value, such as a string, list, dictionary, or boolean. This is
        used to ensure that a default value is returned when the value is
        empty or false.
    @param d <any>:
        The default value to return if v is empty or false. This can be
        any type of value, such as a string, list, or dictionary.
    @return v <any>:
        Returns v if it is not empty or false, otherwise returns d.
    """
    if not v:
        v = d
    return v


def biotypes(gtf):
    """Creates dictionary to map each gene to its biotype. Biotype features
    listed as mRNA will be converted to protein_coding.
    @param gtf <str>:
        The path to the GTF file to parse for gene biotypes. This should be a
        GTF file that has been pre-processed by AGAT or another tool to ensure
        that it is in the correct format and contains the necessary metadata.
    @return gene2type <dict>:
        Returns a dictionary mapping gene IDs to their biotypes where the keys
        are gene IDs and the values are the biotypes. If a gene does not have a
        biotype, it will be set to "unknown". This is used to ensure that all
        genes in the GTF file have a biotype associated with them, which is
        required for some downstream applications.
    """
    gene2type = {}
    open_func = gzip.open if gtf.endswith('.gz') else open
    with open_func(gtf, 'rt') as file:
        for line in file:
            if line.startswith('#'):
                # Skip over comments in header section
                continue
            linelist = line.strip().split('\t')
            metadata = parse(linelist)
            # Get gene and biotype
            gene = lookup('gene_id', metadata)
            # Setting biotype to unknown as default
            # value, then checking if metadata contains
            # any fields with biotype as a sub string,
            # then if biotype is not an empty string
            # set it to whatever is in the gtf file
            if gene not in gene2type:
                gene2type[gene] = "unknown"
            biotype = contains('biotype', metadata)
            if default(biotype, 'unknown') != 'unknown':
                if biotype.lower() == 'mrna':
                    # agat_convert_sp_gff2gtf.pl does
                    # not set this value correct even
                    # when it is in the original GTF
                    # file, fixing the problem for
                    # RSeQC TIN reference file
                    biotype = 'protein_coding'
                gene2type[gene] = biotype
    return gene2type


def formatted(metadata):
    """Reformats key, value metadata to be written into the 9th column.
    @param metadata <dict>:
        The dictionary containing key, value pairs parsed from the 9th column
        of a GTF file. This should be the output of the parse() function.
    @return out <str>:
        Returns a formatted string of key, value pairs in the format of
        'key "value";'. This is used to ensure that the metadata is in the
        correct format for the 9th column of a GTF file and to remove any
        nested quotes that may cause issues with other programs or tools
        that parse GTF files downstream of this script.
    """
    out = ''
    for k,v in metadata.items():
        out += '{} "{}"; '.format(k,v)
    out = out.rstrip(' ')
    return out


def main():
    """Main function to run the script."""
    # Parse command line arguments
    args = parse_cli_arguments()
    input_gtf = args.input_gtf
    output_gtf = args.output_gtf
    gene_name_attribute = args.gene_name_attribute
    log("Running {0} script with the following options: ".format(_NAME), args)
    log("Parsing biotypes from input GTF file: ", input_gtf)
    g2b = biotypes(input_gtf)

    # Create output directory if
    # it does not exist
    output_dir = os.path.abspath(os.path.dirname(args.output_gtf))
    if not os.path.exists(output_dir):
        try: os.makedirs(output_dir)
        except OSError as e:
            fatal(
                "Fatal error: Failed to create output directory: {0}\n{1}".format(
                    output_dir, e
                )
            )

    # Parse the input GTF file and write cleaned
    # output to the output GTF file
    log("Started checking input GTF file: ", input_gtf)
    ofh = open(output_gtf, 'w')
    # Handler for opening files, i.e.
    # uncompressed or gzip files
    open_func = gzip.open if input_gtf.endswith('.gz') else open
    line_number = 0  # Used for error reporting
    with open_func(input_gtf, 'rt') as file:
        for line in file:
            line_number += 1
            if line.startswith('#'):
                # Skip over comments in header section
                ofh.write(line.strip()+'\n')
                continue
            linelist = line.strip().split('\t')
            feature = linelist[2]
            metadata = parse(linelist)
            # Should always be in GTF file
            gene_id = lookup('gene_id', metadata)
            if feature == 'gene':
                # May not be in GTF, add as needed
                gene_name = default(lookup(gene_name_attribute, metadata) , gene_id)
                metadata['gene_name'] = gene_name
                gene_biotype = default(lookup('gene_biotype', metadata) , g2b[gene_id])
                metadata['gene_biotype'] = gene_biotype
            elif feature in ['transcript', 'exon']:
                # May not be in GTF, add as needed
                # assumes transcript_id is in gtf
                gene_name = default(lookup(gene_name_attribute, metadata) , gene_id)
                metadata['gene_name'] = gene_name
                gene_biotype = default(lookup('gene_biotype', metadata) , g2b[gene_id])
                metadata['gene_biotype'] = gene_biotype
                transcript_id = lookup('transcript_id', metadata)
                transcript_name = default(lookup('transcript_name', metadata) , transcript_id)
                metadata['transcript_name'] = transcript_name
                transcript_type = default(lookup('transcript_type', metadata) , g2b[gene_id])
                metadata['transcript_type'] = transcript_type
            tags = formatted(metadata)
            linelist[8] = tags
            ofh.write("\t".join(linelist)+'\n')
    ofh.close()
    log("Finished writing output GTF file (contains {0} lines): {1}".format(line_number, output_gtf))


if __name__ == '__main__':
    main()
