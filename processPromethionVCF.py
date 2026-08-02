#!/usr/bin/env python


# example usage: python promethion_vcf.py --input <input_file> --output <output_file> 
import vcf
import re
import copy
import argparse

def main():
    # initiate the parser and handle arguments from command line
    args = parse_args()
    # Open file, this will also read in the header
    vcf_reader = vcf.Reader(open(args.input, 'r'))
    # Open file to write the filterd VCF
    vcf_writer = vcf.Writer(open(args.output, 'w'), vcf_reader)
    i = 0
    for record in vcf_reader:
        if "BND" in record.ID:
            i += i
            # store mate in another variable

            # copy.deepcopy() doesn't seem to generate a valid output, so we are doing this
            record_mate = copy.copy(record)
            record_mate.INFO = copy.copy(record.INFO)
            record_mate.ALT = copy.copy(record.ALT)
            record_mateid = record.ID + "B"

            # get chromosome and start position of mate
            variant_alt= str(record.ALT[0]).split(':')
            chrom = re.findall(r'\d+', variant_alt[0])
            start_pos = re.findall(r'\d+', variant_alt[1])
            # change main entries of mate
            record_mate.CHROM = chrom[0]
            record_mate.POS = start_pos[0]
            record_mate.ID = record_mateid

            # Set orientation of the mate of the mate breakpoint
            # Which side?
            alt_allele = record.ALT[0]
            if alt_allele.startswith("N"):
                # Is the inserted sequence flipped?
                if "[" in alt_allele:
                    comp_alt = "[" +  
            record_mate.ALT[0] = record.CHROM[0] + ":" + str(record.POS)

            # add MATEID to info column for both mates
            record_mate.INFO["MATEID"] = record.ID
            record.INFO["MATEID"] = record_mateid

            vcf_writer.write_record(record)
            vcf_writer.write_record(record_mate)
        else:
            vcf_writer.write_record(record)

    vcf_writer.close()


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input",
                        help="Input vcf file to be filtered",
                        required=True)
    parser.add_argument("--output",
                        help="Resulting vcf file after filtering",
                        required=True)

    args, unknown = parser.parse_known_args()

    return args


if __name__ == '__main__':
    main()
