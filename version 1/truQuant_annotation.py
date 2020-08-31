'''
truQuant_annotation
October 8, 2019
Geoff Collins
'''

import csv
import getopt
import sys
import os

from progress_bar import printProgressBar

regions_dict = {}
gene_tsr_dict = {}
max_tsrs_dict = {}
non_max_tsrs_dict = {}
five_prime_counts_dict = {}

annotations = []

def make_regions_list(regions_filename):
    with open(regions_filename, 'r') as file:
        for line in file:
            chromosome, left, right, strand, gene_name, met_left, met_right = line.split()

            if strand == "+":
                # We want to go from the left position to the met_left position
                if chromosome in regions_dict:
                    regions_dict[chromosome].append([chromosome, int(left) - annotation_extension, met_left, strand, gene_name])
                else:
                    regions_dict[chromosome] = [[chromosome, int(left) - annotation_extension, met_left, strand, gene_name]]
            else:
                if chromosome in regions_dict:
                    regions_dict[chromosome].append([chromosome, met_right, int(right) + annotation_extension, strand, gene_name])
                else:
                    regions_dict[chromosome] = [[chromosome, met_right, int(right) + annotation_extension, strand, gene_name]]

            annotations.append([chromosome, left, right, strand, gene_name])


# Loop through the mapped sequencing data file to load the five_prime_counts_dict
def make_5_prime_counts_dict(sequencing_filename):
    with open(sequencing_filename, 'r') as file:
        num_of_lines = sum(1 for _ in file)
        file.seek(0)
        print("Building the 5' counts dict.")
        printProgressBar(0, num_of_lines, prefix='Progress:', suffix='Complete', length=50)

        for i, line in enumerate(file):
            chromosome, left, right, _, _, strand = line.split()

            if strand == "+":
                five_prime_position = left
            else:
                five_prime_position = right

            if chromosome + '_' + five_prime_position + '_' + strand in five_prime_counts_dict:
                # If that base has a 5' end mapped already, add one to the value
                five_prime_counts_dict[chromosome + '_' + five_prime_position + '_' + strand] += 1
            else:
                # If that base does not have a mapped 5' read yet, set the value to 1
                five_prime_counts_dict[chromosome + '_' + five_prime_position + '_' + strand] = 1

            if i % 100_000 == 0 or i == num_of_lines - 1:
                printProgressBar(i + 1, num_of_lines, prefix='Progress:', suffix='Complete', length=50)


def map_tsrs(tsr_filename):
    with open(tsr_filename, 'r') as file:

        num_of_lines = sum(1 for _ in file)
        file.seek(0)
        print("Mapping TSRs.")
        printProgressBar(0, num_of_lines, prefix='Progress:', suffix='Complete', length=50)

        # Loop through the TSR file
        for i, line in enumerate(file):
            tsr_chromosome, tsr_left, tsr_right, _, tsr_counts, tsr_strand = line.split()
            # Set the correct chromosome to true. It is false when the chromosome in the TSR file is not present in the annotation file
            correct_chromosome = True

            try:
                regions_dict[tsr_chromosome]
            except KeyError:
                # This means the TSR file has a chromosome that is not present in the annotation file
                correct_chromosome = False

            if correct_chromosome:
                for region in regions_dict[tsr_chromosome]:
                    chromosome, left, right, strand, gene_name = region

                    # If there is a TSR as that base pair, add it to the gene_tsr dict
                    if tsr_strand == strand and not (int(tsr_right) < int(left) or int(tsr_left) > int(right)):
                        if gene_name not in gene_tsr_dict:
                            gene_tsr_dict[gene_name] = [line.split()]
                        else:
                            gene_tsr_dict[gene_name].append(line.split())

                    # If the left of the region is past the right side of the TSR, we don't need to search the rest of the file
                    if int(left) > int(tsr_right):
                        break

            if i % 300 == 0 or i == num_of_lines - 1:
                printProgressBar(i + 1, num_of_lines, prefix='Progress:', suffix='Complete', length=50)


def find_max_tsr():
    for gene_name in gene_tsr_dict:
        for tsr in gene_tsr_dict[gene_name]:
            tsr_chromosome, tsr_left, tsr_right, _, tsr_counts, tsr_strand = tsr

            if gene_name not in max_tsrs_dict:
                max_tsrs_dict[gene_name] = tsr
            elif int(tsr_counts) > int(max_tsrs_dict[gene_name][4]):  # If the current max is larger than the previous max
                # Add the old tsr to the non_max
                if gene_name not in non_max_tsrs_dict:
                    non_max_tsrs_dict[gene_name] = [max_tsrs_dict[gene_name]]
                else:
                    non_max_tsrs_dict[gene_name].append(max_tsrs_dict[gene_name])

                # Put in the new one
                max_tsrs_dict[gene_name] = tsr
            else:
                if gene_name not in non_max_tsrs_dict:
                    non_max_tsrs_dict[gene_name] = [tsr]
                else:
                    non_max_tsrs_dict[gene_name].append(tsr)


def rewrite_annotations(rewritten_annotations_filename):
    # Take the original annotation and change the 5' end to the max TSS in the max TSR
    with open(rewritten_annotations_filename, 'w') as output_file:
        tsv_writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

        for annotation in annotations:
            chromosome, left, right, strand, gene_name = annotation

            # If there is a max TSR mapped to the gene
            if gene_name in max_tsrs_dict:
                if strand == "+":
                    # We need to loop through the TSR to find the max TSS. Then set that to the left position
                    max_pos = 0
                    curr_max = -1

                    for i in range(tsr_size):
                        position = max_tsrs_dict[gene_name][0] + "_" + str(int(max_tsrs_dict[gene_name][1]) + i) + "_" + \
                                   max_tsrs_dict[gene_name][5]

                        if position in five_prime_counts_dict:
                            if five_prime_counts_dict[position] > curr_max:
                                curr_max = five_prime_counts_dict[position]
                                max_pos = i

                    left = int(max_tsrs_dict[gene_name][1]) + max_pos
                else:

                    max_pos = 0
                    curr_max = -1

                    for i in range(tsr_size):
                        position = max_tsrs_dict[gene_name][0] + "_" + str(int(max_tsrs_dict[gene_name][1]) + i) + "_" + \
                                   max_tsrs_dict[gene_name][5]

                        if position in five_prime_counts_dict:
                            if five_prime_counts_dict[position] > curr_max:
                                curr_max = five_prime_counts_dict[position]
                                max_pos = i

                    right = int(max_tsrs_dict[gene_name][1]) + max_pos

                tsv_writer.writerow([chromosome, left, right, gene_name, curr_max, strand])


def define_gene_bodies_and_paused_regions(rewritten_annotations_filename, paused_region_filename, gene_body_region_filename):
    with open(rewritten_annotations_filename, 'r') as file:
        with open(paused_region_filename, 'w') as paused_file:
            paused_writer = csv.writer(paused_file, delimiter='\t', lineterminator='\n')
            with open(gene_body_region_filename, 'w') as gene_bodies_file:
                gene_bodies_writer = csv.writer(gene_bodies_file, delimiter='\t', lineterminator='\n')

                for line in file:
                    chromosome, left, right, gene_name, _, strand = line.split()

                    if strand == "+":
                        paused_writer.writerow([chromosome, int(left) - paused_region_radius,
                                                int(left) + paused_region_radius + 1, gene_name, 0, strand])

                        gene_bodies_writer.writerow([chromosome, int(left) + paused_region_radius + 2, int(right),
                                                     gene_name, 0, strand])
                    else:
                        paused_writer.writerow([chromosome, int(right) - (paused_region_radius + 1),
                                                int(right) + paused_region_radius, gene_name, 0, strand])

                        gene_bodies_writer.writerow([chromosome, int(left), int(right) - (paused_region_radius + 2),
                                                     gene_name, 0, strand])


def make_blacklisted_regions(blacklist_filename):
    # We are blacklisting all non max TSRs that are not inside the paused region
    # I expand the TSR 75 bp on each side to give some bore buffer for the paused region
    with open(blacklist_filename, 'w') as output_file:
        tsv_writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

        for gene_name in non_max_tsrs_dict:
            for tsr in non_max_tsrs_dict[gene_name]:
                tsr_chromosome, tsr_left, tsr_right, _, tsr_counts, tsr_strand = tsr

                if tsr_strand == "+":
                    if int(tsr_right) + blacklist_extension < int(max_tsrs_dict[gene_name][1]) - blacklist_extension or \
                            int(tsr_left) - blacklist_extension > int(max_tsrs_dict[gene_name][2]) + blacklist_extension + 1:
                        # If the expanded TSR has no overlap with the paused region

                        if int(tsr_counts) >= int(max_tsrs_dict[gene_name][4]) * percent_for_blacklisting:
                            # If the current TSR is the max seen so far
                            tsv_writer.writerow([tsr_chromosome, int(tsr_left) - blacklist_extension,
                                                 int(tsr_right) + blacklist_extension + 1, "blacklist", 0, tsr_strand])

                else:
                    if int(tsr_right) + blacklist_extension + 1 < int(max_tsrs_dict[gene_name][1]) - (blacklist_extension + 1) \
                            or int(tsr_left) - blacklist_extension > int(max_tsrs_dict[gene_name][2]) + blacklist_extension + 1:
                        # If the expanded TSR has no overlap with the paused region

                        if int(tsr_counts) >= int(max_tsrs_dict[gene_name][4]) * percent_for_blacklisting:
                            tsv_writer.writerow([tsr_chromosome, int(tsr_left) - (blacklist_extension + 1),
                                                 int(tsr_right) + blacklist_extension, "blacklist", 0, tsr_strand])


if __name__ == "__main__":
    arguments = sys.argv[1:]

    annotation_extension = 1000
    paused_region_radius = 75
    blacklist_extension = 75
    tsr_size = 20
    percent_for_blacklisting = 0.3

    tsr_file, sequencing_file, annotation_file = [''] * 3

    try:
        opts, args = getopt.getopt(arguments, "t:s:a:r:n:e:b:p:")
        # T = TSR file
        # S = Sequencing file .bed
        # A = Annotation file
        # N = TSR window size
        # E = Annotation Extension
        # B = Blacklist extension amount
        # P = Minimum percentage for blacklisting TSRs

    except getopt.GetoptError as error:
        print(error)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-t':
            tsr_file = arg
        elif opt == '-s':
            sequencing_file = arg
        elif opt == '-a':
            annotation_file = arg
        elif opt == '-r':
            paused_region_radius = int(arg)
        elif opt == '-n':
            tsr_size = int(arg)
        elif opt == '-e':
            annotation_extension = int(arg)
        elif opt == '-b':
            blacklist_extension = int(arg)
        elif opt == '-p':
            percent_for_blacklisting = float(arg)
        else:
            print("One of the arguments you have entered is not valid.")

    if '' in [tsr_file, sequencing_file, annotation_file]:
        print("You have not provided all of the necessary files. Exiting.")
        sys.exit(2)

    parameters_string = "-r" + str(paused_region_radius) + "n" + str(tsr_size) + "e" + str(annotation_extension) + "b" +\
                        str(blacklist_extension) + "p" + str(percent_for_blacklisting) + "-"

    rewritten_annotations_filename = os.path.basename(sequencing_file).replace('.bed', '') + parameters_string + "TruQuant_annotations.bed"
    paused_region_filename = os.path.basename(sequencing_file).replace('.bed', '') + parameters_string + "paused_regions.bed"
    gene_body_region_filename = os.path.basename(sequencing_file).replace('.bed', '') + parameters_string + "gene_body_regions.bed"
    blacklist_filename = os.path.basename(sequencing_file).replace('.bed', '') + parameters_string + "blacklisted_regions.bed"

    # 1: Make the regions we are going to be searching for max TSSs in max TSRs

    make_regions_list(annotation_file)

    # 2: Define the 5' end of the genes

    make_5_prime_counts_dict(sequencing_file)
    map_tsrs(tsr_file)
    find_max_tsr()
    rewrite_annotations(rewritten_annotations_filename)

    # 3: Make the regions for gene bodies and paused regions. Make the TSR blacklist

    define_gene_bodies_and_paused_regions(rewritten_annotations_filename, paused_region_filename, gene_body_region_filename)
    make_blacklisted_regions(blacklist_filename)
