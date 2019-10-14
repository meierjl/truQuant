'''
truQuant_quantitation
October 8, 2019
Geoff Collins
'''

from progress_bar import printProgressBar
import csv
import getopt
import sys
import os

five_prime_counts_dict = {}
three_prime_counts_dict = {}

def make_5_prime_counts_dict(sequencing_filename):
    with open(sequencing_filename, 'r') as file:
        num_of_lines = sum(1 for _ in file)
        file.seek(0)
        print("Building the 5' and 3' counts dicts.")
        printProgressBar(0, num_of_lines, prefix='Progress:', suffix='Complete', length=50)

        for i, line in enumerate(file):
            chromosome, left, right, _, _, strand = line.split()

            if strand == "+":
                five_prime_position = left
                three_prime_position = right
            else:
                five_prime_position = right
                three_prime_position = left

            if chromosome + '_' + five_prime_position + '_' + strand in five_prime_counts_dict:
                # If that base has a 5' end mapped already, add one to the value
                five_prime_counts_dict[chromosome + '_' + five_prime_position + '_' + strand] += 1
            else:
                # if that base does not have a mapped 5' read yet, set the value to 1
                five_prime_counts_dict[chromosome + '_' + five_prime_position + '_' + strand] = 1

            if chromosome + '_' + three_prime_position + '_' + strand in three_prime_counts_dict:
                # If that base has a 3' end mapped already, add one to the value
                three_prime_counts_dict[chromosome + '_' + three_prime_position + '_' + strand] += 1
            else:
                # if that base does not have a mapped 3' read yet, set the value to 1
                three_prime_counts_dict[chromosome + '_' + three_prime_position + '_' + strand] = 1

            if i % 100_000 == 0 or i == num_of_lines - 1:
                printProgressBar(i + 1, num_of_lines, prefix='Progress:', suffix='Complete', length=50)


def blacklist_regions(blacklisting_file):
    with open(blacklisting_file, 'r') as file:
        num_of_lines = sum(1 for _ in file)
        file.seek(0)
        print("Blacklisting.")
        printProgressBar(0, num_of_lines, prefix='Progress:', suffix='Complete', length=50)
        for curr_line_num, line in enumerate(file):
            chromosome, left, right, name, _, strand = line.split()

            # Loop through each base in the blacklist region and change the number of counts in the 3' dictionary to zero
            for i in range(int(left), int(right)):
                if chromosome + '_' + str(i) + '_' + strand in three_prime_counts_dict:
                    # If that base has a 3' end mapped already, add one to the value
                    three_prime_counts_dict[chromosome + '_' + str(i) + '_' + strand] = 0

            if curr_line_num % 10_000 == 0 or curr_line_num == num_of_lines - 1:
                printProgressBar(curr_line_num + 1, num_of_lines, prefix='Progress:', suffix='Complete', length=50)


def get_counts_in_paused_region(regions_filename, paused_region_output_filename):
    with open(regions_filename, 'r') as file:
        with open(paused_region_output_filename, 'w') as output_file:
            tsv_writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

            num_of_lines = sum(1 for _ in file)
            file.seek(0)
            print("Getting paused region counts.")
            printProgressBar(0, num_of_lines, prefix='Progress:', suffix='Complete', length=50)

            for curr_line_num, line in enumerate(file):
                chromosome, left, right, gene_name, _, strand = line.split()

                curr_total = 0
                # Loop through each base in the region
                for i in range(int(left), int(right)+1):
                    # Get the height and add it to the curr_total
                    if chromosome + '_' + str(i) + '_' + strand in five_prime_counts_dict:
                        curr_total += five_prime_counts_dict[chromosome + '_' + str(i) + '_' + strand]

                tsv_writer.writerow([gene_name, curr_total])

                if curr_line_num % 300 == 0 or curr_line_num == num_of_lines - 1:
                    printProgressBar(curr_line_num + 1, num_of_lines, prefix='Progress:', suffix='Complete', length=50)


def get_counts_in_gene_bodies(regions_filename, gene_bodies_output_filename):
    with open(regions_filename, 'r') as file:
        with open(gene_bodies_output_filename, 'w') as output_file:
            tsv_writer = csv.writer(output_file, delimiter='\t', lineterminator='\n')

            num_of_lines = sum(1 for _ in file)
            file.seek(0)
            print("Getting gene body counts.")
            printProgressBar(0, num_of_lines, prefix='Progress:', suffix='Complete', length=50)

            for curr_line_num, line in enumerate(file):
                chromosome, left, right, gene_name, _, strand = line.split()

                curr_total = 0
                # Loop through each base in the region
                for i in range(int(left), int(right)+1):
                    # Get the height and add it to the curr_total
                    if chromosome + '_' + str(i) + '_' + strand in three_prime_counts_dict:
                        curr_total += three_prime_counts_dict[chromosome + '_' + str(i) + '_' + strand]

                tsv_writer.writerow([gene_name, curr_total])
                printProgressBar(curr_line_num + 1, num_of_lines, prefix='Progress:', suffix='Complete', length=50)

if __name__ == '__main__':

    arguments = sys.argv[1:]

    paused_region_file, gene_body_region_file, sequencing_file, blacklisting_file = ['']*4

    try:
        opts, args = getopt.getopt(arguments, "p:g:s:b:")
        # P = Paused region file
        # G = Gene body region file
        # S = Sequencing data file
        # B = Blacklist file

    except getopt.GetoptError as error:
        print(error)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-p':
            paused_region_file = arg
        elif opt == '-g':
            gene_body_region_file = arg
        elif opt == '-s':
            sequencing_file = arg
        elif opt == '-b':
            blacklisting_file = arg
        else:
            print("One of the arguments you have entered is not valid.")

    if '' in [paused_region_file, gene_body_region_file, sequencing_file, blacklisting_file]:
        print("You have not provided all of the necessary files. Exiting.")
        sys.exit(2)

    paused_region_output_filename = os.path.basename(sequencing_file).replace('.bed', '') + "_paused_counts.bed"
    gene_bodies_output_filename = os.path.basename(sequencing_file).replace('.bed', '') + "_gene_bodies_counts.bed"

    make_5_prime_counts_dict(sequencing_file)
    blacklist_regions(blacklisting_file)
    get_counts_in_paused_region(paused_region_file, paused_region_output_filename)
    get_counts_in_gene_bodies(gene_body_region_file, gene_bodies_output_filename)
