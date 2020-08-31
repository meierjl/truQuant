import os

from utils.make_random_filename import generate_random_filename


def run_coverage(regions_filename, sequencing_filename, output_filename='', flags=None):
    """
    Runs strand specific bedtools coverage to get the number of counts in the sequencing file in the regions file.

    :param regions_filename: filename of the regions of the genome to quantify
    :type regions_filename: str
    :param sequencing_filename: filename of the sequencing data collected
    :type sequencing_filename: str
    :param output_filename: optional name of the output file (will be random if not provided)
    :type output_filename: str
    :param flags: optional flags (like -d for depth at each base)
    :type flags: list
    :return: filename of the resultant bedtools coverage output
    :rtype: str
    """
    if output_filename == '':
        output_filename = generate_random_filename()

    # Convert the flags to string
    if flags != None:
        flag_string = ' '.join(flags)
    else:
        flag_string = ''

    os.system("bedtools coverage -s " + flag_string + " -a " + regions_filename + " -b " + sequencing_filename + " > " + output_filename)

    return output_filename
