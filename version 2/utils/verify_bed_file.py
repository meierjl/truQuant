class InvalidBedFormatException(Exception):
    def __init__(self, message, file, line):
        super().__init__(message)
        self.message = message
        self.file = file
        self.line = line

    def __str__(self):
        return "InvalidBedFormatException: " + self.message + " " + self.file + ".\n" + "Line: " + self.line.rstrip()


def __verify_file(filename):
    with open(filename) as file:
        first_line = file.readlines()[0]

    if len(first_line.split()) != 6:
        raise InvalidBedFormatException("Not all columns found in ", filename, first_line)

    chromosome, left, right, name, score, strand = first_line.split()

    try:
        left = int(left)
    except ValueError as e:
        raise InvalidBedFormatException("Left value is not an integer in ", filename, first_line)

    try:
        right = int(right)
    except ValueError as e:
        raise InvalidBedFormatException("Right value is not an integer in ", filename, first_line)

    if strand != "+" and strand != "-":
        raise InvalidBedFormatException("Strand is not + or - in ", filename, first_line)


def verify_bed_files(*filenames):
    """
        Raises a InvalidBedFormatException error if one of the files is not in bed format (chromosome, left, right,
        name, score, strand)

        :param files: files to be deleted
        :type files: str or list
        """
    for filename in filenames:
        if isinstance(filename, list):
            for sub_filename in filename:
                __verify_file(sub_filename)
        else:
            __verify_file(filename)
