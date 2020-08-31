from shutil import which


class ToolNotFoundException(Exception):
    def __init__(self, tool):
        super().__init__(tool)
        self.tool = tool

    def __str__(self):
        return self.tool + " is not installed."


def _cli_tool_exists(tool):
    """
    Checks if CLI tool exists. Throws ToolNotFoundException if tool does not exist.

    :param tool: name of the tool the program will use
    :type tool: str
    :exception ToolNotFoundException: if CLI tool does not exist
    """

    if which(tool) == None:
        raise ToolNotFoundException(tool)


def check_dependencies(*program_names):
    """
    Checks if the system has bedtools and tsrFinder
    :return:
    """
    # Check for bedtools and tsrFinder
    for program_name in program_names:
        if program_name == "bedtools":
            _cli_tool_exists("bedtools")

        if program_name == "tsrFinderM1I":
            _cli_tool_exists("tsrFinderM1I")

