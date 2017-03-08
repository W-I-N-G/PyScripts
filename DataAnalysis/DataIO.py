"""!
@file DataAnalysis/DataIO.py
@package DataAnalysis

@defgroup DataIO DataIO

@brief Routines support data input and output for the Bevins PyScripts package.

@author James Bevins

@date 7Mar17
"""

import re

#------------------------------------------------------------------------------#
def readDelimitedDataFile(path, delimiter=" +", header=0, breakText=""):
    """!
    @ingroup DataIO
    Reads in a column major dataset and and returns the results.

    @param path: \e string \n
        Absolute path to the file \n
    @param delimiter: \e string \n
        Indicator for the character used to separate columns.  Uses the
        re.split delimiter definitions. \n
    @param header: \e integer \n
        The number of header lines to skip \n
    @param breakText: \e string \n
        Text indicating the end of the data to be read. \n

    @return <em> list of lists <\em>: A nxm list where n is the number of
        columns of data in the input file and m is the number of rows \n
    """

    assert header >= 0, "Valid specifications for the number of header lines \
                      to skip must be positive."
    assert type(delimiter) == str, "Valid specifications for the delimiters \
                                 must be strings."

    data = []

    # Open file
    try:
        f = open(path, 'r')
        # Skip n header lines
        for i in range(0, header):
            line = f.next()

        # Read the file line by line and store the values
        for line in f:
            if line.rstrip() == breakText:
                break
            splitList = re.split(delimiter, line.strip())
            for i in range(0, len(splitList)):
                if len(data) < i+1:
                    data.append([])
                data[i].append(float(splitList[i]))

        # Close the file
        f.close()
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)

    return data
