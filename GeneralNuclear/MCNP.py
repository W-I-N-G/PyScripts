"""!
@file GeneralNuclear/MCNP.py
@package GeneralNuclear

@defgroup MCNP MCNP

@brief Routines and functions to interact with MCNP input and output.

@author James Bevins

@date 16April17
"""

import pandas as pd

#------------------------------------------------------------------------------#
def read_tally(path, tallyNum, readGroups=True, splitTally=False):
    """!
    @ingroup MCNP
    Reads a specified tally from a standard MCNP output and returns the tally
    results in the form of a data frame (if group-wise) and/or total value.

    @param path: \e string \n
        The path, including filename, to the MCNP output file to be read \n
    @param tallyNum: <em> string, integer, or float </em> \n
        The number of the tally to be read \n
    @param readGroups: \e boolean \n
        Optional specifier to return group-wise tally information. If false,
        only the total is returned. If true, but no group-wise information is
        present, only the total is returned and a warning is thrown. \n
    @param splitTally: \e boolean \n
        Optional specifier indicating that the tally is split into multiple
        sub tallies.  This is the case when multiple group-wise specifiers
        (i.e. energy and angle) are used.  This specifier is only valid if
        readGroups==True. If used. this returns a dictionary of pandas
        dataframes where the value is a tuple of the sub-bin specifier
        and a corresponding dataframe. \n
    @return <em> pandas dataframe: </em> Pandas dataframe containing the
        results and uncertainties. Only returned if readGroups==True and
        splitTally==False. The columns are ['bin', 'tally', 'uncertainty']. \n
    @return <em> dictionary of tuples (string, pandas dataframe, total,
        uncertainty): </em> A dictionary of pandas dataframes where the value
        is a tuple of the sub-bin specifier, a corresponding dataframe, the
        sub-tally total, and the sub-tally uncertainty. Only used if
        readGroups==True and splitTally==True. The dataframe columns are
        ['bin', 'tally', 'uncertainty']. \n
    @return \e float: The tally total. \n
    @return \e float: The tally uncertainty. \n
    """

    assert isinstance(path, str) == True, 'Path must be of type str.'

    # Convert tallyNum type
    if type(tallyNum) == int or type(tallyNum) == float:
        tallyNum = str(tallyNum)
    else:
        print "ERROR: tallyNum type unknown."
        return

    # Initialize data structures
    tally = False
    colNames = ['bin', 'tally', 'uncertainty']
    if readGroups == True and splitTally == False:
        df = pd.DataFrame(columns=colNames)
    elif readGroups == True and splitTally == True:
        tallyDict = {}
        df = pd.DataFrame(columns=colNames)

    # Determine number of header lines for tally
    if float(tallyNum) -round(float(tallyNum), -1) == 1:
        headerLines = 6
    if float(tallyNum) - round(float(tallyNum), -1) == 4:
        headerLines = 10

    # Create and open input file
    try:
        with open(path, "r") as f:

            # Read the output file line by line
            for line in f:

                # Find key word for start of flux array
                splitList = line.strip().split()
                if len(splitList) >= 3:
                    if splitList[0].strip() == "1tally" and \
                          splitList[1].strip() == tallyNum and \
                          splitList[2].strip() == "nps":
                        tally = True

                        # Skip header lines
                        for i in range(0, headerLines):
                            line = f.next().strip()
                        if readGroups == True and splitTally == True:
                            subTallyName = line
                            line = f.next().strip()
                        splitList = f.next().split()

                # Fill data structure
                if tally == True:
                    # If blank line, skip
                    if len(splitList) != 0:
                        # If end of tally found, consolidate gains
                        if splitList[0].strip() == "total":
                            if readGroups == True and splitTally == True:
                                tallyDict[len(tallyDict)] = (subTallyName, df,
                                                           float(splitList[1]),
                                                           float(splitList[2]))
                                df = pd.DataFrame(columns=colNames)
                                for i in range(0, 3):
                                    splitList = f.next().split()
                                if splitList[0][0:4].strip() == '====':
                                    return tallyDict
                                else:
                                    subTallyName = " ".join(splitList)
                                    splitList = f.next().split()

                            else:
                                total = splitList[1]
                                uncert = splitList[2]
                        # Defintely nothing left to gain, return
                        elif splitList[0][0:4].strip() == '====':
                            if readGroups == True and splitTally == False:
                                return df, total, uncert
                            else:
                                return total, uncert

                        # Store group-wise results
                        elif readGroups == True:
                            df = df.append(pd.Series(
                                                  [float(splitList[0].strip()),
                                                   float(splitList[1].strip()),
                                                   float(splitList[2].strip())],
                                             index=colNames), ignore_index=True)

        # Close the file
        f.close()

    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)
        print "File not found was: {0}".format(path)

    # Test that the file closed
    assert f.closed == True, "File ({}) did not close properly.".format(path)
