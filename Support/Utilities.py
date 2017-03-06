"""!
@file Support/Utilities.py
@package Support

@defgroup Utilies Utilitites

@brief General support utilities.

@author James Bevins

@date 24Feb17
"""

import os

#------------------------------------------------------------------------------#
def pause():
    """!
    @ingroup Utilities

    Generic pause funtion to await user input before proceeding.

    @return none
    """

    try:
        input("Press enter to continue")
    except SyntaxError:
        pass
    
#------------------------------------------------------------------------------#
def remove_file(fileName):
    """!
    @ingroup Utilities
    Test for files existance; deletes that file if it exists.

    @param fileName: \e string \n
        The path to the file. \n

    @return none
    """

    try:
        os.remove(fileName)
    except OSError:
        print "WARNING: {} not found and cannot be deleted.".format(fileName)