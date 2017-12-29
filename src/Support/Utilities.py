"""!
@file Support/Utilities.py
@package Support

@defgroup Utilities Utilities

@brief General support utilities.

@author James Bevins

@date 20Sep17
"""

import os
import os.path

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
def removeFile(fileName):
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

#------------------------------------------------------------------------------#
class PDF(object):
    """!
    @ingroup Utilities
    Method to display PDF files imbedded into Jupyter notebooks.
    """
    def __init__(self, pdf, size=(200,200)):
        """!
        Constructor to build the PDF class.

        @param self: <em> object pointer </em> \n
            The object pointer. \n
        @param pdf: \e string \n
            The path to the pdf file to be displayed. In a multi-page pdf,
            the image on the second page can be displayed by using '<path>[1]'
            and so forth for other pages in the document. \n
        @param self: <em> tuple of integers </em> \n
            The integer value of the number of pixels to display for the
            image. \n
        """
        ## @var pdf: \e string
        # A path for the pdf image location.
        self.pdf = pdf

        ## @var size: \e string
        # The integer value of the number of pixels to display for the
        # image.
        self.size = size

    def _repr_html_(self):
        """!
        PDF HTML print function.

        @param self: <em> PDF pointer </em> \n
            The PDF pointer. \n
        """
        return '<iframe src={0} width={1[0]} height={1[1]}></iframe>'.format(
                                                         self.pdf, self.size)

    def _repr_latex_(self):
        """!
        PDF Latex print function.

        @param self: <em> PDF pointer </em> \n
            The PDF pointer. \n
        """
        return r'\includegraphics[width=1.0\textwidth]{{{0}}}'.format(
                                                         self.pdf)
#------------------------------------------------------------------------------#
def checkPath(path):
    """!
    @ingroup Utilities
    Test for files existance. Prints the results of the test to the screen.

    @param path: \e string \n
        The path to the file. \n

    @return none
    """
    if os.path.isfile(path): 
        print 'The file exists at: {}'.format(path)
    else:
        print 'ERROR: The file DOES NOT exist at: {}'.format(path)