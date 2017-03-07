"""!
@file Unfolding/Root.py
@package Unfolding

@defgroup Root Root

@brief Routines to support unfolding and python to root interfaces.

@author James Bevins

@date 6Mar17
"""

#------------------------------------------------------------------------------#
class calibParams(object):
    """!
    @ingroup Root
    This class stores the calibration parameters used for energy claibration
    and gausssian smearing of response matrices for root and unfolding
    algorithms.
    """


    ##
    def __init__(self, path=''):
        """!
        Constructor to build the calibParams class. If path is specified, the
        object attributes are populated.

        @param self: <em> object pointer </em> \n
            The object pointer. \n
        @param path: \e string \n
            The path to the Root calibration obput .txt file \n
        """

        if path != '':
            params = self.readParams(path)
            ## @var a: \e float
            # The linear energy calibration parameter. \f$E=a*Ch+b\f$
            self.a = params[0]

            ## @var b: \e float
            # The offset energy calibration parameter. \f$E=a*Ch+b\f$
            self.b = params[1]

            ## @var alpha: \e float
            # The first resolution function parameter. \f$\Gamma(e)=
            # \sqrt{\alpha^2 + \frac{\beta^2}{E} + \frac{\gamma^2}{E^2}}\f$
            self.alpha = params[2]

            ## @var beta: \e float
            # The second resolution function parameter. \f$\Gamma(e)=
            # \sqrt{\alpha^2 + \frac{\beta^2}{E} + \frac{\gamma^2}{E^2}}\f$
            self.beta = params[3]

            ## @var gamma: \e float
            # The third resolution function parameter. \f$\Gamma(e)=
            # \sqrt{\alpha^2 + \frac{\beta^2}{E} + \frac{\gamma^2}{E^2}}\f$
            self.gamma = params[4]

    def __repr__(self):
        """!
        Object print function.

        @param self: <em> object pointer </em> \n
            The object pointer. \n
        """

        return "ADVANTG Settings({}, {}, {}, {}, {})".format(self.a, self.n,
                                                             self.alpha,
                                                             self.beta,
                                                             self.gamma)

    def __str__(self):
        """!
        Human readable object print function.

        @param self: <em> object pointer </em> \n
            The object pointer. \n
        """
        header = ["\nCalibration Parameters"]
        header += ["a = {}".format(self.a)]
        header += ["b = {}".format(self.b)]
        header += ["alpha = {}".format(self.alpha)]
        header += ["beta = {}".format(self.beta)]
        header += ["gamma = {}".format(self.gamma)]
        header = "\n".join(header)+"\n"
        return header

    def readParams(self, path):
        """!
        Reads in the calibration parameters from a Root calibration formatted
        text file.

        @param path: \e string \n
            The path to the Root calibration obput .txt file \n

        @return <em> list of floats </em>
            A list containing [a, b, \f$\alpha\f$, \f$\beta\f$, \f$\gamma\f$]
        """

        params = []

        # Open file
        try:
            f = open(path, 'r')

            # Skip 1 header lines
            f.next()

            # Read the file line by line and store the values
            for i in range(5):
                params.append(f.next().rstrip())

            # Close the file
            f.close()
        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)

        print params
        return params
