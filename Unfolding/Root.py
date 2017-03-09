"""!
@file Unfolding/Root.py
@package Unfolding

@defgroup Root Root

@brief Routines to support unfolding and python to root interfaces.

@author James Bevins

@date 8Mar17
"""

import os
import sys

from datetime import datetime


sys.path.insert(0, os.path.abspath(
 '/home/pyne-user/Dropbox/UCB/Computational_Tools/Scripts/Python/GeneralNuclear'))
from BasicNuclearCalcs import solid_angle_approx

#------------------------------------------------------------------------------#
class CalibParams(object):
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
            The path to the Root calibration output .txt file \n
        """

        if path != '':
            params = self.read_params(path)
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

    def read_params(self, path):
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

        return params


#------------------------------------------------------------------------------#
class FluxNormalization(object):
    """!
    @ingroup Root
    This class does the calculations necessary to normalize a measures spectra
    and obtain absolute results.
    """

    ##
    def __init__(self, path='', **kwargs):
        """!
        Constructor to build the calibParams class. If path is specified, the
        object attributes are populated.

        @param self: <em> fluxNormalization pointer </em> \n
            The object pointer. \n
        @param path: \e string \n
            The path to the current monitor output \n
        @param path: \e kwargs \n
            Optional variables for the setCurrentMonitor function \n
        """

        if path != '':
            self.set_current_monitor(path, **kwargs)
        else:
            ## @var currentMonitor: \e float
            # The integrated current according to the current monitor in microA
            self.currentMonitor = 1

        ## @var currentIntegrator: \e float
        # The integrated current according to the current integrator in microA
        self.currentIntegrator = 1

        ## @var solidAngle: \e float
        # The solid angle subtended by the detector
        self.solidAngle = 1

    def __repr__(self):
        """!
        Object print function.

        @param self: <em> fluxNormalization pointer </em> \n
            The object pointer. \n
        """

        return "ADVANTG Settings({}, {}, {})".format(self.currentMonitor,
                                                     self.currentIntegrator,
                                                     self.solidAngle)

    def __str__(self):
        """!
        Human readable object print function.

        @param self: <em> fluxNormalization pointer </em> \n
            The object pointer. \n
        """
        header = ["\nNormalization Parameters"]
        header += ["Current Monitor Inegrated Current = {} microA"\
                   .format(self.currentMonitor)]
        header += ["Current Integrator Reading = {} microA"\
                   .format(self.currentIntegrator)]
        header += ["Solid Angle = {} sr".format(self.solidAngle)]
        header = "\n".join(header)+"\n"
        return header

    def set_current_monitor(self, path, startTime=None, stopTime=None, scale=1):
        """!
        Reads in a current monitor file and set the currentMonitor attibute
        for the class.  The return value is in microA.

        @param self: <em> fluxNormalization pointer </em> \n
            The object pointer. \n
        @param path: \e string \n
            Absolute path to the file. \n
        @param startTime: <em> integer or float <\em> \n
            The run start time. Used when the file contains multiple runs. \n
        @param stopTime: <em> integer or float <\em> \n
            The run stop time. Used when the file contains multiple runs. \n
        @param scale: <em> integer or float <\em> \n
            Speciies the scaling to perform on the current readings. This is
            required if the readings are not in nA, and the scale factor is
            a ratio of nA to the scale used. \n
        """

        # Initialize variables
        self.currentMonitor = 0
        if startTime == None:
            startTime = datetime(1984, 8, 2, 12, 12, 0)
        if stopTime == None:
            stopTime = datetime(2084, 8, 2, 12, 12, 0)

        # Open file
        try:
            f = open(path, 'r')

            # Store first time slice
            line = f.next().rstrip().split('\t')
            prevTime = datetime.strptime(line[0]+line[1], '%m/%d/%Y%I:%M:%S %p')

            for line in f:
                line = line.rstrip().split('\t')
                curTime = datetime.strptime(line[0]+line[1],
                                            '%m/%d/%Y%I:%M:%S %p')

                # Reading occurred during the experiment and isn't noise
                if curTime > startTime and curTime < stopTime and \
                   abs(float(line[2])) > 0.000000001:
                    self.currentMonitor += float(line[2])\
                                        *(curTime-prevTime).total_seconds()
                prevTime = curTime

            # Close the file and set currentMonitor
            self.currentMonitor = self.currentMonitor*1E6/scale
            f.close()
        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)

    def set_solid_angle(self, dist, area):
        """!
        Sets the solid angle using the approximation of \f$ \frac{A}{d^2} \f$
        valid for d >> A.

        @param self: <em> fluxNormalization pointer </em> \n
            The object pointer. \n
        @param dist: <em> integer or float <\em> \n
            The distance from the source to detector in cm. \n
        @param area: <em> integer or float <\em> \n
            The cross-sectional area of the source in cm\f$^2\f$. \n
        """

        self.solidAngle = solid_angle_approx(dist, area)
