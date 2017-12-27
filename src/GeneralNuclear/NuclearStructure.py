"""!
@file GeneralNuclear/NuclearStructure.py
@package GeneralNuclear

@defgroup NuclearStructure NuclearStructure

@brief Routines and functions to perform basic nuclear structure calculations

@author James Bevins

@date 22Sep17
"""

from math import ceil
#------------------------------------------------------------------------------#
class Atom(object):
    """!
    @ingroup NuclearStructure
    This class creates an object that stores basic nuclear properties of an
    atom.
    """

    ##
    def __init__(self, z, a, bindingEnergy=None, mass=None,
                 bindingEnergyPerA=None, massPerA=None):
        """!
        Constructor to build the Atom class.

        @param self: <em> object pointer </em> \n
            The object pointer. \n
        @param z: \e integer \n
            The number of protons in the atom. \n
        @param a: \e integer \n
            The number of nucleons in the atom. \n
        @param bindingEnergy: \e float \n
            The atom's total binding energy. \n
        @param mass: \e float \n
            The mass of the atom. \n
        @param bindingEnergyPerA: \e float \n
            The atom's binding energy per nucleon. \n
        @param mass: \e float \n
            The mass per nucleon of the atom. \n
        """
        assert z > 0 and z <=120, ("The number of protons must be "
                      "between 1 and 120 inclusive.")
        assert a >= z and a <=400, ("The number of nucleons must be "
                      "between 1 and 400 inclusive.")

        ## @var z: \e integer
        # The number of protons in the atom.
        self.z = int(z)
        ## @var a: \e integer
        # The number of nuclons in the atom.
        self.a = int(a)
        ## @var bindingEnergy: \e float
        # The atom's total binding energy.
        if bindingEnergy != None:
            self.bindingEnergy = bindingEnergy
        else:
            self.calcBE()
        ## @var mass: \e float
        # The mass of the atom.
        if mass != None:
            self.mass = mass
        else:
            self.calcSEMF()
        ## @var bindingEnergyPerA: \e float
        # The atom's total binding energy per nucleon.
        if bindingEnergyPerA != None:
            self.bindingEnergyPerA = bindingEnergyPerA
        ## @var massPerA: \e float
        # The mass of the atom per A.
        if massPerA != None:
            self.massPerA = massPerA

    def __repr__(self):
        """!
        Atom print function.

        @param self: <em> atom pointer </em> \n
            The atom pointer. \n
        """
        return "Atom({}, {}, {}, {}, {}, {})".format(self.z, self.a,
                                                 self.bindingEnergy,
                                                 self.mass,
                                                 self.bindingEnergyPerA,
                                                 self.massPerA)

    def __str__(self):
        """!
        Human readable Atom print function.

        @param self: <em> atom pointer </em> \n
            The atom pointer. \n
        """

        header = ["\nAtom:"]
        header += ["# of Protons = {}".format(self.z)]
        header += ["# of Nucleons = {}".format(self.a)]
        header += ["Binding Energy of Atom = {}".format(self.bindingEnergy)]
        header += ["Mass of Atom = {}".format(self.mass)]
        header += ["BE per Nucleon = {}".format(self.bindingEnergyPerA)]
        header += ["Mass of Atom per Nucleon = {}".format(self.massPerA)]
        header = "\n".join(header)+"\n"
        return header

    def calcBE(self):
        """!
        @ingroup NuclearStructure
        Calculates the binding energy given a A and Z of an atom according to
        Krane Eq 3.28 given by. \n\n

        Eqn:        \f$ B(Z, A)=a_vA - a_sA^{\frac{2}{3}} - 
                                a_cZ(Z-1)A^{\frac{-1}{3}} - 
                                a_{sym} \frac{(A-2Z)^2}{A} + \delta \f$ \n                   

        @param self: <em> Atom pointer </em> \n
            The Atom pointer. \n
        """

        # Constants
        a_v = 15.5 # MeV
        a_s = 16.8 # MeV
        a_c = 0.72 # MeV
        a_sym = 23.0 # MeV

        if self.a == 1:
            self.bindingEnergy =  0.0
            self.bindingEnergyPerA = 0.0 
        else:
            self.bindingEnergy =  a_v*self.a - a_s*self.a**(2/3.) - \
                a_c*self.z*(self.z-1)*self.a**(-1/3.) - \
                a_sym*(self.a-2*self.z)**2/self.a + self.calcPairing()
            self.bindingEnergyPerA =  self.bindingEnergy/self.a

    def calcPairing(self):
        """!
        @ingroup NuclearStructure
        Calculates the pairing term of the binding energy given a A and Z of
        an atom according to Krane Eq 3.28.  The pairing term is given by. \n\n

        Eqn:        \f$ \delta=a_pA^{\frac{-3}{4}} \f$ for Z and N even, \n
                    \f$ \delta=-a_pA^{\frac{-3}{4}} \f$ for Z and N odd, \n
                    and \f$ \delta=0 \f$ for A odd. \n                   

        @param self: <em> Atom pointer </em> \n
            The Atom pointer. \n

        @return \e float: The binding energy of the atom. \n
        """

        # Constants
        a_p = 34 # MeV

        # Calculate the pairing correction:
        if self.a%2 == 1:
            return 0.0
        elif self.z%2 == 1:
            return -a_p*self.a**(-3/4.)
        elif self.z%2 == 0:
            return a_p*self.a**(-3/4.)
        else:
            print "ERROR: Somehow you have specified a noncovered case in \
                   Atom.calcPairing()."

    def calcSEMF(self):
        """!
        @ingroup NuclearStructure
        Calculates the mass of an atom according to the semi-empirical mass
        formula given by Krane eq. 3.29: \n\n

        Eqn:        \f$ M(Z, A) = Zm(^1H)+Nm_n-\frac{B(Z, A)}{c^2} \f$ \n                   

        @param self: <em> Atom pointer </em> \n
            The Atom pointer. \n
        """

        # Constants
        m_H1 = 1.007825 # u
        m_n = 1.00866501 # u
        c2 = 931.502 # MeV/u

        # Initialize BE
        if self.bindingEnergy == None:
            self.calcBE()
        
        # Calculate SEMF
        self.mass = self.z*m_H1 + (self.a-self.z)*m_n - self.bindingEnergy/c2
        self.massPerA = self.mass/self.a

#------------------------------------------------------------------------------#
def calcMinZ(A):
    """!
    @ingroup NuclearStructure
    Calculated the minimum Z for a given mass chain as given by Krane
    eq. 3.30: \n\n

    Eqn:        \f$ Z_min = \frac{[m_n - m(^1H)] + a_cA^{-frac{1}{3}} +
                4*a_{sym}}{2a_cA^{-frac{1}{3}} + 8A_{sym}A^{-1}}   \f$ \n                   

    @param A: \e integer \n
        The number of nucleons. \n

    @return \e integer: The minimum Z of the mass chain. \n
    """

    # Constants
    m_H1 = 938.7910032 # MeV
    m_n = 939.573 # MeV
    a_c = 0.72 # MeV
    a_sym = 23.0 # MeV
   
    # Calculate Min_Z
    minZ = round(((m_n-m_H1) + a_c*A**(-1/3.) + 4*a_sym)/ \
             (2*a_c*A**(-1/3.) + 8*a_sym*A**(-1)))
    if minZ == 0:
        return 1
    else:
        return int(minZ)