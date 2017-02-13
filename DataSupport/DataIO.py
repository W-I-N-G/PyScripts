#######################################################################################################
#
# Module : DataIO.py
#
# Contains : Routines to read in and print out data 
#
# Author : James Bevins
#
# Last Modified: 24Jan17
#
#######################################################################################################

import re

#-------------------------------------------------------------------------------------------------------------#
def readDelimitedDataFile(path,delimiter=" +",header=0,breakText=""):
    """
    Reads in a set of columated data and returns the results.  
   
    Parameters
    ==========
    path : str
        Absolute path to the file

    Optional
    ========
    delimiter : string
        Indicator for the character used to separate columns.  Uses the re.split delimiter definitions.
        [Deafult = " +"]
    header : integer
        The number of header lines to skip
        [Deafult = 0]
    breakText : string
        Text indicating the end of the data to be read
        [Deafult = ""]
        
    Returns
    =======
    data : list of lists
        A nxm list where n is the number of colums of data in the input file and m is the number of rows
    """ 
    
    assert header>=0, "Valid specifications for the number of header lines to skip must be positive."
    assert type(delimiter)==str, "Valid specifications for the delimiters must be strings."
    
    data=[]
    
    # Open file
    try: 
        f = open(path, 'r') 
        # Skip n header lines
        for i in range(0,header):
            line=f.next()
        
        # Read the file line by line and store the values 
        for line in f:
            if line.rstrip()==breakText:
                break
            split_list=re.split(delimiter,line.strip())   
            for i in range(0,len(split_list)):
                if len(data)<i+1:
                    data.append([])
                data[i].append(float(split_list[i]))

        # Close the file
        f.close()
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror) 
    
    return data