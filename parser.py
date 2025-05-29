# parser.py - will take user inputs/options from the command line as opposed to using the standard settings.
#
#
#     Copyright (c) 2023  Kurtis Stanistreet-Welsh
#
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#     Code written by Kurtis Stanistreet-Welsh,
#     all corrispondence should be directed to
#     k.stanistreet.welsh@gmail.com
#
#     Thanks goes to Dr Andrew Kerridge for guidance and helpful discussions.

#=========
# IMPORTS:
#=========
#General imports:
import os
import sys
import numpy as np
#m2m specific imports:
import classes

#---------
#=========
# CLASSES:
#=========
#---------
class track_inputs:
    '''
    Class: track_inputs:
    Description: Keeps track of the code settings

    track_inputs:
                |input_dict()

    - Function: input_dict()
      Description: This function takes the settings stored in __init__() and creates a dictionary which is returned for use in m2m.py
    '''
    def __init__(
    self, HELP=None, INPORB_NAME=None,
    CLEAN=None, NMOS=None,
    H5=None, H5_UPDATE=None,
    SYM=None, FIND_NBAS=True,
    NBAS_ARRAY=None, PRIM=False, NAME=None,
    RECURSIVE_MODE=False, INPORB_SET=None,
    CALC_DENSITY=False, SET_ENERGY_TO_OCC=False,
    ENERGIES_FROM_HD5=False, SYM_LABELS=True,
    KEEP_HDF5=False
    ):
        self.HELP = False
        self.INPORB_NAME = None
        self.CLEAN = False
        self.NMOS = None
        self.H5 = None
        self.H5_UPDATE = False
        self.SYM = False
        self.FIND_NBAS = True
        self.NBAS_ARRAY = None
        self.PRIM = False
        self.NAME='m2m'
        self.RECURSIVE_MODE = False
        self.INPORB_SET = None # will contain a list of INPORB files
        self.CALC_DENSITY = False
        self.SET_ENERGY_TO_OCC = False
        self.ENERGIES_FROM_HD5 = False
        self.SYM_LABELS = True
        self.KEEP_HDF5 = False
    def input_dict(self):
        input_dict = {
        'HELP': self.HELP,
        'INPORB_NAME': self.INPORB_NAME,
        'CLEAN': self.CLEAN,
        'NMOS': self.NMOS,
        'H5': self.H5,
        'H5_UPDATE': self.H5_UPDATE,
        'FIND_NBAS': self.FIND_NBAS,
        'SYM': self.SYM,
        'NBAS_ARRAY': self.NBAS_ARRAY,
        'PRIM': self.PRIM,
        'NAME': self.NAME,
        'RECURSIVE_MODE': self.RECURSIVE_MODE,
        'INPORB_SET': self.INPORB_SET,
        'CALC_DENSITY': self.CALC_DENSITY,
        'SET_ENERGY_TO_OCC': self.SET_ENERGY_TO_OCC,
        'ENERGIES_FROM_HD5': self.ENERGIES_FROM_HD5,
        'SYM_LABELS': self.SYM_LABELS,
        'KEEP_HDF5' : self.KEEP_HDF5,
        }
        return input_dict

#-----------
#===========
# FUNCTIONS:
#===========
#-----------
def no_args():
    '''
    Function: no_args()
    Description: Generates code settings from the settings.py file if no command line options given.
    '''
    #initialise class
    store_inputs = track_inputs()
    #stores command line inputs
    arg_array = sys.argv[:]
    use_settings = True # Hard-set so the code makes use of the settings.py file
    import settings
    classes.calc_log().settings()
    #sets the various code settings using the settings.py file:
    if use_settings == True:
        store_inputs.HELP = settings.HELP
        store_inputs.INPORB_NAME = settings.INPORB_NAME
        store_inputs.CLEAN = settings.CLEAN
        store_inputs.NMOS = settings.NMOS
        store_inputs.H5 = settings.H5
        store_inputs.H5_UPDATE = settings.H5_UPDATE
        store_inputs.SYM = settings.SYM
        store_inputs.FIND_NBAS = settings.FIND_NBAS
        store_inputs.NBAS_ARRAY = settings.NBAS_ARRAY
        store_inputs.PRIM = settings.PRIM
        store_inputs.NAME = settings.NAME
        store_inputs.RECURSIVE_MODE = settings.RECURSIVE_MODE
        store_inputs.INPORB_SET = settings.INPORB_SET
        store_inputs.CALC_DENSITY = settings.CALC_DENSITY
        store_inputs.SET_ENERGY_TO_OCC = settings.SET_ENERGY_TO_OCC
        store_inputs.ENERGIES_FROM_HD5 = settings.ENERGIES_FROM_HD5
        store_inputs.SYM_LABELS = settings.SYM_LABELS
        store_inputs.KEEP_HDF5 = settings.KEEP_HDF5
        classes.calc_log().message(screen=True, msg='</> No Command Line Options Given: Settings Loaded')
    else:
        classes.calc_log().message(screen=True, msg='</> Settings Not Used.')
        pass

    #Next block of code examines the working dir for INPORB and HDF5 files:
    cwd = os.getcwd()
    files = os.listdir(cwd)
    INPORB_COUNT = 0
    H5_COUNT = 0
    INPORB_SET = []
    for file in files:
        if '.RasOrb' in file:
            INPORB_COUNT += 1
            classes.calc_log().message(screen=True, msg='</> Detected RasOrb File Type: ' + str(file))
            store_inputs.INPORB_NAME = str(file)
            INPORB_SET.append(str(file))
        #04/09/2023 added support for Localised Orbital files
        if '.LocOrb' in file:
            INPORB_COUNT += 1
            classes.calc_log().message(screen=True, msg='</> Detected LocOrb File Type: ' + str(file))
            store_inputs.INPORB_NAME = str(file)
            INPORB_SET.append(str(file))
        if '.ScfOrb' in file:
            INPORB_COUNT += 1
            classes.calc_log().message(screen=True, msg='</> Detected ScfOrb File Type: ' + str(file))
            store_inputs.INPORB_NAME = str(file)
            INPORB_SET.append(str(file))
        if '.SOOrb' in file:
            INPORB_COUNT += 1
            classes.calc_log().message(screen=True, msg='</> Detected SOOrb File Type: ' + str(file))
            store_inputs.INPORB_NAME = str(file)
            INPORB_SET.append(file)
        if '.h5' in file:
            H5_COUNT += 1
            classes.calc_log().message(screen=True, msg='</> Detected HDF5 File: ' + str(file))
            store_inputs.H5 = str(file)

    #Determines whether there is a single INPORB or multiple INPORB files (RECURSIVE_MODE):
    if INPORB_COUNT > 1:
        classes.calc_log().message(screen=True, msg='</> NOTICE: More than One Acceptable INPORB File in Current DIR.')
        store_inputs.RECURSIVE_MODE = True
        store_inputs.INPORB_SET = INPORB_SET
    else:
        store_inputs.RECURSIVE_MODE = False
        store_inputs.INPORB_SET = INPORB_SET
    if H5_COUNT > 1:
        classes.calc_log().message(screen=True, msg='</> ERROR: More than one acceptable HD5 (.h5) file in current DIR.')
        exit()
    #returns the code settings:
    return store_inputs

def obtain_user_inputs():
    '''
    Function: obtain_user_inputs()
    Description: Determines the settings used for conversion via:
                              - User inputs from the command line
                              - From the settings.py file
    '''
    store_inputs = track_inputs()
    arg_array = sys.argv[:]
    use_settings = True # Hard-set so that the settings file is used initially.
    import settings
    classes.calc_log().settings()
    #Sets the code settings from the settings.py file:
    if use_settings == True:
        store_inputs.HELP = settings.HELP
        store_inputs.INPORB_NAME = settings.INPORB_NAME
        store_inputs.CLEAN = settings.CLEAN
        store_inputs.NMOS = settings.NMOS
        store_inputs.H5 = settings.H5
        store_inputs.H5_UPDATE = settings.H5_UPDATE
        store_inputs.SYM = settings.SYM
        store_inputs.FIND_NBAS = settings.FIND_NBAS
        store_inputs.NBAS_ARRAY = settings.NBAS_ARRAY
        store_inputs.PRIM = settings.PRIM
        store_inputs.NAME = settings.NAME
        store_inputs.RECURSIVE_MODE = settings.RECURSIVE_MODE
        store_inputs.INPORB_SET = settings.INPORB_SET
        store_inputs.CALC_DENSITY = settings.CALC_DENSITY
        store_inputs.SET_ENERGY_TO_OCC = settings.SET_ENERGY_TO_OCC
        store_inputs.ENERGIES_FROM_HD5 = settings.ENERGIES_FROM_HD5
        store_inputs.SYM_LABELS = settings.SYM_LABELS
        store_inputs.KEEP_HDF5 = settings.KEEP_HDF5
        classes.calc_log().message(screen=True, msg='</> Settings Loaded')
    else:
        classes.calc_log().message(screen=True, msg='</> Settings Not Used.')
        pass

    #=============
    # USER INPUTS:
    #=============
    # Examines the user inputs to detect a specific INPORB file and HDF5 file:
    i = -1
    for argument in arg_array:
        i += 1
        INPORB_SET = []
        print('USER INPUT: ' + str(argument))
        if '.RasOrb' in argument:
            file_name = arg_array[i]
            classes.calc_log().message(screen=True, msg='</> Detected RasOrb File Type: ' + str(file_name))
            store_inputs.INPORB_NAME = str(file_name)
            INPORB_SET.append(str(file_name))
            store_inputs.INPORB_SET = INPORB_SET
        #04/09/2023 added support for Localised Orbital files
        if '.LocOrb' in argument:
            file_name = arg_array[i]
            classes.calc_log().message(screen=True, msg='</> Detected LocOrb File Type: ' + str(file_name))
            store_inputs.INPORB_NAME = str(file_name)
            INPORB_SET.append(str(file_name))
            store_inputs.INPORB_SET = INPORB_SET
        if '.ScfOrb' in argument:
            file_name = arg_array[i]
            classes.calc_log().message(screen=True, msg='</> Detected ScfOrb File Type: ' + str(file_name))
            store_inputs.INPORB_NAME = str(file_name)
            INPORB_SET.append(str(file_name))
            store_inputs.INPORB_SET = INPORB_SET
        if '.SOOrb' in argument:
            file_name = arg_array[i]
            classes.calc_log().message(screen=True, msg='</> Detected SOOrb File Type: ' + str(file_name))
            store_inputs.INPORB_NAME = str(file_name)
            INPORB_SET.append(str(file_name))
            store_inputs.INPORB_SET = INPORB_SET
        if '.GssOrb' in argument:
            file_name = arg_array[i]
            classes.calc_log().message(screen=True, msg='</> Detected GssOrb File Type: ' + str(file_name))
            store_inputs.INPORB_NAME = str(file_name)
            INPORB_SET.append(str(file_name))
            store_inputs.INPORB_SET = INPORB_SET
        if '.h5' in argument:
            file_name = arg_array[i]
            classes.calc_log().message(screen=True, msg='</> Detected HDF5 File: ' + str(file_name))
            store_inputs.H5 = str(file_name)

        #Examines the user command line inputs to change specific settings:
        if 'nmos' in argument:
            nMOs = arg_array[i]
            try:
                nMOs = nMOs.split('=')
                nMOs = int(nMOs[1])
                store_inputs.NMOS = int(nMOs)
            except:
                print(
                "ERROR: Inproper formatting" + "\n"
                "     Please follow this example:" + "\n"
                "          nmos=5" + "\n"
                "No spaces in the input"
                )

        if 'nbas' in argument:
            nbas = arg_array[i]
            store_inputs.FIND_NBAS = False
            try:
                nbas = nbas.split('=')
                nbas = str(nbas[1])
                nbas = nbas.split(',')
                n_bas = []
                for element in nbas:
                    n_bas.append(element)
                n_bas = np.asarray(n_bas)
                n_bas = n_bas.astype(np.int64)
                print(n_bas)
                store_inputs.NBAS_ARRAY = n_bas
            except:
                print(
                "ERROR: Inproper formatting" + "\n"
                "     Please follow this example:" + "\n"
                "          nbas=7,3,3,1,7,3,3,1" + "\n"
                "No spaces in the input"
                )
                exit()

            if 'nmos' in argument:
                pass
            else:
                nMOs = sum(store_inputs.NBAS_ARRAY)
                store_inputs.NMOS = int(nMOs)
                classes.calc_log().message(screen=True, msg="</> Number of MOs = " + str(nMOs))

        if '-help' in argument:
            store_inputs.HELP = True
            print('Coming soon...'
            )
            exit()

        if '-clean' in argument:
            store_inputs.CLEAN = True

            os.system('rm lines')
            os.system('rm MO_*')
            os.system('rm *.tmp')
            exit()

        if '-h5' in argument:
            store_inputs.H5_UPDATE = True

        if '-prim' in argument:
            store_inputs.PRIM = True

        if '-sym' in argument:
            store_inputs.SYM = True

        if '-nosym' in argument:
            store_inputs.SYM = False

        if '-noden' in argument:
            store_inputs.CALC_DENSITY = False

        if '-enocc' in argument:
            store_inputs.SET_ENERGY_TO_OCC = True

        if '-enh5' in argument:
            store_inputs.ENERGIES_FROM_HD5 = True

        if '-sym_labels' in argument:
            store_inputs.SYM_LABELS = True

    #Return the user specified settings for use in m2m.py
    return store_inputs
