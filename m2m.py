# m2m.py - Main script which calls and runs the other various functions in the proper order to convert molcas orbital files to MOLDEN format.
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

#========
# IMPORTS
#========
# General imports:
#04/09/2023 removed numba import (not needed).
import os
import sys
from itertools import islice
import numpy as np
import h5py
from scipy import linalg as la
# m2m specific imports
import parser
import classes

#=================
# DEFINE FUNCTIONS
#=================
def argsort(lst, rank=None):
    """
    NOTE: This function was taken from: https://github.com/steabert/molpy/blob/master/molpy/basis.py (2023)
    Molpy is an open source project under the terms of the GNU General Public License v2.0.
    Listed Molpy contributors are: Steven Vancoillie, Felix Plasser and Marcus Johansson.

    This function sorts indices of a list """
    return np.array(sorted(np.arange(len(lst)), key=lst.__getitem__))

def sanity_check_sym_mode():
    """ Detects if there is a single or multiple irreps, corrisponding to nosym and sym respectively """
    nbas_vec = InpOrbClass.nbas
    dimension = len(nbas_vec)
    original_setting = input_dictionary['SYM']
    if dimension > 1:
        input_dictionary['SYM'] = True
    if dimension == 1:
        input_dictionary['SYM'] = False
    if not original_setting == input_dictionary['SYM']:
        print("NOTICE: No Symmetry Detected in INPORB file: Will Run in NOSYM Mode.")

#============
# USER INPUTS
#============
# Obtain and manage user inputs:
user_inputs = parser.obtain_user_inputs()
input_dictionary = user_inputs.input_dict()

# Block of code detects if a specific file is given, a single file is in the DIR, or multiple (RECURSIVE_MODE):
if input_dictionary['INPORB_NAME'] == None: #when no files are specified in command line
    user_inputs = parser.no_args()
    input_dictionary = user_inputs.input_dict()
    print("</> RUNNING IN RECURSIVE MODE: " + str(input_dictionary['RECURSIVE_MODE']))
    if input_dictionary['RECURSIVE_MODE'] == False:
        print(input_dictionary['INPORB_SET'])
        pass
    else:
        print(input_dictionary['INPORB_SET'])
        command1 = 'mkdir MOLDENS'
        os.system(command1)
        if input_dictionary['KEEP_HDF5'] == True:
            command2 = 'mkdir HD5'
            os.system(command2)

#--------------------------------------------------
#==================================================
# CONVERT MOLCAS TO MOLDEN OVER ALL FILES PROVIDED:
#==================================================
#--------------------------------------------------
run = 0 #This value should not be changed.

for input_orbital_file in input_dictionary['INPORB_SET']:
    run = classes.calc_log().progress(run, total=len(input_dictionary['INPORB_SET'])) # Tracks number of files converted
    classes.calc_log().message(screen=True, msg='</> Converting ' + str(input_orbital_file) + ' to MOLDEN Format.')
    if input_dictionary['RECURSIVE_MODE'] == True:
        state_number = input_orbital_file.split('.')[-1]
        input_dictionary['NAME'] = 'STATE_' + str(state_number)

    #===================
    #Instancing classes:
    #===================
    classes.calc_log().calc_time(calcstate=0)
    h5_class = classes.h5_data()
    InpOrbClass = classes.InpOrb()
    MoldenClass = classes.molden()
    ToolBoxClass = classes.ToolBox()
    #Defining some class variables:
    InpOrbClass.input_file = str(input_orbital_file)
    InpOrbClass.nMOs = str(input_dictionary['NMOS'])
    InpOrbClass.h5_file = str(input_dictionary['H5'])
    InpOrbClass.nbas = input_dictionary['NBAS_ARRAY']

    #=============================
    # NUMBER OF IRREPS INFORMATION:
    #=============================
    if input_dictionary['FIND_NBAS'] == True:
        #Code finds the number of orbitals per irreps:
        nbas_vec = InpOrbClass.obtain_nbas_data(InpOrbClass.input_file)
        InpOrbClass.nbas = nbas_vec
        print(nbas_vec)
        nMOs = sum(nbas_vec)
        InpOrbClass.nMOs = int(nMOs)
        classes.calc_log().message(screen=True, msg="</> Number of MOs = " + str(nMOs))

    #=========================================
    # OBTAIN AND MANAGE BASIS SET INFORMATION:
    #=========================================
    if input_dictionary['PRIM'] == True:
        #Code obtains basis set information:
        float_formatter = "{:.8e}".format
        np.set_printoptions(formatter={'float_kind':float_formatter})
        prims = np.asarray(h5_class.obtain_primitives(InpOrbClass.h5_file))
        prim_ids = h5_class.obtain_primitive_ids(InpOrbClass.h5_file)

        _basis_ = ""
        _line_breaks_groups_ = []
        _line_breaks_centers_ = []
        group = 0
        center = 0
        group_members = 1
        line_count = 0
        for element in prims:
            line_count += 1
        for i in range(line_count):
            grp = prim_ids[i][2]
            cntr = prim_ids[i][0]
            if cntr == center:
                pass
            else:
                center = cntr
                _line_breaks_centers_.append(i)
                group = 0
            if grp == group:
                pass
            else:
                group = grp
                _line_breaks_groups_.append(i)

        _A_ = _line_breaks_groups_
        _A_.append(len(prim_ids))
        orbs_in_groups = [y-x for x, y in zip(_A_[:-1], _A_[1:])]

        orb_set = {
        '0': "s",
        '1': "p",
        '2': "d",
        '3': "f",
        '4': "g",
        '5': "h"
        }

        _basis_ = ""
        line_count = 0
        break_counter_groups = 0
        break_counter_centers = 0
        orb_type_counter = -1
        for element in prims:
            line_count += 1
            orb_type_counter += 1
            angl = prim_ids[orb_type_counter][1]
            if (line_count-1) in _line_breaks_centers_:
                break_counter_centers += 1
                _ID_INSERT_ = "\n" + "  " + str(prim_ids[line_count-1][0])
                _basis_ = _basis_ + _ID_INSERT_ + "\n"
            if (line_count-1) in _line_breaks_groups_:
                break_counter_groups += 1
                orb_type = str(angl)
                _ID_INSERT_ = "  " + str(orb_set[orb_type]) + "   " + str(orbs_in_groups[break_counter_groups-1])
                _basis_ = _basis_ + _ID_INSERT_ + "\n"
                if orb_set[orb_type] == 's':
                    MoldenClass.s_orb = True
                if orb_set[orb_type] == 'p':
                    MoldenClass.p_orb = True
                if orb_set[orb_type] == 'd':
                    MoldenClass.d_orb = True
                if orb_set[orb_type] == 'f':
                    MoldenClass.f_orb = True
                if orb_set[orb_type] == 'g':
                    MoldenClass.g_orb = True
                if orb_set[orb_type] == 'h':
                    MoldenClass.h_orb = True
            center_number = prim_ids[line_count-1][0]
            angl = prim_ids[line_count-1][1]
            grp = prim_ids[line_count-1][2]

            line = ""
            for item in element:
                item = ToolBoxClass.format_number(item)
                line = line + " " + str(item)
            _basis_ = _basis_ + line + "\n"

        MoldenClass.molden_basis_set = _basis_

    #===========================
    # OBTAIN THE MO COEFFICENTS:
    #===========================
    InpOrbClass.obtain_coeff_data_h5py(InpOrbClass.input_file,InpOrbClass.nMOs)

    #=============================================================
    # SANITY CHECK WHETHER TO RUN IN SYMMETRY OR NO-SYMMETRY MODE:
    #=============================================================
    sanity_check_sym_mode()

    #--------------------------------------------------------
    #========================================================
    #   RUNNING MOLCAS TO MOLDEN CONVERSION IN SYMMETRY MODE:
    #========================================================
    #--------------------------------------------------------
    if input_dictionary['SYM'] == True:
        classes.calc_log().message(screen=True, msg='</> Running in SYM Mode.')

        #======================
        #OBTAIN MO INFORMATION:
        #======================
        #Obtain mo-coefficents
        C_mo,h5_class_mo_vectors = ToolBoxClass.form_matricies_sym(InpOrbClass.nMOs,InpOrbClass.nbas)
        h5_class.mo_vectors = h5_class_mo_vectors
        #Obtain the MO occupancies:
        A_vec = InpOrbClass.obtain_occ_data(InpOrbClass.input_file)

        # Detects if MO energies are given, if not, energies are set to zero:
        if input_dictionary['ENERGIES_FROM_HD5'] == True:
            name = InpOrbClass.h5_file
            name = name.split('.h5')[0]
            hdf5_mo_energy_data = h5_class.read_data(name, key='MO_ENERGIES')
            mo_energies = hdf5_mo_energy_data
        else:
            try:
                with open(InpOrbClass.input_file, 'r') as f:
                    data = f.readlines()
                file_found = True
            except:
                print("ERROR: Check for mo energies failed to read INPORB FILE.")
                file_found = False

            energies_present = False
            if file_found == True:
                for line in data:
                    if "#ONE" in line:
                        energies_present = True
                        classes.calc_log().message(screen=True, msg='</> MO Energies are in INPORB')
                    else:
                        pass

            if energies_present == True:
                mo_energies = InpOrbClass.obtain_energy_data(InpOrbClass.input_file)
            else:
                classes.calc_log().message(screen=True, msg='</> Will build energy data setting MO energies to zero.')
                mo_energies = InpOrbClass.forge_energy_data(InpOrbClass.input_file)
                mo_energies = InpOrbClass.obtain_energy_data(InpOrbClass.input_file)

        #Store MO energies and occupations:
        h5_class.mo_energies = mo_energies
        h5_class.mo_occupations = A_vec
        h5_class.calc_type_indicies()

        #obtain the symmetry adapted linear combination matrix and apply to mo-coefficents:
        salc_matrix = h5_class.obtain_salcs(InpOrbClass.h5_file,InpOrbClass.nbas)
        C_SYM = np.matmul(salc_matrix, C_mo)
        mo_vectors = np.asarray(C_SYM)

        #==================================================
        # SORT THE MO-COEFFICENT MATRIX INTO MOLDEN FORMAT:
        #==================================================

        #========================== sorting the coefficent matrix to that of the molden file format ===========
        desym_basis_function_ids = h5_class.obtain_desym_basis_function_ids(InpOrbClass.h5_file)
        '''
        Start of code snippet.

        NOTE: The following code snippet is taken from: https://github.com/steabert/molpy/blob/master/molpy/basis.py (2023)
              specifically taken from molpy basis.py:     def _idtuples_updown_order(self):
              Molpy is an open source project under the terms of the GNU General Public License v2.0.
              Listed Molpy contributors are: Steven Vancoillie, Felix Plasser and Marcus Johansson.
        '''
        idtuples = []
        for id_tuple in desym_basis_function_ids:
            center, n, l, m, = id_tuple
            if l == 1 and m == 0:
                m = 2
            if m < 0:
                m = -m + 0.5
            idtuples.append((center, l, n, m))
        '''
        End of code snippet.
        '''
        nmos = int(InpOrbClass.nMOs)
        ids = np.arange(0,nmos)
        cgto_molden_rank = np.argsort(argsort(idtuples)) #combination of np.argsort and an argsort func define at the top of script
        rank = cgto_molden_rank[ids]
        rank = np.argsort(rank)
        mo_molden_vecs = mo_vectors[rank,:]
        molden_vectors = mo_molden_vecs.T
        molden_occs = InpOrbClass.obtain_occ_data(InpOrbClass.input_file)
        molden_nbas = InpOrbClass.nbas
        molden_energies = h5_class.mo_energies
        molden_label = h5_class.obtain_desym_center_labels(InpOrbClass.h5_file)
        molden_coords = h5_class.obtain_desym_center_coordinates(InpOrbClass.h5_file)
        molden_atnum = h5_class.obtain_desym_center_atnum(InpOrbClass.h5_file)
        molden_orb_set_dict = MoldenClass.orb_set_dict()
        molden_basis_set = MoldenClass.molden_basis_set
        molden_name = input_dictionary['NAME']
        enocc = input_dictionary['SET_ENERGY_TO_OCC']
        sym_labels = input_dictionary['SYM_LABELS']
        MoldenClass.gen_molden(molden_vectors, molden_occs, molden_nbas,
                               molden_energies, molden_label, molden_coords,
                               molden_atnum, molden_orb_set_dict,
                               molden_basis_set, enocc, sym_labels, molden_name)

        #==========================================================
        # GENERATE A DENSITY MATRIX FROM MO-COEFFS AND OCCUPANCIES:
        #==========================================================
        if input_dictionary['CALC_DENSITY'] == True:
            A_vec = np.asarray(A_vec)
            nMOs = int(InpOrbClass.nMOs)

            #Obtain the overlap matrix:
            with h5py.File(str(InpOrbClass.h5_file), 'r') as f:
                dset = f['AO_OVERLAP_MATRIX']
                arr = np.array(dset[:])
                n_bas = InpOrbClass.nbas
                overlap = ToolBoxClass.reshape_square(arr, n_bas)
                f.close()

            dens = ToolBoxClass.density_from_momatrix(mo_vectors, A_vec)
            square = ToolBoxClass.density_as_square(denvec=dens)
            overlap = np.dot(np.dot(salc_matrix, overlap), salc_matrix.T)
            MO_Den, nonzero = ToolBoxClass.MO_Dens_SYM(mo_vectors, square, overlap,A_vec)
            h5_class.format_density_matrix(density=nonzero)
            np.set_printoptions(suppress=False)
            np.set_printoptions(precision=10)

    #-----------------------------------------------------------
    #===========================================================
    #   RUNNING MOLCAS TO MOLDEN CONVERSION IN NO-SYMMETRY MODE:
    #===========================================================
    #-----------------------------------------------------------
    if input_dictionary['SYM'] == False:
        classes.calc_log().message(screen=True, msg='</> Running in NOSYM Mode.')

        #===========================
        #OBTAIN ORBITAL INFORMATION:
        #===========================
        #Obtain orbital-coefficents
        C_mo = ToolBoxClass.form_matricies(InpOrbClass.nMOs)
        h5_class.mo_vectors = C_mo
        #Obtain orbital occupancies
        A_vec = InpOrbClass.obtain_occ_data(InpOrbClass.input_file)
        #re-format
        C_mo = np.asarray(C_mo)
        A_vec = np.asarray(A_vec)
        #save:
        InpOrbClass.C_mo = C_mo
        InpOrbClass.A_vec = A_vec

        # Trys to open the provided INPORB file:
        try:
            with open(InpOrbClass.input_file, 'r') as f:
                data = f.readlines()
            file_found = True
        except:
            print("ERROR: Check for mo energies failed to read INPORB FILE.")
            file_found = False

        # Detects if MO energies are given, if not, energies are set to zero:
        energies_present = False
        if file_found == True:
            for line in data:
                if "#ONE" in line:
                    energies_present = True
                    classes.calc_log().message(screen=True, msg='</> MO Energies are in INPORB')
                else:
                    pass

        if energies_present == True:
            mo_energies = InpOrbClass.obtain_energy_data(InpOrbClass.input_file)
        else:
            classes.calc_log().message(screen=True, msg='</> Will build energy data setting MO energies to zero.')
            mo_energies = InpOrbClass.forge_energy_data(InpOrbClass.input_file)
            mo_energies = InpOrbClass.obtain_energy_data(InpOrbClass.input_file)

        #Store MO energies and occupations:
        h5_class.mo_energies = mo_energies
        h5_class.mo_occupations = A_vec
        h5_class.calc_type_indicies()

        #==========================================================
        # GENERATE A DENSITY MATRIX FROM MO-COEFFS AND OCCUPANCIES:
        #==========================================================
        if input_dictionary['CALC_DENSITY'] == True:
            #Obtain the overlap matrix:
            nMOs = int(InpOrbClass.nMOs)
            with h5py.File(str(InpOrbClass.h5_file), 'r') as f:
                dset = f['AO_OVERLAP_MATRIX']
                ncols = nMOs
                ao_overlap = np.reshape(dset, (-1, ncols))
                f.close()

            dens = ToolBoxClass.density_from_momatrix(C_mo,A_vec)
            square = ToolBoxClass.density_as_square(dens)
            MO_Den = ToolBoxClass.MO_Dens(C_mo, square, ao_overlap,A_vec)
            h5_class.format_density_matrix(density=MO_Den)

        #==================================================
        # SORT THE MO-COEFFICENT MATRIX INTO MOLDEN FORMAT:
        #==================================================
        C_mo = InpOrbClass.C_mo
        A_vec = InpOrbClass.A_vec
        mo_vectors = C_mo
        #========================== sorting the coefficent matrix to that of the molden file format ===========
        desym_basis_function_ids = h5_class.obtain_basis_function_ids(InpOrbClass.h5_file)
        '''
        Start of code snippet.

        NOTE: The following code snippet is taken from: https://github.com/steabert/molpy/blob/master/molpy/basis.py (2023)
              specifically taken from molpy basis.py:     def _idtuples_updown_order(self):
              Molpy is an open source project under the terms of the GNU General Public License v2.0.
              Listed Molpy contributors are: Steven Vancoillie, Felix Plasser and Marcus Johansson.
        '''
        idtuples = []
        for id_tuple in desym_basis_function_ids:
            center, n, l, m, = id_tuple
            if l == 1 and m == 0:
                m = 2
            if m < 0:
                m = -m + 0.5
            idtuples.append((center, l, n, m))
        '''
        End of code snippet
        '''
        nmos = int(InpOrbClass.nMOs)
        ids = np.arange(0,nmos)
        cgto_molden_rank = np.argsort(argsort(idtuples)) #combination of np.argsort and an argsort func define at the top of script
        rank = cgto_molden_rank[ids]
        rank = np.argsort(rank)
        mo_molden_vecs = mo_vectors[rank,:]
        molden_vectors = mo_molden_vecs.T
        molden_occs = InpOrbClass.obtain_occ_data(InpOrbClass.input_file)
        molden_nbas = InpOrbClass.nbas
        molden_energies = h5_class.mo_energies
        molden_label = h5_class.obtain_center_labels(InpOrbClass.h5_file)
        molden_coords = h5_class.obtain_center_coordinates(InpOrbClass.h5_file)
        molden_atnum = h5_class.obtain_center_atnum(InpOrbClass.h5_file)
        molden_orb_set_dict = MoldenClass.orb_set_dict()
        molden_basis_set = MoldenClass.molden_basis_set
        molden_name = input_dictionary['NAME']
        enocc = input_dictionary['SET_ENERGY_TO_OCC']
        sym_labels = input_dictionary['SYM_LABELS']
        MoldenClass.gen_molden(molden_vectors, molden_occs, molden_nbas,
                               molden_energies, molden_label, molden_coords,
                               molden_atnum, molden_orb_set_dict,
                               molden_basis_set, enocc, sym_labels, molden_name)

    #------------------------------------------
    #==========================================
    # END OF CONVERSION, DATA HANDLING SECTION:
    #==========================================
    #------------------------------------------

    #============================
    # GENERATES m2m HDF5 FILES:
    #============================
    if input_dictionary['H5_UPDATE'] == True:
        name = input_dictionary['NAME'] # will create a final .h5 file with the name given in settings.py
        h5_class.update(InpOrbClass.h5_file, name, input_dictionary)

    #code runs a series of system commands to manage molden and HDF5 files:
    if input_dictionary['RECURSIVE_MODE'] == True:
        name = str(input_dictionary['NAME'])
        command = 'mv m2m_DATA.h5 ' + str(name) + '_DATA.h5'
        os.system(command)

        name = str(input_dictionary['NAME'])
        command1 = 'mv ' + str (name) + '.molden MOLDENS'
        os.system(command1)

        if input_dictionary['KEEP_HDF5'] == False:
            command2 = 'rm ' + str (name) + '.h5'
        if input_dictionary['KEEP_HDF5'] == True:
            command2 = 'mv ' + str (name) + '.h5 HD5'
        os.system(command2)

        if input_dictionary['KEEP_HDF5'] == False:
            command3 = 'rm ' + str(name) + '_DATA.h5'
        if input_dictionary['KEEP_HDF5'] == True:
            command3 = 'mv ' + str (name) + '_DATA.h5 HD5'
        os.system(command3)

    classes.calc_log().message(screen=True, msg="</> Finished!.")
    classes.calc_log().calc_time(calcstate=1)

exit()
