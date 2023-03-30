# classes.py - contains the main classes used overall by m2m.py
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
#General imports
import os
import h5py
import datetime
#Required for Class functions:
from itertools import islice
import numpy as np
from scipy import linalg as la
#m2m specific imports:
import settings #This loads in the settings file

#---------
#=========
# CLASSES:
#=========
#---------

class calc_log:
    '''
    Class: calc_log:
    Description: This class manages information either displayed to screen or written to the log file.

    calc_log:
            |message()
            |write_to_log()
            |calc_time()
            |settings()
            |progress()

    - Function: message()
      Description: Prints messages to screen.

    - Function: write_to_log()
      Description: Writes messages to the log file.

    - Function: calc_time()
      Description: Writes the start and finish times to the log file.

    - Function: settings()
      Description: Writes the contents of the settings.py file to the log file.

    - Progress()
      Description: Keeps track of overall progress

    '''
    def __init__(self, NAME=None, TIME=None):
        self.NAME='m2m.log'
        self.TIME=None

    def message(self, screen=True, msg=None):
        '''
        Function: message()
        Description: Prints messages to screen.
        '''
        if screen == True:
            print(msg)
            self.write_to_log(msg + '\n')
        else:
            self.write_to_log(msg + '\n')

    def write_to_log(self, log):
        '''
        Function: write_to_log()
        Description: Writes messages to the log file.
        '''
        with open(str(self.NAME), 'a') as f:
            f.write(str(log))

    def calc_time(self, calcstate=3):
        '''
        Function: calc_time()
        Description: Writes the start and finish times to the log file.
        '''
        self.TIME=str(datetime.datetime.now())
        if calcstate == 0:
            log = (
            '\n' +
            'Molcas2Molden Started: ' + self.TIME + '\n'
            '\n'
            )
        elif calcstate == 1:
            log = (
            '\n' +
            'Molcas2Molden Finished: ' + self.TIME + '\n'
            '\n'
            )
        else:
            log = str(
            '\n' +
            'Molcas2Molden FAILED: ' + self.TIME + '\n' +
            '\n'
            )
        self.write_to_log(log)

    def settings(self):
        '''
        Function: settings()
        Description: Writes the contents of the settings.py file to the log file.
        '''
        log = (
        'Settings are taken from settings.py :' + '\n'
        '::SETTINGS::' + '\n'
        )
        for attr in vars(settings):
            var = attr
            if '_' in var:
                pass
            else:
                val = getattr(settings,attr)
                log = log + str(attr) + ' = ' + str(val) + '\n'
        self.write_to_log(log)

    def progress(self, run, total):
        '''
        Function: progress()
        Description: Keeps track of overall progress
        '''
        run = run + 1
        log = '===================' + '\n'
        log = log + 'Converting: [' + str(run) + ' | ' + str(total) + ']' + '\n'
        log = log + '==================='
        self.write_to_log(log)
        self.message(screen=True, msg=log)
        return run

class InpOrb:
    '''
    Class: InpOrb:
    Description: Manages the reading of data from the INPORB file supplied by the user.

    InpOrb:
          |obtain_coeff_data_h5py()
          |forge_energy_data()
          |obtain_occ_data()
          |obtain_energy_data()
          |unique()
          |obtain_nbas_data()

    - Function: obtain_coeff_data_h5py()
      Description: This function opens the user-supplied INPORB file, reads its contents, extracts the MO coefficient values,
                   and then stores these values in a new HDF5 file (m2m.h5) for further use.

    - Function: forge_energy_data()
      Description: When no energy data is present in the user-supplied INPORB file
                   (this is the case for SOOrb files), this function generates energy data by setting
                   all MOs to have energy values of zero and then appends the generated data to the INPORB file.

    - Function: obtain_occ_data()
      Description: This function opens the user-supplied INPORB file, reads the MO occupation values from it,
                   stores them in a list, and returns the list as the output.

    - Function: obtain_energy_data()
      Description: This function opens the user-supplied INPORB file, reads the MO energy values from it,
                   stores them in a list, and returns the list as the output.

    - Function: unique()
      Description: Takes the supplied list and returns a new list containing only unique elements from
                   the original list. Effectively removes duplicates from the original list and returns this.

    - Function: obtain_nbas_data()
      Description: Opens the supplied INPORB file, reads the number of basis functions (NBAS) line,
                   stores it as a list and returns it.
    '''
    def __init__(
    self, input_file=None, nMOs=None,
    nMO_coeffs=None, C_mo=None, A_vec=None,
    h5_file=None, NatOccs = None, NatOrbs = None,
    nbas=None, primitives = None
    ):
        self.input_file = None
        self.nMOs = None
        self.nMO_coeffs = None
        self.C_mo = None
        self.A_vec = None
        self.h5_file = None
        self.NatOccs = None
        self.NatOrbs = None
        self.nbas = None
        self.primitives = None

    def obtain_coeff_data_h5py(self,InpOrbClass_input_file,InpOrbClass_nMOs):
        """
        Function: obtain_coeff_data_h5py()
        Description: This function opens the user-supplied INPORB file, reads its contents, extracts the MO coefficient values,
        and then stores these values in a new HDF5 file (m2m.h5) for further use.
        """
        with open(InpOrbClass_input_file, 'r') as f:
            desym_data = f.readlines()
        orbital_lines = []
        nline = 0
        for line in desym_data:
            nline += 1
            if "* ORBITAL" in line:
                orbital_lines.append(nline)
            if "#OCC" in line:
                orbital_lines.append(nline)
        pairs = []
        start_orb = orbital_lines[0]
        end_orb = orbital_lines[1] - 1
        pair = str(start_orb) + ";" + str(end_orb)
        MOs = int(InpOrbClass_nMOs)
        for i in range(int(MOs)):
            start_orb = orbital_lines[i]
            end_orb = orbital_lines[i+1] - 1
            pair = str(start_orb) + ";" + str(end_orb)
            pairs.append(pair)

        h5_data().save_data(name='m2m_DATA', key='m2m_PAIRS', contents=pairs)

        calc_log().message(screen=True, msg="</> Extracting Orbital Data From INPORB File.")
        for i in range(len(pairs)):
            percentage = (int(i)/len(pairs)) * 100
            print(str(round(percentage)) + '%', end='\r')
            line_pair = pairs[i]
            line_pair = line_pair.split(';')
            start_orb = line_pair[0]
            end_orb = line_pair[1]
            file_name_save_mo = "MO_" + str(i+1)
            h5_line_data = []
            with open(InpOrbClass_input_file) as fin:
                lines = islice(fin, int(start_orb), int(end_orb)) # or whatever ranges
                for line in lines:
                    h5_line_data.append(line)
            h5_data().save_data(name='m2m_DATA', key=str(file_name_save_mo), contents=h5_line_data)
        calc_log().message(screen=False, msg=(str(round(percentage)) + '%'))
        calc_log().message(screen=True, msg="</> Storing Orbital Data to HDF5 (m2m_DATA.h5).")

        for i in range(len(pairs)):
            percentage = (int(i)/len(pairs)) * 100
            print(str(round(percentage)) + '%', end='\r')
            file_name = "MO_" + str(i+1)

            key_name = file_name

            h5_list_data = []

            dset_byte_array = h5_data().read_data(name='m2m_DATA', key=file_name)
            dset_string_list = [s.decode().strip() for s in dset_byte_array]

            mo_coeffs = dset_string_list
            for line in mo_coeffs:
                line = line.split()
                for data in line:
                    h5_list_data.append(data)
                    file_name_coeff = "MO_LIST_" + str(i+1)
            h5_data().save_data(name='m2m_DATA', key=str(file_name_coeff), contents=h5_list_data)
        calc_log().message(screen=False, msg=(str(round(percentage)) + '%'))

    def forge_energy_data(self,InpOrbClass_input_file):
        '''
        Function: forge_energy_data()
        Description: When no energy data is present in the user-supplied INPORB file
        (this is the case for SOOrb files), this function generates energy data by setting
        all MOs to have energy values of zero and then appends the generated data to the INPORB file.
        '''
        calc_log().message(screen=True, msg="</> Appending the INPORB FILE to add zero energies for MO energies.")
        with open(InpOrbClass_input_file, 'r') as f:
            data = f.readlines()

        end_line = len(data)
        start_line = None
        line_counter = 0
        for line in reversed(data):
            line_counter += 1
            if 'HUMAN' in line:
                start_line = end_line - line_counter + 1

        data = []
        with open(InpOrbClass_input_file) as fin:
            lines = islice(fin, int(start_line), int(end_line)) # or whatever ranges
            for line in lines:
                data.append(line)

        energy_build = '#ONE' + "\n" + '* ONE ELECTRON ENERGIES"'
        for line in data:
            energy_build = energy_build + '\n'
            this_line = line.split()
            for element in this_line:
                energy_build = energy_build + " " + "0.0000" + ""

        energy_build = energy_build + "\n" + '#INDEX'

        with open(InpOrbClass_input_file, 'a') as f:
            f.write(energy_build)

    def obtain_occ_data(self,InpOrbClass_input_file):
        '''
        Function: obtain_occ_data()
        Description: This function opens the user-supplied INPORB file, reads the MO occupation values from it,
        stores them in a list, and returns the list as the output.
        '''
        calc_log().message(screen=True, msg="</> Obtaining MO Occupation data.")
        with open(InpOrbClass_input_file, 'r') as f:
            desym_data = f.readlines()
        occ_lines = []
        nline = 0
        for line in desym_data:
            nline += 1
            if "#OCC" in line:
                occ_lines.append(nline)
            if "#OCHR" in line:
                occ_lines.append(nline)
        pairs = []
        start_occ = occ_lines[0] + 2
        end_occ = occ_lines[1] - 1

        A_vec = []
        with open(InpOrbClass_input_file) as fin:
            lines = islice(fin, int(start_occ-1), int(end_occ)) # or whatever ranges
            for line in lines:
                mo_occs = line.split()
                for data in mo_occs:
                    A_vec.append(float(data))
        return A_vec

    def obtain_energy_data(self,InpOrbClass_input_file):
        '''
        Function: obtain_energy_data()
        Description: This function opens the user-supplied INPORB file, reads the MO energy values from it,
        stores them in a list, and returns the list as the output.
        '''
        calc_log().message(screen=True, msg="</> Obtaining MO Energy data.")
        with open(InpOrbClass_input_file, 'r') as f:
            desym_data = f.readlines()
        energy_lines = []
        nline = 0
        for line in desym_data:
            nline += 1
            if "#ONE" in line:
                energy_lines.append(nline)
            if "#INDEX" in line:
                energy_lines.append(nline)
        pairs = []
        start_energy = energy_lines[0] + 2
        end_energy = energy_lines[1] - 1

        energy_vec = []
        with open(InpOrbClass_input_file) as fin:
            lines = islice(fin, int(start_energy-1), int(end_energy)) # or whatever ranges
            for line in lines:
                mo_energies = line.split()
                for data in mo_energies:
                    energy_vec.append(float(data))
        return energy_vec

    def unique(self,list1):
        '''
        Function: unique()
        Description: Takes the supplied list and returns a new list containing only unique elements from
        the original list. Effectively removes duplicates from the original list and returns this.
        '''
        unique_list = []
        for x in list1:
            if x not in unique_list:
                unique_list.append(x)
        return unique_list

    def obtain_nbas_data(self, InpOrbClass_input_file):
        '''
        Function: obtain_nbas_data()
        Description: Opens the supplied INPORB file, reads the number of basis functions (NBAS) line,
        stores it as a list and returns it.
        '''
        calc_log().message(screen=True, msg="</> Obtaining NBAS Data:")
        with open(InpOrbClass_input_file, 'r') as f:
            inputorb_data = f.readlines()
        irreps = []
        nbas = []
        for line in inputorb_data:
            if "ORBITAL" in line:
                line = line.split()
                sym = int(line[2])
                irreps.append(sym)
        sym = self.unique(irreps)
        length = len(sym)
        nbas_floats = np.zeros(length)
        nbas = []
        for line in inputorb_data:
            if "ORBITAL" in line:
                line = line.split()
                sym = int(line[2])
                element = nbas_floats[sym-1] + 1
                nbas_floats[sym-1] = element
        for float in nbas_floats:
            float = int(float)
            nbas.append(float)
        return nbas

class molden:
    '''
    Class: molden:
    Description: Handles functions needed to write molecular orbitals to a molden file.

    molden:
          |orb_set_dict()
          |gen_molden()

    - Function: orb_set_dict()
      Description: Returns a dictionary that is used to know whether [5D], [7F] and [9G]
                   is needed in the molden file.

    - Function: gen_molden()
      Description: Description: Writes molecular orbitals to a molden file.
    '''
    def __init__(self, molden_vectors=None, molden_occs=None,
                    molden_nbas=None, molden_energies=None,
                    molden_basis_set=None,
                    s_orb=False, p_orb=False, d_orb=False,
                    f_orb=False, g_orb=False, h_orb=False):
        self.molden_vectors = None
        self.molden_occs = None
        self.molden_nbas = None
        self.molden_energies = None
        self.molden_basis_set = None
        self.s_orb = False
        self.p_orb = False
        self.d_orb = False
        self.f_orb = False
        self.g_orb = False
        self.h_orb = False

    def orb_set_dict(self):
        '''
        Function: orb_set_dict()
        Description: Returns a dictionary that is used to know whether [5D], [7F] and [9G]
                     is needed in the molden file.
        '''
        orb_set_dict = {
        's': self.s_orb,
        'p': self.p_orb,
        'd': self.d_orb,
        'f': self.f_orb,
        'g': self.g_orb,
        'h': self.h_orb,
        }
        return orb_set_dict

    def gen_molden(self, molden_vectors=None, molden_occs=None, molden_nbas=None,
                   molden_energies=None, label=None, coords=None, atnum=None,
                   orb_set_dict=None, molden_basis_set=None, enocc=False, sym_labels=True,
                   molden_name='m2m'):
        '''
        Function: gen_molden()
        Description: Writes molecular orbitals to a molden file.
        '''
        calc_log().message(screen=True, msg="</> Writing Orbitals to MOLDEN Format.")
        molden_name = molden_name + '.molden'
        _molden_ = []
        _molden_.append('[Molden Format]' + '\n')
        _molden_.append('[Atoms] (AU)' + '\n')

        atom_label = []
        for atom in label:
            atom = atom.decode('UTF-8')
            atom_label.append(atom)

        atom_position = []
        for k in range(len(atnum)):
            k = k + 1
            atom_position.append(k)

        coords = coords.tolist()

        for atom, apos, z, position in zip(atom_label, atom_position, atnum, coords):
            apos = format(apos, '7.0f')
            z = format(z, '7.0f')
            entry = atom + apos + z
            for element in position:
                entry = entry + format(element, '15.7f')
            _molden_.append(entry + '\n')

        if 'd' in orb_set_dict:
            _molden_.append('[5D]' + '\n')
        if 'f' in orb_set_dict:
            _molden_.append('[7F]' + '\n')
        if 'g' in orb_set_dict:
            _molden_.append('[9G]' + '\n')
        _molden_.append('[GTO] (AU)')

        _molden_.append(molden_basis_set)

        _molden_.append('\n' + '[MO]' + '\n')

        with open(molden_name, 'a') as f:
            f.writelines(_molden_)

        sym = []
        for i in range(int(len(molden_nbas))):
            output_sym = "Sym = " + str(i + 1)
            number_of_irreps = molden_nbas[i]
            for k in range(number_of_irreps):
                sym.append(int(i+1))

        if sym_labels == True:
            calc_log().message(screen=True, msg="</> Generating Symmetry Labels.")
            sym_labels = []
            for unique_irrep in np.unique(sym):
                occurance = sym.count(unique_irrep)
                for i in range(occurance):
                    sym_labels.append('(' + str(unique_irrep) + ',' + str(i+1) + ')')
            sym = sym_labels

        occ = []
        for element in molden_occs:
            element = format(element, '11.5f')
            occ.append(element)

        en = []
        for element in molden_energies:
            element = format(element, '11.4f')
            en.append(element)

        neg_occ = []
        for element in molden_occs:
            element = -abs(element)
            element = format(element, '11.5f')
            neg_occ.append(element)

        set_en_2_neg_occ = enocc
        count = -1
        mo_num = 0
        nmos = sum(molden_nbas)
        calc_log().message(screen=True, msg="</> Writing the MOLDEN File.")
        for molecular_orbital in molden_vectors:
            _molden_ = []
            count = count + 1
            percentage = (int(count)/int(nmos)) * 100
            print(str(round(percentage)) + '%', end='\r')
            mo_sym = sym[count]
            mo_en = en[count]
            mo_occ = occ[count]
            neg_mo_occ = neg_occ[count]

            _molden_.append("Sym = " + str(mo_sym) + "\n")
            if set_en_2_neg_occ == False:
                _molden_.append("Ene =" + str(mo_en) + "\n")
            else:
                _molden_.append("Ene =" + str(neg_mo_occ) + "\n")
            _molden_.append("Spin = alpha" + "\n")
            _molden_.append("Occup =" + str(mo_occ) + "\n")
            i = 0
            for element in molecular_orbital:
                i = i + 1
                element = format(element, '19.8f')
                indx = format(i, '4d')
                entry = str(indx) + str(element)
                _molden_.append(entry + "\n")

            with open(molden_name, 'a') as f:
                f.writelines(_molden_)

class h5_data:
    '''
    Class: h5_data:
    Description: Contains a variety of functions useful for writing or extracting data from HDF5 files.

    h5_data:
          |obtain_salcs()
          |obtain_primitives()
          |obtain_primitive_ids()
          |obtain_desym_basis_function_ids()
          |obtain_desym_center_labels()
          |obtain_desym_center_coordinates()
          |obtain_desym_center_atnum()
          |obtain_basis_function_ids()
          |obtain_center_labels()
          |obtain_center_coordinates()
          |obtain_center_atnum()
          |calc_type_indicies()
          |format_density_matrix()
          |update()
          |save_data()
          |read_data()

    - Function: obtain_salcs()
      Description: Reads DESYM_MATRIX data from the HDF5 file, returns as a list.

    - Function: obtain_primitives()
      Description: Reads PRIMITIVES data from the HDF5 file, returns as a list.

    - Function: obtain_primitive_ids()
      Description: Reads PRIMITIVE_IDS data from the HDF5 file, returns as a list.

    - Function: obtain_desym_basis_function_ids()
      Description: Reads DESYM_BASIS_FUNCTION_IDS data from the HDF5 file, returns as a list.

    - Function: obtain_desym_center_labels()
      Description: Reads DESYM_CENTER_LABELS data from the HDF5 file, returns as a list.

    - Function: obtain_desym_center_coordinates()
      Description: Reads DESYM_CENTER_COORDINATES data from the HDF5 file, returns as a list.

    - Function: obtain_desym_center_atnum()
      Description: Reads DESYM_CENTER_ATNUMS data from the HDF5 file, returns as a list.

    - Function: obtain_basis_function_ids()
      Description: Reads BASIS_FUNCTION_IDS data from the HDF5 file, returns as a list.

    - Function: obtain_center_labels()
      Description: Reads CENTER_LABELS data from the HDF5 file, returns as a list.

    - Function: obtain_center_coordinates()
      Description: Reads CENTER_COORDINATES data from the HDF5 file, returns as a list.

    - Function: obtain_center_atnum()
      Description: Reads CENTER_ATNUMS data from the HDF5 file, returns as a list.

    - Function: calc_type_indicies()
      Description: Generates type indicies data based on the molecular orbital occupations.

    - Function: format_density_matrix()
      Description: Formats the density provided to return only the section of the density matrix
                     that contains values that deviate from 2.0. Effectively taking the full density,
                     removing the core, and returning the density that corrisponds to the valence space.

    - Function: update()
      Description: Copies the provided HDF5 file to m2m.h5, removes un-needed data,
                   then creates the following datasets:
                            DENSITY_MATRIX
                            MO_ENERGIES
                            MO_OCCUPATIONS
                            MO_TYPEINDICES
                            MO_VECTORS
                    and stores the density as well as the MO energies, occupations,
                    vectors (MO Coefficents) and type indicies data.

    - Function: save_data()
      Description: Appends supplied contents to a new dataset in the m2m.h5 file.

    - Function: read_data()
      Description: Reads data from a dataset contained in a HDF5 file.
    '''
    def __init__(self, density_matrix = None,
                        mo_energies = None,
                        mo_occupations = None,
                        mo_vectors = None,
                        mo_typeindices = None,
                        center_atnums = None,
                        center_charges = None,
                        center_coordinates = None,
                        desym_matrix = None,
                        desym_basis_function_ids = None,
                        desym_center_charges = None,
                        desym_center_coordinates = None,
                        desym_center_labels = None,
                        primitives = None,
                        primitive_ids = None,
                        ao_overlap_matrix = None
                        ):
        self.density_matrix = None #this will be updated in operations once built.
        self.mo_energies = None
        self.mo_occupations = None
        self.mo_vectors = None
        self.mo_typeindices = None
        self.center_atnums = None
        self.center_charges = None
        self.center_coordinates = None
        self.desym_matrix = None
        self.desym_basis_function_ids = None
        self.desym_center_charges = None
        self.desym_center_coordinates = None
        self.desym_center_labels = None
        self.primitives = None
        self.primitive_ids = None
        self.ao_overlap_matrix = None

    def obtain_salcs(self,InpOrbClass_h5_file,InpOrbClass_nbas):
        '''
        Function: obtain_salcs()
        Description: Reads DESYM_MATRIX data from the HDF5 file, returns as a list.
        '''
        calc_log().message(screen=True, msg="</> Obtaining SALCS.")
        with h5py.File(str(InpOrbClass_h5_file), 'r') as f:
            salcs = f['DESYM_MATRIX']
            n_bas = InpOrbClass_nbas
            n_sym = len(n_bas)
            n_bast = np.sum(n_bas)
            desym_matrix = salcs[:].reshape((n_bast,n_bast), order='F')
            f.close()
            return desym_matrix

    def obtain_primitives(self,InpOrbClass_h5_file):
        '''
        Function: obtain_primitives()
        Description: Reads PRIMITIVES data from the HDF5 file, returns as a list.
        '''
        calc_log().message(screen=True, msg="</> Obtaining PRIMITIVES.")
        with h5py.File(str(InpOrbClass_h5_file), 'r') as f:
            prims = f['PRIMITIVES']

            return prims[:]

    def obtain_primitive_ids(self,InpOrbClass_h5_file):
        '''
        Function: obtain_primitive_ids()
        Description: Reads PRIMITIVE_IDS data from the HDF5 file, returns as a list.
        '''
        calc_log().message(screen=True, msg="</> Obtaining PRIMITIVE IDS.")
        with h5py.File(str(InpOrbClass_h5_file), 'r') as f:
            prim_ids = f['PRIMITIVE_IDS']

            return prim_ids[:]

    def obtain_desym_basis_function_ids(self,InpOrbClass_h5_file):
        '''
        Function: obtain_desym_basis_function_ids()
        Description: Reads DESYM_BASIS_FUNCTION_IDS data from the HDF5 file, returns as a list.
        '''
        calc_log().message(screen=True, msg="</> Obtaining DESYM BASIS FUNCTION IDS.")
        with h5py.File(str(InpOrbClass_h5_file), 'r') as f:
            desym_basis_function_ids = f['DESYM_BASIS_FUNCTION_IDS']

            return desym_basis_function_ids[:]

    def obtain_desym_center_labels(self,InpOrbClass_h5_file):
        '''
        Function: obtain_desym_center_labels()
        Description: Reads DESYM_CENTER_LABELS data from the HDF5 file, returns as a list.
        '''
        calc_log().message(screen=True, msg="</> Obtaining DESYM CENTER LABELS.")
        with h5py.File(str(InpOrbClass_h5_file), 'r') as f:
            desym_center_labels = f['DESYM_CENTER_LABELS']
            label = []
            for center in desym_center_labels:
                label.append(center.split()[0])
            return label

    def obtain_desym_center_coordinates(self,InpOrbClass_h5_file):
        '''
        Function: obtain_desym_center_coordinates()
        Description: Reads DESYM_CENTER_COORDINATES data from the HDF5 file, returns as a list.
        '''
        calc_log().message(screen=True, msg="</> Obtaining DESYM CENTER COORDS.")
        with h5py.File(str(InpOrbClass_h5_file), 'r') as f:
            coord = f['DESYM_CENTER_COORDINATES']
            coord = coord[:]
            return coord

    def obtain_desym_center_atnum(self,InpOrbClass_h5_file):
        '''
        Function: obtain_desym_center_atnum()
        Description: Reads DESYM_CENTER_ATNUMS data from the HDF5 file, returns as a list.
        '''
        calc_log().message(screen=True, msg="</> Obtaining DESYM CENTER ATOM NUMS.")
        with h5py.File(str(InpOrbClass_h5_file), 'r') as f:
            atnum = f['DESYM_CENTER_ATNUMS']
            atnum = atnum[:]
            return atnum

    def obtain_basis_function_ids(self,InpOrbClass_h5_file):
        '''
        Function: obtain_basis_function_ids()
        Description: Reads BASIS_FUNCTION_IDS data from the HDF5 file, returns as a list.
        '''
        calc_log().message(screen=True, msg="</> Obtaining BASIS FUNCTION IDS.")
        with h5py.File(str(InpOrbClass_h5_file), 'r') as f:
            desym_basis_function_ids = f['BASIS_FUNCTION_IDS']

            return desym_basis_function_ids[:]

    def obtain_center_labels(self,InpOrbClass_h5_file):
        '''
        Function: obtain_center_labels()
        Description: Reads CENTER_LABELS data from the HDF5 file, returns as a list.
        '''
        calc_log().message(screen=True, msg="</> Obtaining CENTER LABELS.")
        with h5py.File(str(InpOrbClass_h5_file), 'r') as f:
            desym_center_labels = f['CENTER_LABELS']
            label = []
            for center in desym_center_labels:
                label.append(center.split()[0])
            return label

    def obtain_center_coordinates(self,InpOrbClass_h5_file):
        '''
        Function: obtain_center_coordinates()
        Description: Reads CENTER_COORDINATES data from the HDF5 file, returns as a list.
        '''
        calc_log().message(screen=True, msg="</> Obtaining CENTER COORDS.")
        with h5py.File(str(InpOrbClass_h5_file), 'r') as f:
            coord = f['CENTER_COORDINATES']
            coord = coord[:]
            return coord

    def obtain_center_atnum(self,InpOrbClass_h5_file):
        '''
        Function: obtain_center_atnum()
        Description: Reads CENTER_ATNUMS data from the HDF5 file, returns as a list.
        '''
        calc_log().message(screen=True, msg="</> Obtaining CENTER ATOM NUMS.")
        with h5py.File(str(InpOrbClass_h5_file), 'r') as f:
            atnum = f['CENTER_ATNUMS']
            atnum = atnum[:]
            return atnum

    def calc_type_indicies(self):
        '''
        Function: calc_type_indicies()
        Description: Generates type indicies data based on the molecular orbital occupations.
        '''
        A_vec = self.mo_occupations
        type_indicies = []
        for occ in A_vec:
            if occ == 2.0:
                entry = b'I'
            elif occ == 0.0:
                entry = b'S'
            else:
                entry = b'2'
            type_indicies.append(entry)
        self.mo_typeindices = type_indicies

    def format_density_matrix(self, density):
        '''
        Function: format_density_matrix()
        Description: Formats the density provided to return only the section of the density matrix
                     that contains values that deviate from 2.0. Effectively taking the full density,
                     removing the core, and returning the density that corrisponds to the valence space.
        '''
        import numpy as np
        np.set_printoptions(suppress=False)
        np.set_printoptions(precision=10)
        F_Den = np.delete(density, np.where(density == 2.0), axis=1)
        data = np.asarray(F_Den)
        data = data[~np.all(data == 0, axis=1)]
        nonzero = np.matrix(data)
        self.density_matrix = nonzero

    def update(self, InpOrbClass_h5_file, name, input_dictionary):
        '''
        Function: update()
        Description: Copies the provided HDF5 file to m2m.h5, removes un-needed data,
                     then creates the following datasets:
                            DENSITY_MATRIX
                            MO_ENERGIES
                            MO_OCCUPATIONS
                            MO_TYPEINDICES
                            MO_VECTORS
                    and stores the density as well as the MO energies, occupations,
                    vectors (MO Coefficents) and type indicies data.
        '''
        filename = str(InpOrbClass_h5_file)
        newfile = str(name) + '.h5'
        command = 'cp ' + str(filename) + ' ' + str(newfile)
        os.system(command)
        with h5py.File(newfile, 'a') as f:
            try:
                del f['AO_FOCKINT_MATRIX']
            except:
                print(".", end='\r')
            try:
                del f['AO_MLTPL_X']
            except:
                print("..", end='\r')
            try:
                del f['AO_MLTPL_XX']
            except:
                print(".", end='\r')
            try:
                del f['AO_MLTPL_XY']
            except:
                print("..", end='\r')
            try:
                del f['AO_MLTPL_XZ']
            except:
                print("...", end='\r')
            try:
                del f['AO_MLTPL_Y']
            except:
                print(".", end='\r')
            try:
                del f['AO_MLTPL_YY']
            except:
                print("..", end='\r')
            try:
                del f['AO_MLTPL_YZ']
            except:
                print("...", end='\r')
            try:
                del f['AO_MLTPL_Z']
            except:
                print(".", end='\r')
            try:
                del f['AO_MLTPL_ZZ']
            except:
                print("..", end='\r')
            try:
                del f['AO_OVERLAP_MATRIX']
            except:
                print("...", end='\r')
            try:
                del f['BASIS_FUNCTION_IDS']
            except:
                print(".", end='\r')
            try:
                del f['CENTER_LABELS']
            except:
                print("..", end='\r')
            try:
                del f['CI_VECTORS']
            except:
                print("...", end='\r')
            try:
                del f['DESYM_CENTER_ATNUMS']
            except:
                print(".", end='\r')
            try:
                del f['MLTPL_ORIG']
            except:
                print("..", end='\r')
            try:
                del f['ROOT_ENERGIES']
            except:
                print("...", end='\r')
            try:
                del f['SPINDENSITY_MATRIX']
            except:
                print(".", end='\r')
            try:
                del f['SUPSYM_IRREP_INDICES']
            except:
                print("..", end='\r')
        f = h5py.File(newfile, 'r+')
        MO_Den_SET = []
        MO_Den_SET.append(self.density_matrix)
        density_matrix = MO_Den_SET
        mo_energies = self.mo_energies
        mo_occupations = self.mo_occupations
        mo_typeindices = self.mo_typeindices
        mo_vectors = self.mo_vectors
        try:
            del f['DENSITY_MATRIX']
        except:
            print(".", end='\r')
        try:
            del f['MO_ENERGIES']
        except:
            print("..", end='\r')
        try:
            del f['MO_OCCUPATIONS']
        except:
            print("...", end='\r')
        try:
            del f['MO_TYPEINDICES']
        except:
            print(".", end='\r')
        try:
            del f['MO_VECTORS']
        except:
            print("..", end='\r')
        if input_dictionary['CALC_DENSITY'] == True:
            f.create_dataset('DENSITY_MATRIX', data=density_matrix)
        f.create_dataset('MO_ENERGIES', data=mo_energies)
        f.create_dataset('MO_OCCUPATIONS', data=mo_occupations)
        f.create_dataset('MO_TYPEINDICES', data=mo_typeindices)
        f.create_dataset('MO_VECTORS', data=mo_vectors)

    def save_data(self, name, key, contents):
        '''
        Function: save_data()
        Description: Appends supplied contents to a new dataset in the m2m.h5 file.
        '''
        filename = str(name) + '.h5'
        f = h5py.File(filename, 'a')
        f.create_dataset(key, data=contents)
        f.close()

    def read_data(self, name, key):
        '''
        Function: read_data()
        Description: Reads data from a dataset contained in a HDF5 file.
        '''
        filename = str(name) + '.h5'
        f = h5py.File(filename, 'r')
        dset = f.get(key)[:]
        f.close()
        return dset


class ToolBox:
    '''
    Class: ToolBox:
    Description: This class manages information either

    ToolBox:
            |format_number()
            |form_matricies()

    - Function: format_number()
      Description: Converts number to scientific notation.

    - Function: form_matricies()
      Description: Reads lists of coefficents for each MO from the MO_LIST_X datasets
                   contained within the m2m.h5 file. Then generates a coefficent matrix for all MO's.

    - Function: density_from_momatrix()
      Description: Takes the coefficent matrix and generates the density matrix in the AO basis.

    - Function: MO_Dens()
      Description: Converts the provided AO density to the MO basis.

    - Function: MO_Dens_SYM()
      Description: Converts the provided AO density to the MO basis.

    - Function: density_as_square()
      Description: Takes the provided density, which contains only one half of the values on one side of the diagonal, and generates
                   a square density.

                     Provided:
                     ---------
                     d11
                     d21 d21
                     d31 d32 d33

                     After conversion:
                     -----------------
                     d11 (d12=d21) (d13=d31)
                     d21 d22 (d23=d32)
                     d31 d32 d33

    - Function: reshape_square()
      Description:
        This code defines a method reshape_square for an object self that takes two arguments,
        an array arr and a list of dimensions dims.
        The method first checks if the length of dims is 1, if so, it assigns the first element
        of dims to dim and reshape the array arr in to a square matrix with dimensions dim*dim using 'F' order.

        Otherwise, it creates an empty list lst and assigns the value 0 to offset.
        Then it enters a for loop, iterating over the elements of dims.
        For each iteration, it takes a slice of arr starting at index offset and ending at index offset+dim^2.
        It reshapes this slice into a square matrix with dimensions dim*dim and using 'F' order, and appends it to lst.
        Finally, it increases the value of offset by dim^2.
        The method returns the block diagonal matrix made from the matrices in lst using numpy's block_diag method.

    - Function: form_matricies_sym()
      Description: Reads lists of coefficents for each MO from the MO_LIST_X datasets
                   contained within the m2m.h5 file. Then generates a coefficent matrix for all MO's.
    '''
    def format_number(self,number):
        '''
        Function: format_number()
        Description: Converts number to scientific notation.
        '''
        number = "{:e}".format(number)
        return number

    def form_matricies(self,InpOrbClass_nMOs):
        '''
        Function: form_matricies()
        Description: Reads lists of coefficents for each MO from the MO_LIST_X datasets
                     contained within the m2m.h5 file. Then generates a coefficent matrix for all MO's.
        '''
        calc_log().message(screen=True, msg="</> Building Coeff matrix.")
        dset_byte_array = h5_data().read_data(name='m2m_DATA', key='MO_LIST_1')
        dset_string_list = [s.decode().strip() for s in dset_byte_array]
        coeffs = dset_string_list
        nMO_coeffs = int(len(coeffs))
        InpOrbClass_nMO_coeffs = nMO_coeffs
        nMOs = int(InpOrbClass_nMOs)
        Cmat = np.zeros((nMOs,nMOs))

        for i in range(int(InpOrbClass_nMOs)):
            file_name = "MO_LIST_" + str(i+1)
            dset_byte_array = h5_data().read_data(name='m2m_DATA', key=file_name)
            dset_string_list = [s.decode().strip() for s in dset_byte_array]
            coeffs = dset_string_list
            Cmat[:,i] = coeffs[:]
        return Cmat

    def density_from_momatrix(self,cmat, occvec):
        '''
        Function: density_from_momatrix()
        Description: Takes the coefficent matrix and generates the density matrix in the AO basis.
        '''
        calc_log().message(screen=True, msg="</> Bulding Density in AO Basis.")
        nbas = len(occvec)
        arlen = nbas * (nbas + 1) // 2
        dens = np.empty(arlen, dtype=np.float64)
        cnt = 0
        for i in range(nbas):
            for j in range(i + 1):
                dens[cnt] = (cmat[i,:] * cmat[j,:] * occvec).sum()
                cnt += 1
        return dens

    def MO_Dens(self,cmat, dens, ao_overlap,occvec):
        '''
        Function: MO_Dens()
        Description: Converts the provided AO density to the MO basis.
        '''
        calc_log().message(screen=True, msg="</> Transforming the AO Density to MO Basis.")
        Sao = np.asmatrix(ao_overlap)
        Pao = np.asmatrix(dens)
        Cmo = np.asmatrix(cmat)
        Vt = (Sao * Cmo).T
        V = Sao * Cmo
        D_mo = Vt * Pao * V
        D_mo = D_mo[:,occvec>0]
        D_mo = D_mo.round(decimals=9)
        np.set_printoptions(suppress=True)
        np.set_printoptions(precision=10)
        A_len = int(len(occvec[occvec>0]))
        F_Den = np.ones((A_len, A_len))
        for i in range(A_len):
        	for j in range(A_len):
        		F_Den[i,j] = D_mo[i,j]
        return F_Den

    def MO_Dens_SYM(self,cmat, dens, ao_overlap,occvec):
        '''
        Function: MO_Dens_SYM()
        Description: Converts the provided AO density to the MO basis.
        '''
        calc_log().message(screen=True, msg="</> Transforming the AO Density to MO Basis.")
        Sao = np.asmatrix(ao_overlap)
        Pao = np.asmatrix(dens)
        Cmo = np.asmatrix(cmat)
        Vt = (Sao * Cmo).T
        V = Sao * Cmo
        D_mo = Vt * Pao * V
        D_mo = D_mo[:,occvec>0]
        D_mo = D_mo.round(decimals=9)
        np.set_printoptions(suppress=True)
        np.set_printoptions(precision=10)
        nonzero = D_mo[~np.all(D_mo == 0, axis=1)]
        F_Den = nonzero
        F_Den = np.delete(F_Den, np.where(F_Den == 2.0), axis=1)
        return D_mo, nonzero

    def density_as_square(self,denvec):
        '''
        Function: density_as_square()
        Description: Takes the provided density, which contains only one half of the values on one side of the diagonal, and generates
                     a square density.

                     Provided:
                     ---------
                     d11
                     d21 d21
                     d31 d32 d33

                     After conversion:
                     -----------------
                     d11 (d12=d21) (d13=d31)
                     d21 d22 (d23=d32)
                     d31 d32 d33
        '''
        calc_log().message(screen=True, msg="</> Converting triangle density to square.")
        nbas = int((-1 + np.sqrt(1 - 4 * -2 * len(denvec))) / 2)
        square = np.empty((nbas, nbas), dtype=np.float64)
        cnt = 0
        for i in range(nbas):
            for j in range(i + 1):
                square[i, j] = denvec[cnt]
                square[j, i] = denvec[cnt]
                cnt += 1
        return square

    def reshape_square(self,arr, dims):
        '''

        NOTE: This function was taken from: https://github.com/steabert/molpy/blob/master/molpy/wfn.py (2023)
        An open source project under the terms of the GNU General Public License v2.0.
        Listed Molpy contributors are: Steven Vancoillie, Felix Plasser and Marcus Johansson.

        Function: reshape_square()
        Description:
        This code defines a method reshape_square for an object self that takes two arguments,
        an array arr and a list of dimensions dims.
        The method first checks if the length of dims is 1, if so, it assigns the first element
        of dims to dim and reshape the array arr in to a square matrix with dimensions dim*dim using 'F' order.

        Otherwise, it creates an empty list lst and assigns the value 0 to offset.
        Then it enters a for loop, iterating over the elements of dims.
        For each iteration, it takes a slice of arr starting at index offset and ending at index offset+dim^2.
        It reshapes this slice into a square matrix with dimensions dim*dim and using 'F' order, and appends it to lst.
        Finally, it increases the value of offset by dim^2.
        The method returns the block diagonal matrix made from the matrices in lst using numpy's block_diag method.
        '''
        if len(dims) == 1:
            dim = dims[0]
            return arr.reshape((dims[0],dims[0]), order='F')
        lst = []
        offset = 0
        for dim in dims:
            slice_ = arr[offset:offset+dim**2]
            lst.append(slice_.reshape((dim,dim), order='F'))
            offset += dim**2
        return la.block_diag(*lst)

    def form_matricies_sym(self,InpOrbClass_nMOs,InpOrbClass_nbas):
        '''
        Function: form_matricies_sym()
        Description: Reads lists of coefficents for each MO from the MO_LIST_X datasets
                     contained within the m2m.h5 file. Then generates a coefficent matrix for all MO's.
        '''
        calc_log().message(screen=True, msg="</> Building COEFF Matrix.")
        MO_LIST = []
        for i in range(int(InpOrbClass_nMOs)):
            file_name = "MO_LIST_" + str(i+1)
            dset_byte_array = h5_data().read_data(name='m2m_DATA', key=file_name)
            dset_string_list = [s.decode().strip() for s in dset_byte_array]
            coeffs = dset_string_list
            for coeff in coeffs:
                MO_LIST.append(coeff)
        MO_ARRAY = np.asarray(MO_LIST)
        MO_ARRAY = MO_ARRAY.astype(np.float)
        h5_class_mo_vectors = MO_ARRAY
        n_bas = InpOrbClass_nbas
        MO_ARRAY = self.reshape_square(MO_ARRAY, n_bas)
        return MO_ARRAY, h5_class_mo_vectors
