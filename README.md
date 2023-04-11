# Molcas2Molden
Python3 code to convert Molcas INPORB files to MOLDEN format.

# Introduction:
This code was written to convert SOOrb files generated in RASSI calculations to MOLDEN files. SOOrb files contain the orbital coefficients for the spin-orbit natural orbitals for spin-orbit coupled states. The code takes these coefficients and basis set information from the HDF5 file to generate the MOLDEN file. The code was only written to handle the conversion of SOOrb files, which contain the orbital coefficients in a 'restricted' format (only a single Alpha set of orbitals). As such the code will work with converting other restricted orbital sets contained within RasOrb and ScfOrb files. In the future, additional features will be added such as the ability to generate MOLDEN files that contain unrestricted set of orbitals (Alpha and Beta set).

# Installation:

Simply download all the files and unzip to an appropriate location.

# Quick start:
The program can convert either single or multiple files during a run.

For single file runs, ensure you have both the INPORB and HDF5 file from your molcas calculation in a directory. Then run the m2m.py script:
'''
python3 m2m.py
'''
The program will detect the files and run the conversion.

For multiple file runs, ensure you have all the INPORB files and the single HDF5 file from your molcas calculation in a directory. Then run the m2m.py script:
'''
python3 m2m.py
'''
The program will detect all the files and run the conversion for each one at a time. Once a SOOrb.x file is converted to MOLDEN format, 'STATE_X.molden', it is stored in a new directory called MOLDENS.

# Acknowledgments:

Thanks goes to Dr Andrew Kerridge for the helpful discussions and guidance.

I'd like to acknowledge the molpy program, https://github.com/felixplasser/molpy, for inspiration to write this conversion tool. Some code was utilized from this program and this is indicated within scripts when done.
