# Calculate_LIPID_parameters_using_Voronoi_Method
#########################################################################################################
###### README file for Usage of "calc_lipid_parameters_voronoi.py" script ###############################
#########################################################################################################

#################################################################
####################### DESCRIPTION #############################
#################################################################

Python Script to Calculate LIPID parameters using Voronoi tessellation and Monte Carlo integration methods.

By: Prajwal Nandekar
Email: prajwal.pharm07@gmail.com
Updated on: 1st May 516

Command:
python calc_lipid_parameters_voronoi.py -i frame -start 1 -stop 5

Requirements:

Python: Scientific Python packages

Required modules
########################
import sys,os
import itertools
from scipy import stats
import numpy as np
import argparse
import subprocess
from tabulate import tabulate
#######################

Sample Input files:

Directory: ~/Examples

Description:

This script calculate the LIPID parameters using Voronoi tessellation and Monte Carlo integration methods.
 Prior to running this script an you need PDB files exatracted from MD trajectory in sequences.
For example: frame_0.pdb, frame_1.pdb ... frame_5.pdb.

Essential Note:
1) Make sure the box information should be present as header line in PDB files.
You can used VMD to save PDBs from MD trajectory.
Sample VMD script is in Examples directory 'stride_and_save_frames.vmd'
2) Make sure the "vtmc" binary file is in working directory

        :Arguments:
            *pdbfile*
                Input prefix of PDB filenames. (Example: "-i frame", "if PDB filenames are frame_1.pdb").
                            *start_frame*
                Starting frame number. (Example: "-start 1", "If PDB filenames sequence starts from frame_1.pdb")
            *stop_frame*
                Last frame number . (Example: "-stop 5", "If last PDB filenames in sequence is frame_5.pdb").

        :Returns:
            *AREA PER LIPID VALUES FOR TOP LAYER, BOTTOM LAYER, AND AVERAGE VALUES FOR ALL LIPIDS, BOUND LIPIDS< AND UNBOUND LIPIDS.*

Command:
        python calc_lipid_parameters_voronoi.py -i prefix -start frame_number -stop frame_number

#################################################################
########################### E N D ###############################
#################################################################

