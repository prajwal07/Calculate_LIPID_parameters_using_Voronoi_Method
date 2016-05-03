#!/sw/mcm/app/anaconda/2.3.0/envs/mdaenv/bin/python
import sys,os
import itertools
from scipy import stats
import numpy as np
import argparse
import subprocess
from tabulate import tabulate

######################################################

parser = argparse.ArgumentParser(description='Calculate LIPID parameters using Voronoi tessellation and Monte Carlo integration methods. \n Command: "python calc_lipid_parameters_voronoi.py -i prefix -start frame_number -stop frame_number" ')
parser.add_argument('-i', '--input1', dest = "pdbfile_prefix", help='Input prefix of PDB filenames. (Example: "-i frame", "if PDB filenames are frame_1.pdb").')
parser.add_argument('-start', '--input2', dest = "start_frame", help='Starting frame number. (Example: "-start 1", "If PDB filenames sequence starts from frame_1.pdb").')
parser.add_argument('-stop', '--input3', dest = "stop_frame", help='Last frame number . (Example: "-stop 5", "If last PDB filenames in sequence is frame_5.pdb").')
args = parser.parse_args()

#####################################################
"""
Application:

Python Script to Calculate LIPID parameters using Voronoi tessellation and Monte Carlo integration methods.

Description:

This script calculate the LIPID parameters using Voronoi tessellation and Monte Carlo integration methods. 
Prior to running this script an you need PDB files exatracted from MD trajectory in sequences. 
For example: frame_0.pdb, frame_1.pdb ... frame_5.pdb.
Make sure the box information should be present as header line in PDB files. (You can used VMD to save PDBs from MD trajectory.

Essential Note:
1) Make sure the box information should be present as header line in PDB files.
You can used VMD to save PDBs from MD trajectory.
2) Make sure the "vtmc" binary file is in working directory

        :Arguments:
            *pdbfile*
                Input prefix of PDB filenames. (Example: "-i frame", "if PDB filenames are frame_1.pdb").
            *start_frame*
                Starting frame number. (Example: "-start 1", "If PDB filenames sequence starts from frame_1.pdb")                
            *stop_frame*
                Last frame number . (Example: "-stop 20", "If last PDB filenames in sequence is frame_5.pdb").
            
        :Returns:
            *AREA PER LIPID VALUES FOR TOP LAYER, BOTTOM LAYER, AND AVERAGE VALUES FOR ALL LIPIDS, BOUND LIPIDS AND UNBOUND LIPIDS.*

Command:
        python calc_lipid_parameters_voronoi.py -i prefix -start frame_number -stop frame_number
"""
#####################################################

PATH = os.getcwd()
print '\nCurrent working directory is: ', PATH
print ""
os.system("rm -rf vtm-analysis")
os.system("mkdir vtm-analysis")
print "\nVTM analysis Directory created \n"

global xdim
global ydim
global zdim
global output_PDB_file_name
global lipid_name
global head_atom
input_file_name_prefix = args.pdbfile_prefix
start_frame_num = int(args.start_frame)
last_frame_num = int(args.stop_frame)

##################################################################################################
def modify_pdbs(frame):
#    DESCRIPTION
#         Python script to generate input file for VTM area per lipid calculator tool for calculation of area per lipid
#    ARGUMENTS
#        It takes total 1 parameters as input
#        1.) Name of PDB file of protein-membrane system
#    USAGE
#         python prepare_pdb_for_VTM_area_per_lipid.py pdb_file
#         pdb_file -   Name of the PDB file for the protein or ligand with .pdb extension
#    EXAMPLE
#        python prepare_pdb_for_VTM_area_per_lipid.py Rivaroxaban.pdb
    input_PDB_file_name = ('%s' % input_file_name_prefix) + "_" + ('%d' % frame) + ".pdb"
#    print input_PDB_file_name
    input_PDB_file = open(input_PDB_file_name, 'r').readlines()

####################### FOR LIPID 14 ############################
    lipid14_flag = 0
    find_lipid14 = " OL "
    global lipid_name
    global head_atom
    lipid_name = "POC"
    head_atom = "P"
    preprocess_filename = "temp_preprocess_lipid14.sh"
    script_data = """#!/bin/sh
##inputfile=%s_%d.pdb
i=%d
file=%s_${i}.pdb
echo $file
output=temp_frame_ren_${i}.pdb
l=`fgrep -n PA $file  | cut -d':' -f 1`
tmp=
cmd=sed
for i in $l; do
        if [ "${i}" != "$((tmp+1))" ] && [ ! -z "$tmp" ]; then
                n=`sed "${tmp}q;d" $file | cut -c23-27`
                cmd="${cmd} -e \\"s/X\\\(\\\ \\\?\\\)$((n+1))/X\\\\\\1$((n))/g\\""
                cmd="${cmd} -e \\"s/X\\\(\\\ \\\?\\\)$((n+2))/X\\\\\\1$((n))/g\\""
        fi
        tmp=$i
done
n=`sed "${tmp}q;d" $file | cut -c23-27`
cmd="${cmd} -e \\"s/X\\\(\\\ \\\?\\\)$((n+1))/X\\\\\\1$((n))/g\\""
cmd="${cmd} -e \\"s/X\\\(\\\ \\\?\\\)$((n+2))/X\\\\\\1$((n))/g\\""

echo Build

echo $cmd $file > temp_comm.sh
chmod +x temp_comm.sh

echo Running ... Please wait ... 

./temp_comm.sh > $output
rm temp_comm.sh
echo "Rename output file"
mv $output $file
echo Done

""" % (input_file_name_prefix, frame, frame, input_file_name_prefix)

    for line in input_PDB_file:
        if find_lipid14 in line:
            print "The input PDBs consist of LIPID molecules for LIPID14 force fields."
            print "Now, generating and running the script to preprocess the PDBs for LIPID14: %s" % preprocess_filename
            with open(preprocess_filename, 'w') as f:
                f.write(script_data)
                lipid14_flag = 1
                break

    if (lipid14_flag == 1):
        preprocess = "chmod +x " + ('%s' % preprocess_filename)
        os.system (preprocess)
        run_preprocess = "./" + ('%s' % preprocess_filename)
        os.system (run_preprocess)
        lipid_name = "PC"
        head_atom = "P31"
    
################################################################################    
    
    
    aa_name = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","HIE","HIP","HID","GLH","ASH","HEC"]
    ligand = ["PA","OL","PC","POC"] # Change Ligand Names

    output_PDB_file = input_PDB_file

    for i in range(len(output_PDB_file)):
        line = output_PDB_file[i]
        residue_name = line[17:20].strip()
        if line.startswith("ATOM"):
            if (residue_name in aa_name):
                new_line = line[:72] + "PRO1" + line[75:]
                output_PDB_file[i] = new_line

            if (residue_name in ligand):
                new_line = line[:72] + "LIP1" + line[75:]
                output_PDB_file[i] = new_line
    global output_PDB_file_name
    output_PDB_file_name = "temp_" + input_PDB_file_name
    open(output_PDB_file_name, 'w').writelines(output_PDB_file)

###################################################################################################

def generate_conf_file(frame):
    input_PDB_file_name = ('%s' % input_file_name_prefix) + "_" + ('%d' % frame) + ".pdb"
#    print input_PDB_file_name
    with open(input_PDB_file_name) as f:
        for line in itertools.islice(f, 0, 1):  # start=0 (1st line), stop=1 (2nd line), This the line contains BOX DIMENSIONS
            split_line=line.split()
            global xdim
            global ydim
            global zdim
            xdim = float(split_line[1]) # X-dimension
            ydim = float(split_line[2]) # Y-dimension
            zdim = float(split_line[3]) # Z-dimension
            grid_dim_line = ('%3f' % xdim) + "  " + ('%3f' % ydim) + "\n"
    conf_data = """# configulation file for POPC 
## system size(x,y,z)  = (%f, %f, %f)

[INPUT]
INPUTPDB  = %s

[OUTPUT]
OUTNAME   = vtm-analysis/%s_%d
OUTPUTLEV = 2

[PROTEIN]
SELECT    = ALA ARG ASN ASP CYS GLN GLU GLY HIS &
            ILE LEU LYS MET PHE PRO SER THR TRP &
            TYR VAL HIE HIP HID GLH ASH HEC ALA
            ARG ASN ASP CYS GLN GLU GLY

HYDROGEN  = INCLUDE

[LIPID]
SELECT    = %s
HEAD_ATOM = %s
HYDROGEN  = INCLUDE

[ATOM_MASS]
H =   1.008
C =  12.01
O =  16.00
N =  14.01
P =  30.97
S =  32.07

[VDW_RADIUS]
H =  1.2
C =  1.7
O =  1.52
N =  1.55
P =  1.8
S =  1.8

[ANALYSIS]
CELL_CENTER =  0.0     0.0
CELL_SIZE   =  %f  %f
RADIUS      =  0.1
FACTOR      =  100
ISEED       = 3141592
RANDOMIZE   = NO
DECIMAL_PLACE =  12
""" % (xdim, ydim, xdim, output_PDB_file_name, input_file_name_prefix, frame, lipid_name, head_atom, xdim, ydim)
    
    conf_file_name = "temp_" + ('%d' % frame) + "_vtmc.conf"
    conf_log_file_name = "temp" + ('%d' % frame) + "_vtmc.log"
    with open(conf_file_name, 'w') as fconf:
        fconf.write(conf_data)
                                
    vtmc_cmd = "./vtmc " + ('%s' % conf_file_name) + " > " + ('%s' % conf_log_file_name)
    os.system(vtmc_cmd)

###################################################################################################

####################################################################################################
output_data_file_name = "LIPIDS_PROPERTIES.DAT"
output_data_file = open(output_data_file_name, 'w')
output_data_file_lines = []

output_header2 = "FILENAME | X_DIMENSION | Y_DIMENSION | TOP_ALL_LIPIDS | TOP_APL_ALL | TOP_BOUND_LIPIDS | TOP_APL_BOUND | TOP_UNBOUND_LIPIDS | TOP_APL_UNBOUND | DOWN_ALL_LIPIDS | DOWN_APL_ALL | DOWN_BOUND_LIPIDS | DOWN_APL_BOUND | DOWN_UNBOUND_LIPIDS | DOWN_APL_UNBOUND | ALL_LIPIDS | AVG_APL_ALL | ALL_BOUND_LIPIDS | AVG_APL_BOUND | ALL_UNBOUND_LIPIDS | AVG_APL_UNBOUND " + "\n"

output_data_file_lines.append(output_header2)

#print output_header1
print output_header2


def extract_data(frame):

    input_data_filename = "temp" + ('%d' % frame) + "_vtmc.log"

    find_layer1 = "Layer 1"
    find_layer2 = "Layer 2"
    find_summary = "Summary of the results"
    flag=0

####################################################################################################

    with open(input_data_filename, 'r') as input_data_file_lines:
        for line in input_data_file_lines:
            if find_layer1 in line:
                flag=1
            elif find_layer2 in line:
                break
            elif (flag==1):
                find = "RESULTS"
                if find in line:
                    line1 = input_data_file_lines.next()
                    line2 = input_data_file_lines.next()
                    split_line2=line2.split()
                    num_lipids_all = (split_line2[2])
                    num_lipids_bound = (split_line2[3])
                    num_lipids_unbound = (split_line2[4])
                    line3 = input_data_file_lines.next()
                    split_line3=line3.split()
                    apl_lipids_all = split_line3[4]
                    apl_lipids_bound = split_line3[5]
                    apl_lipids_unbound = split_line3[6]
        
                    up_output = ('%s' % num_lipids_all) + " | " + ('%s' % apl_lipids_all) + " | " + ('%s' % num_lipids_bound) + " | " + ('%s' % apl_lipids_bound) + " | " + ('%s' % num_lipids_unbound) + " | " + ('%s' % apl_lipids_unbound)
                    
####################################################################################################

    with open(input_data_filename, 'r') as input_data_file_lines:
        for line in input_data_file_lines:
            if find_layer2 in line:
                flag=2
            elif find_summary in line:
                break
            elif (flag==2):
                find = "RESULTS"
                if find in line:
                    line1 = input_data_file_lines.next()
                    line2 = input_data_file_lines.next()
                    split_line2=line2.split()
                    num_lipids_all = (split_line2[2])
                    num_lipids_bound = (split_line2[3])
                    num_lipids_unbound = (split_line2[4])
                
                    line3 = input_data_file_lines.next()
                    split_line3=line3.split()
                    apl_lipids_all = split_line3[4]
                    apl_lipids_bound = split_line3[5]
                    apl_lipids_unbound = split_line3[6]
                
                    down_output = ('%s' % num_lipids_all) + " | " + ('%s' % apl_lipids_all) + " | " + ('%s' % num_lipids_bound) + " | " + ('%s' % apl_lipids_bound) + " | " + ('%s' % num_lipids_unbound) + " | " + ('%s' % apl_lipids_unbound)

####################################################################################################

    with open(input_data_filename, 'r') as input_data_file_lines:
        for line in input_data_file_lines:
            if find_summary in line:
                line1 = input_data_file_lines.next()
                line2 = input_data_file_lines.next()
                split_line2=line2.split()
                avg_apl_lipids_all = (split_line2[6])
                num_lipids_all = (split_line2[8])

                line3 = input_data_file_lines.next()
                split_line3=line3.split()
                avg_apl_lipids_bound = (split_line3[6])
                num_lipids_bound = (split_line3[8])

                line4 = input_data_file_lines.next()
                split_line4=line4.split()
                avg_apl_lipids_unbound = (split_line4[6])
                num_lipids_unbound = (split_line4[8])

                avg_output = ('%s' % num_lipids_all) + " | " + ('%s' % avg_apl_lipids_all) + " | " + ('%s' % num_lipids_bound) + " | " + ('%s' % avg_apl_lipids_bound) + " | " + ('%s' % num_lipids_unbound) + " | " + ('%s' % avg_apl_lipids_unbound)

####################################################################################################

    output_data = "frame_" + ('%d' % frame) + " | " + ('%3f' % xdim) + " | " + ('%3f' % ydim) + " | " + ('%s' % up_output)  + " | " + ('%s' % down_output) + " | " + ('%s' % avg_output) + "\n"
    print output_data
    output_data_file_lines.append(output_data)

####################################################################################################

for frame in range(start_frame_num, last_frame_num+1):
    modify_pdbs(frame)
    generate_conf_file(frame)
    extract_data(frame)

output_data_file.writelines(output_data_file_lines)
remove_temp = "rm -rf temp*"
os.system(remove_temp)
print "Please find the OUTPUT data is in file '%s'. \n" % (output_data_file_name)

output_data_file.close()
 
