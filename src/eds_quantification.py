"""
Description:
    This script quantifies the Electron Density Shift (EDS) calculated with 
    Gaussian and contained in a cube file.

    The EDS is calculated using the eds_gaussian.sh script

Dependencies:
    numpy
    pandas
    glob
    itertools

Usage:
    python3 eds_quatification -f <filename> -c <cutoff> -a

Arguments:
    <filename> (str): Path to the EDS cube file.
    <cutoff> (float): Cutoff for quantifying EDS in
    <a> (boolean): Select for quantify all the cub files in the working directory

Author information:
"""

__author__ = 'Inigo Iribarren'
__email__ = 'iribarri@tcd.ie'
__version__ = '1.0'
__license__ = 'MIT'

# Import necessary libraries
import argparse
import os
import numpy as np
import pandas as pd
import glob
import itertools
from eds_cubetools import *


def cube_parser(cube_file):
    """
    Extracts all the needed information to quantify the EDS form a 
    selected cube file 

    Args:
        cube_file (str): Name of the cube file to extract the data from 

    Returns:
        grid_df (pandas.DataFrame): dataframe with all the grid coordinatees and rho values 
        cub_vol (float): volume of individual cube forming the grid
    """

    # Read the cube file and storage the date in two different ways
    # Store EDS at each point as 1 dimmension array
    rho, meta_data, nx, ny, nz = read_cube(cube_file)

    # Extracting cage and grid size
    org_list = list(meta_data['org']) 
    xorg, yorg, zorg = org_list[0], org_list[1], org_list[2] # Origin of grid

    xstep = list(meta_data['xvec'])[0] 
    ystep = list(meta_data['yvec'])[1] 
    zstep = list(meta_data['zvec'])[2] 

    cub_vol = xstep * ystep * zstep # Volume of a single cube in the grid

    # Coordinates of all the grid points and their density
    xcoord = np.arange(xorg, xorg + xstep * nx, xstep)
    ycoord = np.arange(yorg, yorg + ystep * ny, ystep)
    zcoord = np.arange(zorg, zorg + zstep * nz, zstep)

    # This is needed to delete extra points that can be created at the end due 
    # to accuracy issues 
    if ( xcoord.size != nx ):
        xcoord = xcoord[:nx-xcoord.size]
    if ( ycoord.size != ny ):
        ycoord = ycoord[:ny-ycoord.size]
    if ( zcoord.size != nz ):
        zcoord = zcoord[:nz-zcoord.size]
    
    # Extract the grid points coordinates and the rho value
    grid = np.asarray([[x, y, z, rho[i]] for i, (x, y, z) 
            in enumerate(itertools.product(xcoord, ycoord, zcoord))])
    grid_df = pd.DataFrame(grid, columns=['x','y','z','rho'])

    return grid_df, cub_vol

def calc_eds_total(rho_df, cub_vol):
    """
    Calculates the total value of the EDS in the file

    Args:
        rho_df (pandas.DataFrame): Pandas DataFrame containing the grid points
        coordinates and the rho values. Columns:
            - x (float): X coordinates
            - y (float): Y coordinates
            - z (float): Z coordinates
            - rho (float): EDS values at the gridpoint
        cub_vol (float): Volume of each individual cube in the grid

    Returns:
        eds_total_pos (float): total value of all the positive EDS points
        eds_total_neg (float): total value of all the positive EDS points
        eds_total (float): total value of all the EDS points
        eds_error (float): error between EDS pos and neg
    """

    eds_total_pos = rho_df[rho_df['rho'] > 0]['rho'].to_numpy().cumsum()[-1] * cub_vol
    eds_total_neg = rho_df[rho_df['rho'] < 0]['rho'].to_numpy().cumsum()[-1] * cub_vol
    eds_total = rho_df['rho'].to_numpy().cumsum()[-1] * cub_vol
    eds_error = abs(eds_total / ( eds_total_pos - eds_total_neg ) * 100)

    return eds_total_pos, eds_total_neg, eds_total, eds_error

def calc_eds_cutoff(rho_df, cub_vol, cutoff):
    """
    Calculates the value of the EDS in the file at a certain cutoff

    Args:
        rho_df (pandas.DataFrame): Pandas DataFrame containing the grid points
        coordinates and the rho values. Columns:
            - x (float): X coordinates
            - y (float): Y coordinates
            - z (float): Z coordinates
            - rho (float): EDS values at the gridpoint
        cub_vol (float): Volume of each individual cube in the grid
        cutoff (float): Cutoff for the EDS calculations.

    Returns:
        eds_cutoff_pos (float): cutoff value of all the positive EDS points
        eds_cutoff_neg (float): cutoff value of all the positive EDS points
    """

    eds_cutoff_pos = rho_df[rho_df['rho'] > cutoff]['rho'].to_numpy().cumsum()[-1] * cub_vol
    eds_cutoff_neg = rho_df[rho_df['rho'] < cutoff]['rho'].to_numpy().cumsum()[-1] * cub_vol

    return eds_cutoff_pos, eds_cutoff_neg

# Multiple EDS 
def multiple_eds(cutoff):
    """
    Quantifies the EDS of all the <name>_final.cub files in the working directory

    Args:
        cutoff (float): Cutoff for the EDS calculations.

    Returns:
        df_eds (pandas.Dataframe): Dataframe with all EDS data from the files in the working directory
    """

    files = glob.glob('./*final.cub')
    files.sort()

    df_columns = [
            'File_name',
            'EDS_total_positive',
            'EDS_total_negative',
            'EDS_total',
            'EDS_error',
            'EDS_cutoff_positive',
            'EDS_cutoff_negative'
            ]
    eds_data = []

    for file in files:
        name = file.split('/')[-1].split('.')[0]

        eds_total_pos, eds_total_neg, eds_total, eds_error, eds_cutoff_pos, eds_cutoff_neg = main(file,cutoff)
        
        eds_data.append([name, eds_total_pos, eds_total_neg, eds_total, eds_error, eds_cutoff_pos, eds_cutoff_neg])

    df_eds = pd.DataFrame(data=eds_data, columns=df_columns)

    return df_eds

# Single EDS function
def main(cube_file, cutoff):
    """
    Quantifies the EDS of a selected cub file in the working directory

    Args:
        cube_file (str): Name of the cub file to calculate EDS from
        cutoff (float): Cutoff for the EDS calculations.

    Returns:
        df_eds (pandas.Dataframe): Dataframe with all EDS data from the files in the working directory
    """
    rho_data, cub_vol = cube_parser(cube_file)

    eds_total_pos, eds_total_neg, eds_total, eds_error = calc_eds_total(rho_data,cub_vol)
    eds_cutoff_pos, eds_cutoff_neg = calc_eds_cutoff(rho_data,cub_vol,cutoff)

    return eds_total_pos, eds_total_neg, eds_total, eds_error, eds_cutoff_pos, eds_cutoff_neg


# Call the main function
if __name__ == "__main__":
    #  Parse arguments
    parser = argparse.ArgumentParser(
            description='This script quantifies the Electron Density Shift \
                    (EDS) calculated with Gaussian and contained in a cube \
                    file. The EDS is calculated using the eds_gaussian.sh script')

    parser.add_argument('-f','--filename',
                        type=str,
                        help='Filename of the cube file containing the EDS')
    parser.add_argument('-c','--cutoff',
                        type=float,
                        help='Cutoff for EDS value',
                        default=0.001)
    parser.add_argument('-a','--all',
                        action='store_true',
                        help='Use it if you want to caculate the EDS of all the\
                            cube files in the current directory.\
                            The files must be called <name>_final.cub',
                        default=False)

    args = parser.parse_args()

    filename = args.filename
    cutoff = args.cutoff
    eds_all = args.all
    
    if not eds_all and filename == None:
        raise ValueError('Filename has to be specified or "-a" option selected ')
    elif eds_all and filename != None:
        raise ValueError('Filename cannot be specified when "-a" option is selected ')
    elif not eds_all and filename != None:
        if not os.path.exists(filename):
            raise FileNotFoundError(f"The file '{filename}' does not exist")
        else:
            eds_total_pos, eds_total_neg, eds_total, eds_error, eds_cutoff_pos, eds_cutoff_neg = main(filename,cutoff)
            print(f'EDS total positive: {eds_total_pos:.6f}')
            print(f'EDS total negative: {eds_total_neg:.6f}')
            print(f'EDS total: {eds_total}')
            print(f'EDS total error: {eds_error:.2f}%')
            print(f'EDS cutoff positive: {eds_cutoff_pos:.6f}')
            print(f'EDS cutoff negativee: {eds_cutoff_neg:.6f}')
    elif eds_all and filename == None:
            eds_df = multiple_eds(cutoff)
            eds_df.to_csv('EDS_data.csv',index=False)
            

