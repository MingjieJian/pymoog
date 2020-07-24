#!/usr/bin/python
# Internal functions for renewing the database of stellar atmosphere model and linlist.
# WARNING: the dependene in this module may not be completly satisified, and functions may can only run on Mingjie's computer.

import numpy as np
import pandas as pd
import os 
from pymoog import model
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import line_data

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_FEB2017/'.format(os.environ['HOME'])
MOOG_run_path = '{}/.pymoog/rundir/'.format(os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(os.environ['HOME'])

def value2pm(value):
    '''
    Transform the metallicity value to Kurucz format.
    Example: -1.0 -> m10
    
    Parameters
    ----------
    value : float
        The value of metallicity.
    '''
    if value < 0:
        return 'm{:02.0f}'.format(np.abs(value)*10)
    else:
        return 'p{:02.0f}'.format(np.abs(value)*10)
    
def split_kurucz_model():
    '''
    Split the Kurucz model into single. Internal function.
    ''' 
    grid_kurucz = pd.read_csv('files/grid_points_kurucz.csv')
    for m_h in grid_kurucz.groupby('m_h').size().index:
        file = open('files/model/kurucz/standard/a{}k2.dat'.format(value2pm(m_h)))
        content = file.readlines()
        is_first = True
        for line in content:
            if 'EFF ' in line:
                if not(is_first):
                    with open('files/model/kurucz/standard/single/teff{:.0f}logg{:.1f}m_h{:+.1f}.dat'.format(teff, logg, m_h), 'w') as w_file:
                        w_file.writelines(model_line)
                teff = float(line[5:13])
                logg = float(line[21:29])  
                model_line = [line]
                is_first = False
            else:
                model_line.append(line)
        with open('files/model/kurucz/standard/single/teff{:.0f}logg{:.1f}m_h{:+.1f}.dat'.format(teff, logg, m_h), 'w') as w_file:
            w_file.writelines(model_line)

def search_grid_point_kurucz():
    '''
    The function to search all the grid points of Kurucz model and save the list to grid_path.
    The search is limit to standard model with microturbulent = 2.
    Internal use
    '''
    teff_range = np.arange(3500, 50001, 250)
    logg_range = np.arange(0, 5.1, 0.5)
    m_h_range = np.concatenate([np.arange(-5, -0.4, 0.5), np.arange(-0.3, 0, 0.1), [0], np.arange(0.1, 0.35, 0.1) ,[0.5, 1]])

    grid_point_kurucz = []
    for m_h in m_h_range:
        for teff in teff_range:
            for logg in logg_range:
                if os.path.isfile('files/model/kurucz/standard/single/teff{:.0f}logg{:.1f}m_h{:+.1f}.dat'.format(teff, logg, m_h)):
                    _, b, _ = model.read_Kurucz_model('files/model/kurucz/standard/single/teff{:.0f}logg{:.1f}m_h{:+.1f}.dat'.format(teff, logg, m_h))
                    length = b.shape[0]
                    column = b.shape[1]
                    if len(grid_point_kurucz) == 0:
                        grid_point_kurucz = np.array([[teff, logg, m_h, length, column]])
                    else:
                        grid_point_kurucz = np.concatenate([grid_point_kurucz, np.array([[teff, logg, m_h, length, column]])])
    grid_kurucz = pd.DataFrame(grid_point_kurucz, columns=['Teff', 'logg', 'm_h', 'length', 'column'])
    return grid_kurucz

def plot_model_grid():
    '''
    Plot the grid of models in each metallicity.
    Internal use.
    '''
    
    grid_kurucz = pd.read_csv('files/grid_points_kurucz.csv')
    for m_h in grid_kurucz.groupby('m_h').size().index:
        plt.figure(figsize=(13,4))
        index = grid_kurucz['m_h'] == m_h
        grid_matrix = np.array(grid_kurucz.loc[index, ['Teff', 'logg']])
        tri = Delaunay(grid_matrix)
        for i in range(len(tri.simplices)-1, -1, -1):
            if min(grid_matrix[tri.simplices[i]][:,0]) >= 35000:
                teff_gap = 5000
            else:
                teff_gap = 1500
            if  np.ptp(grid_matrix[tri.simplices[i]][:,0]) >= teff_gap or np.ptp(grid_matrix[tri.simplices[i]][:,1]) > 0.5:
                tri.simplices = np.concatenate([tri.simplices[:i], tri.simplices[i+1:]])

        plt.triplot(grid_matrix[:,0], grid_matrix[:,1], tri.simplices, zorder=0, lw=1, color='gray',alpha=0.5)
        if m_h < 0.5:
            plt.plot([50000, 42500], [5, 5], color='gray', zorder=0, alpha=0.5, lw=1)
        elif m_h == 0.5:
            plt.plot([45000, 40000], [5, 5], color='gray', zorder=0, alpha=0.5, lw=1)
        elif m_h == 1:
            plt.plot([40000, 37500], [5, 5], color='gray', zorder=0, alpha=0.5, lw=1)

        
        plt.scatter(grid_kurucz.loc[index & (grid_kurucz['length']==72), 'Teff'], grid_kurucz.loc[index & (grid_kurucz['length']==72), 'logg'], s=5, label='Model length: 72')
        plt.scatter(grid_kurucz.loc[index & (grid_kurucz['length']==64), 'Teff'], grid_kurucz.loc[index & (grid_kurucz['length']==64), 'logg'], s=5, c='C3', label='Model length: 64')

        plt.legend()
        plt.xlim((1175, 52325))
        plt.title('[Fe/H] = {:.1f}'.format(m_h))
        plt.xlabel(r'$T_\mathrm{eff}$'); plt.ylabel('logg')
        plt.gca().invert_xaxis(); plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig('../docs/img/grid_points_kurucz/m_h{:+.1f}.png'.format(m_h), dpi=250)
        plt.close()
        
def combine_linelist()
    for ele in ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Fe']:
        if ele == 'H':
            vald = line_data.read_linelist('files/linelist/vald_H')
        else:
            vald = pd.concat([vald, line_data.read_linelist('files/linelist/vald_{}'.format(ele))])

    vald.sort_values('wavelength', inplace=True)
    vald.reset_index(drop=True, inplace=True)

    line_data.save_linelist(vald, 'files/linelist/vald')