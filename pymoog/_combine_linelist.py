import line_data
import pandas as pd

for ele in ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Fe']:
    if ele == 'H':
        vald = line_data.read_linelist('files/linelist/vald_H')
    else:
        vald = pd.concat([vald, line_data.read_linelist('files/linelist/vald_{}'.format(ele))])

vald.sort_values('wavelength', inplace=True)
vald.reset_index(drop=True, inplace=True)

line_data.save_linelist(vald, 'files/linelist/vald')
