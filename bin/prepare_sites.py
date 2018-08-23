import os
import sys
import shutil
import subprocess
import argparse
import pandas as pd

from MOBPred import moe
from mdtools.utility import reader
from mdtools.utility import utils

boxsize = 30.0

parser = argparse.ArgumentParser(description="Prepare sites")

parser.add_argument('-r',
    type=str,
    dest='input_files_r',
    nargs='+',
    required=True,
    help = 'Target files: .pdb')

parser.add_argument('-l',
    type=str,
    dest='input_files_l',
    nargs='+',
    help='PDB files containing the ligand, used to predict the binding site')

parser.add_argument('-sitefinder',
    action='store_true',
    default=False,
    help='Use moe sitefinder to find the binding sites')

parser.add_argument('-csv',
    type=str,
    dest='csvfile',
    default='sites.csv',
    help='CSV file name')

parser.add_argument('-minplb',
    dest='minplb',
    type=float,
    default=1.0,
    help="Minimimum PLB value to consider")

parser.add_argument('-nsitesmax',
    dest='nsitesmax',
    type=int,
    default=0,
    help="Maximum number of sites kept")

parser.add_argument('-minsize',
    type=float,
    dest='minsize',
    default = 30.0,
    help = 'Min size for the box')

args = parser.parse_args()

for file_r in args.input_files_r:
    if not os.path.exists(file_r):
        raise ValueError("File %s not found!"%(file_r))

if args.input_files_l:
    if args.sitefinder:
        raise ValueError('Sitefinder option selected and list of ligand files provided!')
    for file_l in args.input_files_l:
        if not os.path.exists(file_l):
            raise ValueError("File %s not found!"%(file_l))

    nfiles_r = len(args.input_files_r)
    if len(args.input_files_l) != nfiles_r:
        raise ValueError('Number of sites (%s) should match number of target files (%s)'%(len(args.input_files_l),nfiles_r))
else:
    if not args.sitefinder:
        raise ValueError('Sitefinder option not selected and list of ligand files empty!')

curdir = os.getcwd()

info = {}
features = ['target', 'center', 'size']
for ft in features:
    info[ft] = []

pwd = os.getcwd()

for idx, file_r in enumerate(args.input_files_r):
    label = 'target' + (3-len(str(idx+1)))*'0' + str(idx+1)

    file_r_base = os.path.basename(file_r)
    if args.sitefinder:
         dir_r = os.path.dirname(file_r)
         if not dir_r:
             dir_r = '.'
         sitedir_r = dir_r+'/sitefinder'
         shutil.rmtree(sitedir_r, ignore_errors=True)
         os.mkdir(sitedir_r)
         os.chdir(sitedir_r)

         # use moe sitefinder
         moe.write_sitefinder_script('run_site_finder.sh', '../'+file_r_base, args)
         subprocess.check_call('bash run_site_finder.sh', shell=True)
         with open('moebatch.log', 'r') as ff:
            ff.next()
            nsites = 0
            for idx, line in enumerate(ff):
                line_s = line.split()
                plb = float(line_s[1])

                if (plb > args.minplb and idx+1 <= args.nsitesmax) or idx==0:
                    info['target'].append(label)
                    line_s = line.split()
                    if 'plb' not in info:
                        info['plb'] = []
                    info['plb'].append(plb)

                    # get box center
                    center = ', '.join(['%.2f'%float(x) for x in line_s[2:5]])
                    info['center'].append(center)

                    # get boxsize
                    radius = float(line_s[5])
                    boxsize = max(2*radius,args.minsize)
                    info['size'].append(', '.join(['%.2f'%boxsize for jdx in range(3)]))
                    if 'site' not in info:
                        info['site'] = []
                    info['site'].append('%s'%(idx+1))
         os.chdir(pwd)         
    else:
        file_l = args.input_files_l[idx]
        rpdb = reader.open(file_l)
        coords = []

        base, ext = os.path.splitext(file_l)
        for line in rpdb.next()['ATOM']:
            if ext == '.pdb':
                coords.append(map(float,[line[4], line[5], line[6]]))
            elif ext == '.mol2':
                coords.append(map(float,[line[2], line[3], line[4]]))
            else:
                raise ValueError('Format not recognized!')

        cog = utils.center_of_geometry(coords)
        center = ', '.join(['%.2f'%float(x) for x in cog])
        info['target'].append(label)
        info['center'].append(center)
        info['size'].append(', '.join(['%.2f'%boxsize for jdx in range(3)]))

info = pd.DataFrame(info)
info = info[features].to_csv(args.csvfile, index=False)
