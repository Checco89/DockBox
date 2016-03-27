import sys
import os
import shutil
import numpy as np
import fileinput
import subprocess
import glob

import util.pdbtools as pdbt

default_settings = {'program': 'cpptraj', 'forcefield': 'leaprc.ff14SB', 'cutoff': '2.0'}

def find_consensus(config):

    curdir = os.getcwd()
    workdir = 'consensus'

    if os.path.isdir(workdir):
        shutil.rmtree(workdir)
    os.mkdir(workdir)
    os.chdir(workdir)

    nposes = prepare_files_for_consensus(config)
    if config.consensus['program'] == 'cpptraj':
        run_tleap(config, nposes)
        run_cpptraj(config, nposes)
        extract_results(config)

    os.chdir(curdir)

def prepare_files_for_consensus(config):

    curdir = os.getcwd()
    pdbdir = 'PDB'
    os.mkdir(pdbdir)

    # write files containing the number of poses
    # generated by each software
    ff = open(pdbdir+'/info.dat', 'w')
    print >> ff, '#program        nposes'
      
    sh = 0 # shift of model idx
    for program in config.programs:
        prgdir = '../%s'%program
        os.chdir(prgdir)

        # prepare ligand structure
        shutil.copyfile('lig.out.pdb', 'lig-c.out.pdb')
        pdbt.format_lig_file('lig-c.out.pdb')
        pdbt.give_unique_atom_names('lig-c.out.pdb') # give unique atom names

        # remove hydrogens to leave tleap add them
        pdbt.remove_hydrogens('rec.out.pdb', 'rec-c.out.pdb')

        # create the corresponding file both the receptors and ligands
        pdbt.create_reclig_file('rec-c.out.pdb', 'lig-c.out.pdb', 'rec-lig.out.pdb')
        os.chdir(curdir)
        with open(prgdir+'/rec-lig.out.pdb') as pdbfi:
            for line in pdbfi:
                if line.startswith('MODEL'):
                    try:
                        pdbfo.close()
                    except NameError:
                        pass
                    idx = int(line.split()[1])
                    pdbfo = open(pdbdir+'/rec-lig.%s.pdb'%(idx+sh), 'w')
                elif not line.startswith('END'):
                    pdbfo.write(line)
            sh += idx
            print >> ff, '%10s        %s'%(program,idx)
    pdbfo.close()
    nposes = sh
    
    return nposes

def run_tleap(config, nposes):

    curdir = os.getcwd()

    # create antechamber dir
    antchmbdir = 'antchmb'
    os.mkdir(antchmbdir)
    os.chdir(antchmbdir)

    # get acceptable structure of the ligand
    prgdir = '../../%s'%config.programs[0]
    with open('lig.pdb', 'w') as pdbfout:
        with open(prgdir+'/lig-c.out.pdb', 'r') as pdbfin:
            for line in pdbfin:
                pdbfout.write(line)
                if line.startswith('ENDMDL'):
                    break
    subprocess.call('antechamber -i lig.pdb -fi pdb -o lig.mol2 -fo mol2 -at gaff -du y -pf y > antchmb.log', shell=True)
    subprocess.check_call('parmchk -i lig.mol2 -f mol2 -o lig.frcmod', shell=True)
    os.chdir(curdir)

    # create antechamber dir
    leapdir = 'LEaP'
    os.mkdir(leapdir)
    os.chdir(leapdir)

    # prepare tleap input file
    forcefield = config.consensus['forcefield']

    linespdb = ""
    for idx in range(nposes):
        if idx == 1:
            linespdb += """p = loadPdb ../PDB/rec-lig.%s.pdb
saveAmberParm p rec-lig.prmtop rec-lig.inpcrd
savepdb p rec-lig.%s.pdb\n"""%(idx+1,idx+1)
        else:
            linespdb += """p = loadPdb ../PDB/rec-lig.%s.pdb
savepdb p rec-lig.%s.pdb\n"""%(idx+1,idx+1)

    linespdb = linespdb[:-1]

    with open('leap.in', 'w') as file:
        script ="""source %(forcefield)s
source leaprc.gaff
LIG = loadmol2 ../%(antchmbdir)s/lig.mol2
loadamberparams ../%(antchmbdir)s/lig.frcmod
%(linespdb)s
quit"""% locals()
        file.write(script)

    # run tleap
    subprocess.check_call('tleap -f leap.in > leap.log', shell=True, executable='/bin/bash')
    os.chdir(curdir)

def run_cpptraj(config, nposes):

    curdir = os.getcwd()

    # create antechamber dir
    clusterdir = 'clstr'
    os.mkdir(clusterdir)
    os.chdir(clusterdir)

    lines_trajin = ""
    for idx in range(nposes):
        lines_trajin += "trajin ../LEaP/rec-lig.%s.pdb\n"%(idx+1)

    # remove last \n
    lines_trajin = lines_trajin[:-1]
    cutoff = config.consensus['cutoff']

    # write cpptraj config file to cluster frames
    with open('cpptraj.in', 'w') as file:
        script ="""parm ../LEaP/rec-lig.prmtop
%(lines_trajin)s
rms first "@CA,C,N & !:LIG"
cluster ":LIG & !@/H" nofit mass epsilon %(cutoff)s summary summary.dat info info.dat
"""% locals()
        file.write(script)

    subprocess.check_call('cpptraj -i cpptraj.in > cpptraj.log', shell=True)
    os.chdir(curdir)

def extract_results(config):

    poses = []
    prgms = []
    idxprgm = -1
    with open('PDB/info.dat') as fi:
        for line in fi:
            if not line.startswith('#'):
                idxprgm += 1
                prgm = line.split()[0]
                prgms.append(prgm)
                n = int(line.split()[1])
                poses.extend([idxprgm for idx in range(n)])
    nprgms = len(prgms)
    
    clstrs = []
    heterg = []
    ff = open('clstr/info.dat')
    for line in ff:
        # if line does not start with #
        if not line.startswith('#'):
            # indices = numbers of the poses involved in the current cluster
            indices = [i for i, x in enumerate(line.strip()) if x == 'X']
            clstr = []
            for idx in indices:
                clstr.append([poses[idx], idx])
            clstrs.append(clstr)
            idxs = [x[0] for x in clstr]
            # compute heterogeneity factor for the current cluster
            heterg.append(len(list(set(idxs)))*100./nprgms)
    ff.close()
 
    hetergidxs = np.argsort(-1*np.array(heterg))
    hetmax = heterg[hetergidxs[0]]
    hetmaxs = []
    # in case of equal heterogeneity factors, take the cluster with the most of poses involved
    for idx in hetergidxs:
        if heterg[idx] == hetmax:
            hetmaxs.append(len(clstrs[idx]))
    clstridx = hetergidxs[np.argmax(hetmaxs)]
    file = 'rec-lig.' + str(clstrs[clstridx][0][1]+1) + '.pdb'
    prgmfound =  [prgms[idx] for idx in list(set([x[0] for x in clstrs[clstridx]]))]
    
    with open('info.dat','w') as rf:
        print >> rf, "heterg = %.4f, pose found with %s"%(hetmax,', '.join(prgmfound))
    
    for pose in clstrs[clstridx]:
        file = 'rec-lig.' + str(pose[1]+1) + '.pdb'
        shutil.copyfile('PDB/'+file, file)
