import sys
sys.path.append("/home/aleksey/software/COSMIC/")
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.evolve import Evolve
# from cosmic.evolve import Evolve
sys.path.append("/home/aleksey/software/mse")
import mse
import numpy as np
import pickle

from bash_command import bash_command as bc

ms_grid = np.geomspace(10,150)
# ms_grid = [20.602]
# ms_grid = [10]
m2 = 1.019
zz = 0.0126
ecc = 0.721307

with open("comp_mse_cosmic_final", "w") as ff:
    ff.write("#mi mse cosmic\n") 

for jj, mm in enumerate(ms_grid):
    single_binary = InitialBinaryTable.InitialBinaries(m1=mm, m2=m2, porb=1000000,\
                                                       ecc=ecc, tphysf=10.0, kstar1=1, kstar2=1,\
                                                       metallicity=zz)

    BSEDict = {'xi': 0.5, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1,\
               'alpha1': 1.0, 'pts1': 0.05, 'pts2': 0.01, 'pts3': 0.02, 'epsnov': 0.001,\
               'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1,\
               'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 1, 'eddfac': 1, 'ifflag': 0,\
               'bconst': 3000, 'sigma': 0.0, 'gamma': -2.0, 'pisn': -2,\
               'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]],\
               'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90.0,\
               'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],\
               'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 1, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6,\
               'aic' : 1, 'ussn' : 1, 'sigmadiv' :-20.0, 'qcflag' : 5, 'eddlimflag' : 0,\
               'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0],\
               'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 0, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1}


    bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=single_binary, BSEDict=BSEDict, dtp=0.01)
    bpp.to_csv("bpp_{0}.csv".format(jj))
    bcm.to_csv("bcm_{0}.csv".format(jj))


    ms1=[mm, m2]
    metallicities=[zz, zz]
    eccentricities=[ecc]
    inclinations=[0]
    arguments_of_pericentre=[0]
    longitudes_of_ascending_node=[0]
    semimajor_axes=[1e5]
    mse.Tools.evolve_system("fully_nested", 2, ms1, metallicities, semimajor_axes, eccentricities, inclinations,\
    arguments_of_pericentre, longitudes_of_ascending_node, 1e7, 5000, save_data=True,\
    plot_filename='tmp', show_plots=False)

    bc.bash_command("tar -xvf tmp_archive_tar_pickle")
    with open('tmp.pkl', 'rb') as ff:
        dat=pickle.load(ff)
        ms=[dat['log'][ii]['particles'][0].mass for ii in range(len(dat['log']))]
        ts=[dat['log'][ii]['time'] for ii in range(len(dat['log']))]

    with open("comp_mse_cosmic_final", "a") as ff:
        ff.write("{0} {1} {2}\n".format(mm, ms[-1], np.array(bpp['mass_1'])[-1]))
        ff.flush()
    
    print(len(bcm))
