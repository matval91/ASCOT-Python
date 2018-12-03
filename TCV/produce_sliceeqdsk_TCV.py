"""
Given a set of EQDSK from a TCV shot, produces the inputs

input: shot, timeslices
output: set of input.magn_header(bkg) for BBNBI/ASCOT
        set of images of equilibria

"""

import numpy as np
import ascot_Bfield
import sys
import matplotlib.pyplot as plt
shot = str(sys.argv[1])
t = np.empty(len(sys.argv)-2)
eq = np.empty(len(sys.argv)-2, dtype=object)
for i,arg in enumerate(sys.argv):
    if i<2:
        continue
    else:
        arg = float(arg)
        t[i-2] = '{:.2f}'.format(arg)

fname_ini = 'EQDSK_'+shot+'t'
fname_end = '00_COCOS17'

for ind, i in enumerate(t):
    nn = '{:.2f}'.format(i)
    fname=fname_ini+nn+fname_end
    eq[ind]=ascot_Bfield.Bfield_eqdsk(fname, 259, 259, 'TCV', 17)
    eq[ind].eqdsk_checkplot(); 
    eq[ind].build_SN(); eq[ind].write()
    
    f=plt.gcf(); f.savefig(fname_ini+nn+'.png', dpi=400)
