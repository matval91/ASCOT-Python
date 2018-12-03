import numpy as np

def write_wall_2D(W_fname):
    """
    Function to write the wall file for ASCOT from an ascii file in the format (R,Z)
    """
    wall_file= open(W_fname)
    lines = wall_file.readlines()
    R_W = np.zeros(len(lines), dtype = float)
    Z_W = np.zeros(len(lines), dtype = float)
    div_flag = np.zeros(len(lines), dtype = int)

    for i,line in enumerate(lines):
        tmp = line.split()
        R_W[i] = float(tmp[0])
        Z_W[i] = float(tmp[1])


    out_fname = "input.wall_2d"
    outfile = open(out_fname, 'wa')
    outfile.write('{:d} (R,z) wall points & divertor flag (1 = divertor, 0 = wall)\n'.format(len(lines)))
    np.savetxt(outfile, np.c_[R_W, Z_W, div_flag], fmt='%5.4f %5.4f %1d')


