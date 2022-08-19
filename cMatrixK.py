
#!/usr/bin/env python3
# Author : Diana Aguilar modified from Dr. Rocha
'''
cmatrixK.py
'''
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("c_matrix", help="cmatrix")
parser.add_argument("k", type=int, help="component for selection")
parser.add_argument("h", type=int, help="scaling factor")
args = parser.parse_args()

in_mat = np.loadtxt(args.c_matrix, skiprows=1)
out_mat = in_mat

if args.k == 0:
    out_mat += args.h
else:
    out_mat[args.k-1,args.k-1] += args.h
print(out_mat)

np.savetxt(args.c_matrix+"selK"+str(args.k)+"_h"+str(args.h), out_mat, header='{} {}'.format(in_mat.shape[0], in_mat.shape[1]), comme
nts='')
