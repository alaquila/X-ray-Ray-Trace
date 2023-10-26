
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Simple python script for scatter plotting rays.')
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Input filename. Default: stdin.')
parser.add_argument('-c', default='gray', help='Colormap for plot')
parser.add_argument('--x_range', nargs=2, type=float, default=[-0.018, 0.018])
parser.add_argument('--z_range', nargs=2, type=float, default=[.4829, .5189])
parser.add_argument('--bins', type=int, default=256, help='Number of bins for the histogram')
args = parser.parse_args()

ray = pd.read_csv(args.infile, header=0, sep=', ', engine='python' ,encoding='latin-1')
ray.columns = ["col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11"]

n = ray["col1"].values
x = ray["col2"].values
z = ray["col4"].values
i = ray["col10"].values

# indices where i is equal to 0
ind_zero = np.where(i == 0)[0]

# Delete elements of ind_zero from x, z, and i
x_c = np.delete(x, ind_zero)
z_c = np.delete(z, ind_zero)
i_c = np.delete(i, ind_zero)

plt.subplot(1, 3, 1)
foo1 = np.histogram2d(x_c, z_c, bins = args.bins, range = [args.x_range,args.z_range])
plt.gray()
plt.imshow(foo1[0], origin='lower', aspect='equal')
plt.title('classical imaging', fontsize=14)

for j in range(len(n)-1):
        if n[j] % 2 == 0:
            if n[j+1] - n[j] == 1:
                if i[j+1] == 1 and i[j] == 1:
                    i[j+1] = 0
                    i[j] = 0

ind_zero = np.where(i == 0)[0]
n_gq = np.delete(n, ind_zero)
x_gq = np.delete(x, ind_zero)
z_gq = np.delete(z, ind_zero)
i_gq = np.delete(i, ind_zero)

for j in range(len(n)-1):
        if n[j] % 2 == 0:
            if n[j+1] - n[j] == 1:
                if i[j+1] == 1 and i[j] == 1:
                    i[j+1] = 0
                    i[j] = 0
                else:
                    i[j],i[j+1] = i[j+1],i[j]

ind_zero = np.where(i == 0)[0]               
n_q = np.delete(n, ind_zero)
x_q = np.delete(x, ind_zero)
z_q = np.delete(z, ind_zero)
i_q = np.delete(i, ind_zero)

plt.subplot(1, 3, 2)
foo2 = np.histogram2d(x_q, z_q, bins = args.bins, range = [args.x_range,args.z_range] )
plt.imshow(foo2[0], origin='lower', aspect='equal')
plt.title('direct quantum imaging', fontsize=14)


plt.subplot(1, 3, 3)
foo3 = np.histogram2d(x_gq, z_gq, bins = args.bins, range = [args.x_range,args.z_range] )
plt.imshow(foo3[0], origin='lower', aspect='equal')
plt.title('ghost quantum imaging', fontsize=14)


plt.show()
