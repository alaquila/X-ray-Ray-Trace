import sys
import argparse
import matplotlib.pyplot as plt

#added by PS 07/17/2023
parser = argparse.ArgumentParser(description='Simple python script for plotting detector data.')
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Input filename. Default: stdin.')
parser.add_argument('-c', default='gray', help='Colormap for plot')
args = parser.parse_args()


pixel_array = []
for line in args.infile: # read each line
    pixel_array.append([float(x) for x in line.split(',')])

plt.imshow(pixel_array,cmap=args.c)
plt.colorbar()
               
plt.show()

