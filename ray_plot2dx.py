import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Simple python script for scatter plotting rays.')
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Input filename. Default: stdin.')
parser.add_argument('-X', action='store_true', default=False, help='Axis direction for plot')
parser.add_argument('-Y', action='store_true', default=False, help='Axis direction for plot')
parser.add_argument('-Z', action='store_true', default=False, help='Axis direction for plot')
args = parser.parse_args()

if args.X + args.Y + args.Z != 2:
        print('Error you need to select two and only two axis to plot.')
        sys.exit()

df =  pd.read_csv(args.infile, parse_dates=True, index_col=0, header=0, sep=', ')

# added by PS 07/12/2023: fix aspect ratio, plot intensity as grayscale 
if args.X and args.Y:
        ax = df.plot(x='position x', y='position y', c='intensity', kind='scatter')
        ax.set_aspect('equal')
elif args.X and args.Z:
        ax = df.plot(x='position x', y='position z', c='intensity', kind='scatter')
        ax.set_aspect('equal')
else:
        ax = df.plot(x='position y', y='position z', c='intensity', kind='scatter')
        ax.set_aspect('equal')

plt.show()


