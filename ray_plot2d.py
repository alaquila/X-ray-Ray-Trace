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
	print 'Error you need to select two and only two axis to plot.'
	sys.exit()

df =  pd.read_csv(args.infile, parse_dates=True, index_col=0, header=0, sep=', ')

if args.X and args.Y:
	df.plot(x='position x', y='position y', kind='scatter')
elif args.X and args.Z:
	df.plot(x='position x', y='position z', kind='scatter')
else:
	df.plot(x='position y', y='position z', kind='scatter')

plt.show()
