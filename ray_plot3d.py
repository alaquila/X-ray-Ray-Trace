import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

parser = argparse.ArgumentParser(description='Simple python script for scatter plotting rays.')
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Input filename. Default: stdin.')
args = parser.parse_args()

df =  pd.read_csv(args.infile, parse_dates=True, index_col=0, header=0, sep=', ')
threedee = plt.figure().gca(projection='3d')
threedee.scatter(df['position x'], df['position y'], df['position z'])
threedee.set_xlabel('position x')
threedee.set_ylabel('position y')
threedee.set_zlabel('position z')
plt.show()
