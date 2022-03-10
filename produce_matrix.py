#!/usr/bin/python3

import sys
import argparse
import re
import numpy as np

my_parser = argparse.ArgumentParser()
my_parser.add_argument("-i", "--input", help="input file generated by clmdist")
my_parser.add_argument("-o", "--output", help="output matrix file generated")

args = my_parser.parse_args()

if args.input is None:
  print('Error: Input file not specified', file=sys.stderr)
  sys.exit(1)  

if args.output is None:
  print('Error: Output path not specified', file=sys.stderr)
  sys.exit(1)

output_path = args.output
file_path = args.input

d1_list = []
d2_list = []
n1_list = []
n2_list = []

with open(file_path, "r") as infile:
    for line in infile:
        m = re.search("d1=\d+\sd2=\d+", line)
        n = re.search("n1=\S+\sn2=\S+", line)
        if m and n:
            m = m.group().split("\t")
            n = n.group().split("\t")
            d1_list.append(m[0])
            d2_list.append(m[1])
            n1_list.append(n[0])
            n2_list.append(n[1])            
        else:
            print('No match')

d1_list = [re.sub("^.+=", "", i) for i in d1_list]
d2_list = [re.sub("^.+=", "", i) for i in d2_list]
n1_list = [re.sub("^.+/", "", i) for i in n1_list]
n2_list = [re.sub("^.+/", "", i) for i in n2_list]


clusters_list = []
for i in n1_list:
    if i not in clusters_list:
        clusters_list.append(i)
for i in n2_list:
    if i not in clusters_list:
        clusters_list.append(i)

mat = np.empty(shape=[len(clusters_list)+1, len(clusters_list)+1], dtype="U100")

for i in range(len(clusters_list)):
    mat[i+1][0] = clusters_list[i]
    mat[0][i+1] = clusters_list[i]

for i in range(len(d1_list)):
    n1 = n1_list[i]
    n2 = n2_list[i]
    d1 = d1_list[i]
    d2 = d2_list[i]
    n1_index = clusters_list.index(n1) + 1
    n2_index = clusters_list.index(n2) + 1
    mat[n1_index][n2_index] = d1
    mat[n2_index][n1_index] = d2

for i in range(len(clusters_list)+1):
    mat[i][i] = 0


outfile = output_path + "mat.txt"
np.savetxt(outfile, mat, fmt='%s', delimiter='\t')