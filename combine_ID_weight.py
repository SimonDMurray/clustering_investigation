#! /usr/bin/env python

import csv
import argparse
import sys

my_parser = argparse.ArgumentParser()

my_parser.add_argument("-d", "--dist_file", help="the input distance file")
my_parser.add_argument("-c", "--cells_file", help="the input cells file")
my_parser.add_argument("-o", "--output", help="path to output combined file to")

args = my_parser.parse_args()

if args.dist_file is None:
  print('Error: Dist file is not specified')
  sys.exit()

if args.cells_file is None:
  print('Error: Cells file is not specified')
  sys.exit()

if args.output is None:
  outpath="./"
else:
  outpath = args.output

dist_path = args.dist_file
cells_path = args.cells_file

dist_dict = {}

with open(dist_path, "r") as dist_file:
  for dist in dist_file:
    dist = dist.strip("\n")
    dist_list = dist.split()
    dist_dict[dist_list[0].strip('"')] = dist_list[1::]

invert_dict = {}

for key, dist_list in dist_dict.items():
  invert_list = []
  for dist_val in dist_list:
    if dist_val == "0":
      invert_list.append((dist_val))
    else:
      invert_list.append(str(1/float(dist_val)))
  invert_dict[key] = invert_list

cell_dict = {}

with open(cells_path, "r") as cell_file:
  for cell in cell_file:
    cell = cell.strip("\n")
    cell_list = cell.split()
    cell_dict[cell_list[0]] = cell_list[1::]

common_pairs = []

for key in dist_dict.keys():
  if (key in cell_dict.keys()):
    common_pairs.append(key)

assert len(common_pairs) == len(dist_dict.keys())
assert len(common_pairs) == len(cell_dict.keys())

for key, value in dist_dict.items():
  assert len(value) == len(cell_dict[key])

conc_dict = {}
conc_list = []

for key in cell_dict.keys():
  for val in cell_dict[key]:
    index = cell_dict[key].index(val)
    conc_val = val + ":" + invert_dict[key][index]
    conc_list.append(conc_val)
  conc_dict[key] = conc_list
  conc_list = []

output=outpath + "nn_combined.txt"
with open(output, 'w') as conc_file:
  for key in conc_dict.keys():
    out_list = [key]
    for val in conc_dict[key]:
      out_list.append(val)
    tsv_output = csv.writer(conc_file, delimiter='\t')
    tsv_output.writerow(out_list)
