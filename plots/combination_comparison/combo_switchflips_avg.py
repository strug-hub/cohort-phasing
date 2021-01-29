#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description =
                                 "calculates the average switch/flip error for each chromosome and appends it to a file")
parser.add_argument("-c", "--compare", metavar="FILE", help = "Whatshap compare result in tsv format.")
parser.add_argument("-o", "--out", metavar = "FILE", default = "switch_flip_average.txt", help =
                    "Output file to write the average switch/flip error ")
args = parser.parse_args()

def average_switchflip(compare_results, summary_file):
  file = open(compare_results, "r")
  out = open(summary_file, "a+")
  switch_flip = [] #switchflip of every chromosome
  for line in file:
    if line[0] == "#":
      continue
    else:
      switch_flip.append(float(line.split("\t")[12]))
  if len(switch_flip) == 0:
    print("ERROR cannot read switchflip error rate in file: {}".format(compare_results))
  avg_switchflip = sum(switch_flip) / len(switch_flip)
  out.write("{}\t{}\n".format(compare_results, str(avg_switchflip)))
  file.close()
  out.close()
  
average_switchflip(args.compare, args.out)