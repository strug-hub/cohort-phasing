#!/usr/bin/env python3


import argparse

parser = argparse.ArgumentParser(description =
                                 "calculates the average switch/flip error for each chromosome and appends it to a file")
parser.add_argument("-c", "--compare", metavar="FILE", help = "Whatshap compare result in tsv format.")
parser.add_argument("-o", "--out", metavar = "FILE", default = None, help =
                    "")
parser.add_argument("-p", "--pos", metavar = "FILE", default = None, help =
                    "coordinates of plot")
args = parser.parse_args()

if args.out == None:
  out_file = args.compare + "_parsed.txt"
else:
  out_file = args.out
  
if args.pos == None:
  pos_file = args.compare + "_pos.txt"
else:
  pos_file = args.pos
  

#input:
#compare_file - two column tab separated file, first column: file name, second column: average switch/flip error over all chromosomes
def count_datasets_compare(compare_file, out_file, compare_pos_file):
  compare =  open(compare_file, "r")
  out = open(out_file, "w+")
  out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("file_name", "switch_flip_rate", "ONT", "10xG", "PacBio_CCS", "PacBio_CLR", "num_datasets"))
  pos = open(compare_pos_file, "w+")
  pos.write("{}\t{}\t{}\n".format("x", "y", "dataset"))
  order_tech = ["ONT", "10xG", "PacBio_CCS", "PacBio_CLR"]
  order_pos = [[],[0],[-0.05, 0.05], [-0.1, 0, 0.1], [-0.15, -0.05, 0.05, 0.15]]
  col = 0 #filename col
  for line in compare:
    combo = line.strip("\n").split("\t")
    start = combo[col].index("compare_w") + 9
    end = combo[col][start:].index("_")
    code = combo[col][start:][:end] #weights
    code_equal = ""
    ind = 0
    while ind < len(code):
      if code[ind:ind+3] == ".05":
        ind += 3
      elif code[ind] == ".":
        ind += 1
      else:
        code_equal += code[ind]
      ind += 1
    num_dataset = sum([int(n) for n in code_equal])
    out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(line.strip("\n"), code_equal[0], code_equal[1], code_equal[2], code_equal[3], num_dataset))
    count  = 0
    ind = 0
    while count < num_dataset:
      if code_equal[ind] == "1":
        pos.write("{}\t{}\t{}\n".format(num_dataset + order_pos[num_dataset][count], combo[1], order_tech[ind]))
        count += 1
      ind += 1
  pos.close()
  out.close()
  compare.close()

count_datasets_compare(args.compare, out_file, pos_file)