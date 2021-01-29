#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description =
                                 "calculates the average switch/flip error for each chromosome and appends it to a file")
parser.add_argument("-s", "--stats", metavar="FILE", help = 
                    "compiled stats file from WhatsHap of the complete summary, contains the last line of each stats file (ALL chromosome stats)")
parser.add_argument("-o", "--out", metavar = "FILE", default = None, help =
                    "copies stats_file and adds additional columns to the stats_file, 0/1 for each technology present, and total num datasets")
parser.add_argument("-b", "--phase_block", metavar="FILE", default = None, help = 
                    "")
parser.add_argument("-p", "--phased", metavar="FILE", default = None, help = 
                    "")
parser.add_argument("-n", "--n50", metavar="FILE", default = None, help = 
                    "")
args = parser.parse_args()

if args.out == None:
  out_file = args.stats + "_numdataset.tsv"
else:
  out_file = args.out  

if args.phase_block == None:
  phase_block_file = args.stats + "_phase_blocks.txt"
else:
  phase_block_file = args.phase_block
  
if args.phased == None:
  phased_file = args.stats + "_phased.txt"
else:
  phased_file = args.phased
  
if args.n50 == None:
  n50_file = args.stats + "_n50.txt"
else:
  n50_file = args.n50

#parses the stats data with unequal weights 
# stats_file - compiled stats file from WhatsHap of the complete summary, contains the last line of each stats file (ALL chromosome stats)
#outputs:
# out_file - copies stats_file and adds additional columns to the stats_file, 0/1 for each technology present, and total num datasets
# pos_file_phase_block - coordinates to plot the phase block stats
# pos_file_perc_phased_var - coordinates to plot the % phased variant for each data set
# pos_file_n50
def parse_stats(stats_file, out_file, pos_file_phase_block, pos_file_perc_phased_var, pos_file_n50):
  stats =  open(stats_file, "r")
  out = open(out_file, "w+")
  pos_block = open(pos_file_phase_block, "w+")
  pos_block.write("{}\t{}\t{}\t{}\n".format("x", "y", "dataset", "num_datasets"))
  pos_phased = open(pos_file_perc_phased_var, "w+")
  pos_phased.write("{}\t{}\t{}\t{}\n".format("x", "y", "dataset", "num_datasets"))
  pos_n50 = open(pos_file_n50, "w+")
  pos_n50.write("{}\t{}\t{}\t{}\n".format("x", "y", "dataset", "num_datasets"))
  col = 2 #file name column
  order_tech = ["ONT", "10xG", "PacBio_CCS", "PacBio_CLR",]
  order_pos = [[],[0],[-0.05, 0.05], [-0.1, 0, 0.1], [-0.15, -0.05, 0.05, 0.15]]
  for line in stats:
    if line[0] == "#": #first header line
      out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(line.strip("\n"), "ONT", "10xG", "PacBio_CCS", "PacBio_CLR", "num_datasets"))
      continue
    combo = line.strip('\n').split("\t")
    #parse the weights
    start = combo[col].index("PBlong_w") + 8
    end = combo[col][start:].index("_")
    code = combo[col][start:][:end] #weights
    code_equal = "" #equiv code as if it was equal weight
    ind = 0
    while ind < len(code):
      if code[ind:ind+3] == ".05":
        ind += 3
      elif code[ind] == ".":
        ind += 1
      else:
        code_equal += code[ind]
      ind += 1
    print(code_equal)
    num_dataset = sum([int(n) for n in code_equal])
    out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(line.strip("\n"), code_equal[0], code_equal[1], code_equal[2], code_equal[3], num_dataset))
    count  = 0
    ind = 0
    while count < num_dataset:
      if code_equal[ind] == "1":
        pos_block.write("{}\t{}\t{}\t{}\n".format(num_dataset + order_pos[num_dataset][count], combo[7], order_tech[ind], num_dataset))
        pos_phased.write("{}\t{}\t{}\t{}\n".format(num_dataset + order_pos[num_dataset][count], int(combo[4])/int(combo[18]), order_tech[ind], num_dataset))
        pos_n50.write("{}\t{}\t{}\t{}\n".format(num_dataset + order_pos[num_dataset][count], combo[21], order_tech[ind], num_dataset))
        count += 1
      ind += 1
  pos_n50.close()
  pos_phased.close()
  pos_block.close()
  out.close()
  stats.close()
  
  
parse_stats(args.stats, out_file, phase_block_file, phased_file, n50_file)