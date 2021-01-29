#!/usr/bin/env python3

import argparse
import gzip

parser = argparse.ArgumentParser(description =
                                 "Compares adjacent heterozygous variants, to determine if they are phased correctly, phased incorrectly, " +
                                 "or unphased vs the distance between variants. Bins are inclusive of the start position and exlusive of the end position ([start, end))")
parser.add_argument("-t", "--truth", metavar="VCF", help = "truth phased variants")
parser.add_argument("-q", "--query", metavar = "VCF", help = "query phased variants")
parser.add_argument("-e", "--error_rate", metavar = "TSV", help = "error rate of each bin, incorrectly_phase / total_phased",
                    default = "bins_error_rate.tsv")
parser.add_argument("-s", "--bin_stats", metavar = "TSV", help = "statistics of each bin used to make histogram of the distance and phase accuracy distribution",
                    default = "bins_stats.tsv")
parser.add_argument("-b", "--bin_size", metavar = "INT", help = "size of each bins in number of nucleotides, last bin contains everything to infinity [default = 2000]",
                    default = 2000)
parser.add_argument("-n", "--num_bins", metavar = "INT", help = "total number of bins, last bin size goes to infinity [default = 16]",
                    default = 16)
args = parser.parse_args()


#skip header line
#if write_header == True; copy the header to the output file
def skip_header_lines(vcf_file):
  if vcf_file[-3:] == ".gz" or vcf_file[-2:] == ".GZ":
    vcf =  gzip.open(vcf_file, 'rt')
  else:
    vcf = open(vcf_file, "r")
    
  for line in vcf:
    if line[:2] == "##":
      continue
    else:
      return vcf
    
def cis_trans_phase(last_gt, current_gt):
  if last_gt == current_gt:
    return "cis"
  elif last_gt == current_gt[::-1]:
    return "trans"
  else:
    print("problem with cis_trans_phase, last gt: {}; current gt: {}".format(last_gt, current_gt))
    
def create_bins(bin_size, num_bins):
  bins = [] # [[start, end, phase_correct, phased_incorrect, unphased] ... ]
  ind = 0
  last_val = 0
  while ind < num_bins:
    if ind == num_bins - 1:
      bins.append([last_val, float("inf"), 0, 0, 0])
    else:
      bins.append([last_val, last_val + bin_size ,0,0,0])
      last_val += bin_size
    ind += 1
  return bins

'''
#VCF columns
var = sv.split("\t")
chromo = var[0]
pos = int(var[1])
id_ = var[2]
ref = var[3]
alt = var[4]
mapq = var[5]
filt = var [6]
info = var[7]
form = var[8]
sample = var[9]
'''


def dist_phase(truth_vcf_file, query_vcf_file, hist_file, error_rate_file, bin_size, num_bins):
  hist = open(hist_file, "w+")
  hist.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("bin_start", "bins", "bin_size","phase_status", "count", "percent", "total_var"))
  error_rate = open(error_rate_file, "w+")
  error_rate.write("bin_start\tbins\terror_rate\ttotal_phased\ttotal_vars\n")

  count = 0
  last_chromo = 0
  last_query_pos = 0
  last_query_ps = ""
  last_var_switch = False #if True, last variant was a switch, 2 switch in a row = flip error
  
  bins = create_bins(bin_size, num_bins)
    
  truth_vcf = skip_header_lines(truth_vcf_file)
  query_vcf = skip_header_lines(query_vcf_file)
  
  # ps_blocks = {chr: {ps : [start_pos, [switch_errors, ...], [flip_errors, ...]], ...} ...}
  ps_blocks = {}
  
  for truth_line, query_line in zip(truth_vcf, query_vcf):
    truth_var = truth_line.strip("\n").split("\t")
    query_var = query_line.strip("\n").split("\t")
    
    #checking the positions match
    assert truth_var[1] == query_var[1], "VCF variants out of order. {} \n\n{}".format(truth_var, query_var)
    chromo = query_var[0]
    pos = int(truth_var[1])
    query_sample = query_var[9].split(":")
    truth_sample = truth_var[9].split(":")
    query_form = query_var[8].split(":")
    truth_form = truth_var[8].split(":")
    query_gt = query_sample[query_form.index("GT")]
    truth_gt = truth_sample[truth_form.index("GT")]
    
    #skip multiallelic because WhatsHap cannot phase them, for a fair comparison 
    if "2" in query_gt:
      continue
    
    try:
       query_ps = query_sample[query_form.index("PS")]
    except:
       query_ps = "unphased"

    if query_gt[0] != query_gt[2] and truth_gt[0] != truth_gt[2]: #both truth and query called heterozygous genotype
    
      if len(query_sample[0]) == 1:
        continue
      
      if chromo == last_chromo:
        count += 1
        dist = int(query_var[1]) - last_query_pos
        if query_sample[0][1] == "|" and query_ps == last_query_ps: #phased and same ps
          phase_stats = "phased"
          if last_query_gt[1] == "/": #unphased with PS, longranger does this sometimes
            phase_stats = "unphased"
            
        elif query_sample[0][1] == "/" or query_ps != last_query_ps: #unphased or diff ps
          phase_stats = "unphased"
           
        
        ind_bins = dist // bin_size
        #goes into the last bin
        if ind_bins >= num_bins:
          ind_bins = num_bins - 1 #correct for 0-based indexing of ind
          
        if phase_stats == "phased":
          #phasing error
          truth_phase = cis_trans_phase(last_truth_gt, truth_gt)
          query_phase = cis_trans_phase(last_query_gt, query_gt)
          
          #correct phase
          if truth_phase == query_phase:
            bins[ind_bins][2] += 1
            last_var_switch = False

          #incorrect phase
          else:
            bins[ind_bins][3] += 1
            
            if chromo not in ps_blocks.keys():
              ps_blocks[chromo] = {query_ps : [pos, [],[]]}
              last_var_switch = False #first incorrectly phased var in chromosome or phase block
              last_error_pos = pos
            elif query_ps not in ps_blocks[chromo].keys():
              ps_blocks[chromo][query_ps] = [pos, [], []]
              last_var_switch = False
              last_error_pos = pos
              
            dist = pos - ps_blocks[chromo][query_ps][0]
            if last_var_switch == True:
              
               if last_error_pos in ps_blocks[chromo][query_ps][1]:
                  ps_blocks[chromo][query_ps][1].pop()
                  
               ps_blocks[chromo][query_ps][2].append(last_error_pos)
               ps_blocks[chromo][query_ps][2].append(dist)
            else:
              ps_blocks[chromo][query_ps][1].append(dist)
              
            last_var_switch = True
            


        elif phase_stats == "unphased":
          bins[ind_bins][4] += 1

              
        last_query_pos = int(query_var[1])
        last_query_ps = query_ps
        last_query_gt = query_gt
        last_truth_gt = truth_gt
        if phase_stats == "phased":
          last_phased_sample = query_sample
          last_phased_form = query_form
          last_phased_pos = query_var[1]

            
      elif chromo !=  last_chromo:
        last_chromo = chromo
        last_query_pos = int(query_var[1]) 
        last_query_ps = query_ps
        last_truth_gt = truth_gt
        last_query_gt = query_gt
        last_phased_sample = query_sample
        last_phased_form = query_form
        last_phased_pos = query_var[1]

  for b in bins:
    # bins format: [[start,end,phase_correct, unphased, phased_incorrect] ... ]
    total = b[2] + b[3] + b[4]
    if total == 0 :
      print("0 variants in bin {} - {}".format(b[0],b[1]))
    p_phased = (b[2] + b[3]) / total
    p_unphased = b[4] / total
    if b[2] == 0:
      p_incorrect = 0
      p_correct = 0
      e_rate = 'NA'
    else:
      p_incorrect = b[3]/total
      p_correct = b[2] / total
      e_rate = b[3]/ p_phased
    if b[2] < b[3]:
      print("error phased_incorrect > phased_correct".format(b))

    
    error_rate.write("{}\t[{},{})\t{}\t{}\t{}\n".format(b[0], b[0], b[1], e_rate, b[2], total))
  
    #hist: "bin_start", "bins", "bin_size","phase_status", "count", "percent", "total_var"
    hist.write("{}\t[{},{})\t{}\t{}\t{}\t{}\t{}\n".format(b[0], b[0], b[1], b[1] - b[0], "phased_correct", b[2], p_correct, total))
    hist.write("{}\t[{},{})\t{}\t{}\t{}\t{}\t{}\n".format(b[0], b[0], b[1], b[1] - b[0], "phased_incorrect", b[3], p_incorrect, total))
    hist.write("{}\t[{},{})\t{}\t{}\t{}\t{}\t{}\n".format(b[0], b[0], b[1], b[1] - b[0], "unphased", b[4], p_unphased, total))
  
  print("total number of variants: {}".format(count))
  query_vcf.close()
  truth_vcf.close()
  hist.close()
  error_rate.close()

dist_phase(args.truth, args.query, args.bin_stats, args.error_rate, int(args.bin_size), int(args.num_bins))
