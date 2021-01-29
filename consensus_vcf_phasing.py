#!/usr/bin/env python3

import argparse
import sys
import gzip
import os

parser = argparse.ArgumentParser(description =
                                 "Creates a consesnsus phased VCF based on inputted phased VCFs. " +
                                 "The phased VCFs must contain the same variants and genotype in the same order. " +
                                 "If there are ties in phase call, the variant is unphased. " +
                                 "Makes changes to the INFO and sample (GT and PS) column of the first inputted phased VCF file and writes to the output file. " +
                                 "Header is copied from the first inputted VCF. ")
parser.add_argument("vcf", metavar="VCF", nargs = "+", help = "Two or more single-sampled phased VCFs to create a consensus VCF. " +
                    "Can be compressed (.gz) or uncompressed.")
parser.add_argument("-n", "--names", metavar="STRINGS", help = "Names of each dataset, comma separated and in the order they are inputted. [optional]")
parser.add_argument("-o", "--out", metavar = "VCF", default = "consensus.vcf", help =
                    "Output consensus file name. File is unzipped [default = consensus.vcf]")
parser.add_argument("-w", "--weights", metavar='STRING', default = None, help =
                    "Weights of each of the phased VCF inputs, numbers appearing in the same order and commas separated.\n" +
                    "[default: all input weighted 1]")
parser.add_argument("-t", "--threshold", metavar= "NUM", default = 0, help =
                    "Minimum weight required to maintain a phase block.\n" +
                    "When the sum of weights of datasets that supports a phasing subtracted by the weighted votes of datasets that spports the opposite phasing is less than the threshold, phase block will be broken (|Σweight| < threshold). \n" +
                    "[default = 0]")

args = parser.parse_args()

out = open(args.out, "w+")
num_vcf = len(args.vcf)
threshold = float(args.threshold)

if args.weights == None:
    weights = [1] * num_vcf
elif len(args.weights.split(",")) != num_vcf:
    print("ERROR: inputted weights does not correspond to the number of inputted VCFs.")
    print("{} weight value inputted {}; with {} VCFs".format(len(args.weights.split(",")), args.weights, num_vcf))
    print("weights: {}".format(args.weights.split(",")))
    exit(0)
else:
    weights = [float(w) for w in args.weights.split(",")]
    
assert all(w >= 0 for w in weights), "Cannot have weights less than 0; inputted weights: {}".format(weights)

if args.names:
    names = args.names.split(",")
    if len(names) != num_vcf:
        print("ERROR: number of names does not correspond to the number of inputted VCFs")
        print("number of names: {} ; number of VCFs: {}".format(len(names), num_vcf))
        print("names: {}".format(names))
        exit(0)
else:
    names = []
    for i in range(0, num_vcf):
        names.append("file" + str(i))

vcfs = []
for vcf in args.vcf:
    vcfs.append(vcf)

def verbrose():
    print("\n")
    print("Inputed parameters")
    print("-"*50)
    print("Current working directory: {}".format(os.getcwd()))
    print("Output consensus VCF: {}".format(args.out))
    print("Threshold: {}".format(threshold))
    print("Number of inputted VCFs: {}\n".format(num_vcf))
    for n in range(num_vcf):
        print("VCF: {}".format(vcfs[n]))
        print("Name: {}".format(names[n]))
        print("Weight: {}\n".format(weights[n]))
    print("-"*50)

#skip header lines for vcfs
#if write_header == True; copy the header to the output file; this is done to the first inputted VCF
def skip_header_lines(vcf_file, write_header):

  if vcf_file[-2:] == "gz" or vcf_file[-2:] == "GZ":
    vcf =  gzip.open(vcf_file, 'rt')
  else:
    vcf = open(vcf_file, "r")

  included_new_info = False
  for line in vcf:
    if line[:2] == "##":
      if write_header:
        if "##INFO" in line and not included_new_info:
          out.write(line)
          out.write('##INFO=<ID=phase_support,Number=.,Type=String,Description="Datasets that supports the phasing of this site.">\n')
          out.write('##INFO=<ID=anti_phase_support,Number=.,Type=String,Description="Datasets that supports the opposite phasing of this site.">\n')
          included_new_info = True
        else:
          out.write(line)
    else:
      if write_header:
        out.write('##commandline="{}"\n'.format(str(" ".join(sys.argv))))
        out.write(line)
      return vcf


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
'''
last_ps_call = []
for i in range(0, num_vcf):
    last_ps_call.append({})
    '''

#returns two list
# 1) list of the ps value from the sample column
# 2) list of the index to find the PS in the sample column
# if a dataset does not have PS, it's filled in as False
def ps_parse(form, sample):
    ind = 0
    ps = []
    ps_ind = []
    while ind < num_vcf:
        if "PS" in form[ind]:
          ps_format_ind = form[ind].index("PS")
          ps_ind.append(ps_format_ind)
          ps.append(sample[ind][ps_format_ind])
        else:
          ps.append(False)
          ps_ind.append(False)
        ind += 1
    return ps, ps_ind

# determines if two phased genotypes are cis or trans
# returns two variables; phase orientation (cis or trans) and if one of the two genotypes are multiallelic
# phase orientation: 1 if the two alt vars are cis, -1 if they are trans
def cis_trans_phase(gt, prev_gt):
    #determine if either of the genotypes are multiallelelics
    if "2" in gt or "2" in prev_gt:
        multiallelic = True
    else:
        multiallelic = False

    if prev_gt == "" or prev_gt[1] == "/" or gt[1] == "/":
        return False, multiallelic
    elif gt == prev_gt:
        return 1, multiallelic
    #alt_var trans to prev_var
    elif gt == prev_gt[::-1]:
        return -1, multiallelic
    #no prev gt
    elif prev_gt == "":
        return False, multiallelic

    #multi allelic, phases as if the 2nd allele is the refenence
    #current var is multiallelic
    if "2" in gt and "1" in gt:
        #cis
        if gt.index("1") == prev_gt.index("1"):
          return 1, multiallelic
        #trans
        elif gt.index("1") == prev_gt.index("2"):
          return -1, multiallelic
    #last gt was multi allelic
    elif "2" in prev_gt and "1" in prev_gt:
        #cis
        if gt.index("1") == prev_gt.index("1"):
          return 1, multiallelic
        #trans
        elif gt.index("1") == prev_gt.index("2"):
          return -1, multiallelic
    else:
      print("unphased multiallelic")


#main function
def consensus_phase(vcf_list):
  verbrose()

  vcfs = [skip_header_lines(vcf_list[0], True)] #writing this header to the output
  vcfs += [skip_header_lines(vcf, False) for vcf in vcf_list[1:]]

  #set up variables
  last_ps = ["" for x in vcf_list] #ps of the last hetero var
  last_ps_gt = ["" for x in vcf_list] #gt of the last hetero var
  votes = [0 for x in vcf_list] #alt same = 1, alt trans -1
  multiallelic = ["" for x in vcf_list]
  last_var_gt = ""
  var_gt = ""
  last_chromo = ""
  dict_last_gt_ps = [{} for x in vcf_list] #dictionary of the last gt in this ps

  for datasets in zip(*vcfs):
    pos = [data.split("\t")[1] for data in datasets]
    # check VCF is in the right order
    # only checking POS, not checking CHROMO, REF and ALT
    assert all(x == pos[0] for x in pos), "ERROR VCFs variant positions are not the same: {}".format(pos)

    sample = [data.split("\t")[9].strip('\n').split(":") for data in datasets]
    form = [data.split("\t")[8].split(":") for data in datasets]
    gt = [gt[gt_ind.index("GT")] for gt,gt_ind in zip(sample, form)]
    ps, ps_ind = ps_parse(form, sample) #parse for the PS column

    '''heterozygous variant'''
    if gt[0][0] != gt[0][2]:

        is_phased_presumptive = [gt_phase[1] == "|" for gt_phase in gt] #list of bool
        is_phased = []
        #remove the phase from the list and keep only ones that matter (based on weights)
        if 0 in weights:
          for w, phase in zip(weights, is_phased_presumptive):
            if w == 0:
              is_phased.append(False)
            else:
              is_phased.append(phase)
        else:
          is_phased = is_phased_presumptive

        #if at least one dataset phases this variant:
        if True in is_phased:
              ind = 0
              chromo = datasets[0].split("\t")[0]
              '''make the votes'''
              while ind < num_vcf:

                if is_phased[ind]:
                  #same phase set as the previous hetero var, gets to vote on the var GT
                  if ps[ind] == last_ps[ind] and last_chromo == chromo:
                      votes[ind], multiallelic[ind] = cis_trans_phase(gt[ind], last_ps_gt[ind])
                  #new chromosome, doesnt get to vote
                  elif last_chromo != chromo:
                      last_var_ps = [{} for x in vcf_list]
                      votes[ind] = 0
                      last_ps_gt[ind] = ""
                      last_ps[ind] = ""
                  #ps different from last hetero var or last hetero var was unphased
                  elif ps[ind] != last_ps[ind]:
                      votes[ind] = 0

                elif not is_phased[ind]: #not phased
                  votes[ind] = False

                ind += 1


              #all datasets could not make a vote, but there is phase information in the datasets
              if all(x == 0 for x in votes):
                #try to fill in phase information only if the current variant cannot be phased with respect to the previous hetero variant
                #check if ps skipped variants
                sorted_weights = sorted(weights, reverse=True)
                temp_weights = weights.copy()
                ind = 0 
                while ind < num_vcf:
                  #start with one with the heaviest weight, if tied in weight it is based on input order
                  top_weight_ind = temp_weights.index(sorted_weights[ind])
                  if is_phased[top_weight_ind] and ps[top_weight_ind] in dict_last_gt_ps[top_weight_ind]:
                    last_ps_gt[top_weight_ind] = dict_last_gt_ps[top_weight_ind][ps[top_weight_ind]] 
                    last_ps[top_weight_ind] = ps[top_weight_ind]
                    votes[top_weight_ind], multiallelic[top_weight_ind] = cis_trans_phase(gt[top_weight_ind], last_ps_gt[top_weight_ind])
                    break
                  #if there is a tie, make sure we move on to next dataset
                  temp_weights[top_weight_ind] = 0 
                  ind += 1
                
                if all(x == 0 for x in votes):
                    #check if there is new phase information
                    #first var of the chromo or new phase block for all datasets
                    if all(x != y or x == False or x == "." for x,y in zip(ps, last_ps)) or last_chromo != chromo:
                      var_gt = gt[0][0] + "|" + gt[0][2] #assign arbitray for new phase block
                      var_ps = pos[0]
                      support = []
                      anti_support = []
                    else:
                      #print("ERROR chr: {} pos: {} var is unphased but datasets has phase information".format(chromo, pos[0]))
                      pass
                  


              #weigh the votes
              else:
                weighted_votes = [a*b for a,b in zip(weights, votes)]
                '''summation of weighted votes less than threshold'''
                # variant is unphased
                if abs(sum(weighted_votes)) < threshold:
                  print("vote less than threshold at {}:{}".format(chromo, pos[0]))
                  var_gt = gt[0][0] + "/" + gt[0][2]
                  var_ps = ""
                  support = []
                  anti_support = []

                  '''summation of weighted votes equal or above threshold'''
                elif abs(sum(weighted_votes)) >= threshold:
                    #cis genotype
                    if sum(weighted_votes) > 0:
                        if True in multiallelic:
                            var_gt = ["x", "|", "x"]
                            var_gt[last_var_gt.index("1")] = "1"
                            if "0" in gt[0]:
                                var_gt[var_gt.index("x")] = "0"
                            elif "2" in gt[0]:
                                var_gt[var_gt.index("x")] = "2"
                            var_gt = "".join(var_gt)
                        else:
                          var_gt = last_var_gt

                        var_ps = last_var_ps
                        #datasets that agree/disagree with the choice of phasing
                        support = [x > 0 for x in weighted_votes]
                        anti_support = [x < 0 for x in weighted_votes]

                    #trans genotype
                    elif sum(weighted_votes) < 0:
                        if True in multiallelic:
                            var_gt = ["x", "|", "x"]
                            if last_var_gt.index("1") == 0:
                                var_gt[2] = "1"
                            elif last_var_gt.index("1") == 2:
                                var_gt[0] = "1"
                            if "0" in gt[0]:
                                var_gt[var_gt.index("x")] = "0"
                            elif "2" in gt[0]:
                                var_gt[var_gt.index("x")] = "2"
                            var_gt = "".join(var_gt)
                        else:
                          var_gt = last_var_gt[::-1]
                        var_ps = last_var_ps
                        support = [x < 0 for x in weighted_votes]
                        anti_support = [x > 0 for x in weighted_votes]

                    #tie in votes
                    elif sum(weighted_votes) == 0 :
                        print("variant at chr: {} pos: {} unphased due to tie".format(chromo, pos[0]))
                else:
                    print("threshold issue")



              '''write datasets that supporting phasing'''
              info_write = datasets[0].split("\t")[7]
              if True in support:
                info_write += ";phase_support="
                ind = 0
                while ind < num_vcf:
                  if support[ind]:
                    info_write += names[ind] + ","
                  ind += 1
                info_write = info_write[:-1]


              if True in anti_support :
                info_write += ";anti_phase_support="
                ind = 0
                while ind < num_vcf:
                  if anti_support[ind]:
                    info_write += names[ind] + ","
                  ind += 1
                info_write = info_write[:-1]



              format_write = form[0]
              sample_write = sample[0]
              #edit the PS if variant was sucessfully consensus phased
              if var_gt[1] == "|":
                #doesn't have PS already in format column
                if ps[0] == False:
                   format_write.append("PS")
                   sample_write.append(var_ps)
                else:
                  sample_write[ps_ind[0]] = var_ps

              #update the genotype
              sample_write[format_write.index("GT")] = var_gt
              format_write = ":".join(format_write)
              sample_write = ":".join(sample_write)


              #remember the previous variant to phase the next variant
              last_ps_gt = gt
              last_ps = ps
              last_var_gt = var_gt
              last_var_ps = var_ps
              last_chromo = chromo
              #record the last observed gt of the variant in this ps
              ind = 0
              while ind < num_vcf:
                if ps[ind]:
                  dict_last_gt_ps[ind][ps[ind]] = gt[ind]
                else:
                  pass
                ind += 1


              out.write("\t".join(datasets[0].split("\t")[:7]) + "\t" + info_write + "\t" + format_write + "\t" + sample_write + "\n")

        #all datasets did not phase this variant
        else:
          #change the format column
          sample_write = sample[0]
          format_write = form[0]
          if ps[ps_ind[0]] != False: #has PS in format column
            sample_write[ps_ind[0]] = "." #change PS value to .
          format_write = ":".join(format_write)

          if len(sample_write) == 1: #only one parameter in sample column
            sample_write = gt[0][0] + "/" + gt[0][2] + ":".join(sample[0][1:])
          else:
            sample_write = gt[0][0] + "/" + gt[0][2] + ":" + ":".join(sample[0][1:])

          out.write("\t".join(datasets[0].split("\t")[:8]) + "\t" + format_write + "\t" + sample_write + "\n")


    #homozygous variant
    #only make changes to the sample column to remove the PS and change gt to unphased
    else:
      sample_write = sample[0]
      gt = sample_write[form[0].index("GT")]

      sample_write[form[0].index("GT")] = gt[0] + "/" + gt[2] #change the GT to unphased

      if ps[0] != False:
        sample_write[ps_ind[0]] = "." #remove the PS

      out.write("\t".join(datasets[0].split("\t")[:9]) + "\t" + ":".join(sample_write) + "\n")
      continue

  out.close()


consensus_phase(vcfs)
