To run an example of the phase consensus script for NA12878 chr12 

```
python consensus_vcf_phasing.py -n 10x,nanopore,ccs,clr -w 1.2,1,1.1,1.05 -t 0 \
   example/na12878_longranger_chr12.vcf.gz example/na12878_pacbio_ccs_chr12.vcf.gz \
   example/na12878_pacbio_ccs_chr12.vcf.gz example/na12878_pacbio_clr_chr12.vcf.gz

