# CSE284-IBD-Calculation

A tool for calculating IBD probabilities between samples (only 2 at a time for now).

## Running Locally

Running this script requires Python 3.8+ with `numpy` and `cyvcf2` installed. This can be done with the included `requirements.txt` file or manually with pip. You will also need a VCF file and two samples from that file that you would like to calculate the IBD for. Since the VCF files are too large to upload directly to GitHub, you will have to download them yourself. We downloaded the 1000Genomes files located on the class datahub to test our script, so we would suggest using those as well. 

To run our script, you would use the following command in a terminal: 
```text
python relative_finding.py <path_to_vcf_file> <sample1> <sample2>
```
\
If you don't know what samples are available, you can leave out the sample arguments, and the script will generate a text file called `sample_list.txt` that contains the list of samples for the specified VCF file. For example, the command:
```text
python relative_finding.py ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
```
will generate the a list that looks like so:
```text
HG00096
HG00097
HG00099
HG00100
HG00101
HG00102
...
```
\
A filled-out example would be: 
```text
python relative_finding.py ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz HG03750 HG03754
```
For comparison, plink gives an output of `P(IBD=0)=0.0147`, `P(IBD=1)=0.8491`, and `P(IBD=2)=0.1362` for these two samples.

When running, our script will perform 2 passes through the VCF file to first calculate any prior probabilities needed for the IBD calculations and then to get the IBS data for the two samples. The step to calculate probabilities takes a long time to run (about 20 minutes), so we cache these probabilities per VCF file so that we do not need to calculate them multiple times. We've included the cached probabilities for the `ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz` file, so we suggest using this file too if you'd like to avoid waiting. Without the probability calculation, the script should take 2-3 minutes to run (this depends on the length of the VCF file of course).

## Running on Datahub

If you'd like to run this script on Datahub, the same instructions above will apply. Just upload the relevant files and run the script through the terminal. Since there are already VCF files available on Datahub, there's no need to download them. 
