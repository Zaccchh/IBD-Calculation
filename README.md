# CSE284-IBD-Calculation

**Group 29**  
Zachariah Jasser  
Rafferty Chen

## Introduction

A tool for calculating IBD probabilities between samples.

## Installation

To install, follow these steps:

1. Clone this repository to your local machine:
    ```bash
    git clone https://github.com/Zaccchh/IBD-Calculation.git
    cd IBD-Calculation
    ```
2. Download sample files and install requirements:
   ```bash
   ./setup.sh
   ```

## Running 

To run our script, use the following command in a terminal: 
```text
python relative_finding.py <path_to_vcf_file> <sample1> <sample2>
```
\
If you don't know what samples are available, you can leave out the sample arguments, and the script will generate a text file called `sample_list.txt` that contains the list of samples for the specified VCF file. For example, the command:
```text
python relative_finding.py sample/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
```
will generate a list of sample names that looks like so:
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
python relative_finding.py sample/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz HG00096 HG00097
```

When running, our script will perform 2 passes through the VCF file to first calculate any prior probabilities needed for the IBD calculations and then to get the IBS data for the two samples. The step to calculate probabilities takes a long time to run (about 20 minutes), so we cache these probabilities per VCF file so that we do not need to calculate them multiple times. We've included precalculated probabilities for the sample file downloaded by our setup script, so you do not need to wait for the probabilities to calculate. 

## Results

Due to the time it takes for our tool to run, we could not calculate IBD for every pair of samples in our sample VCF file. In our `analysis.ipynb` notebook, we show how we sample 1500 pairs from the sample VCF file and run our tool on them. Our tool was able to acheive a mean squared error of `0.0006712319403097077`, which suggests we achieve IBD probabilities close to what are calculated with plink. 
