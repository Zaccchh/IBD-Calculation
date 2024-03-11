import argparse
import functools
import pysam
import subprocess
import time
import numpy as np
import multiprocessing as mp

from cyvcf2 import VCF
from pathlib import Path


def get_stats(filename):
    """
    uses bcftools to get contig name, contig length, and number of records
    assumes a .tbi file with the correct stats exists
    """
    bcftools_output = subprocess.run(
        ["bcftools", "index", "-t", "-s", filename], capture_output=True
    )
    try:
        ref, contig_length, num_snps = (
            bcftools_output.stdout.decode("utf-8").strip().split("\t")
        )
    except ValueError:
        raise ValueError("No valid tabix file exists")

    contig_length = int(contig_length)
    num_snps = int(num_snps)

    return ref, contig_length, num_snps


def get_probs(filename, debug=False):
    """calculates probabilities needed for IBD calculations"""
    records = VCF(filename, threads=mp.cpu_count() - 1)
    ps = np.zeros((9,))

    start = time.time()

    print("Calculating probs at SNP: ")

    for count, record in enumerate(records):
        if count % 2_000 == 0:
            print(f"{count}\t\t{time.time()-start:.2f}s", end="\r", flush=True)

        if len(record.ALT) > 1:
            continue

        p = record.aaf
        q = 1 - p

        ps[0] += 2 * p**2 * q**2
        ps[1] += 4 * p**3 * q + 4 * p * q**3
        ps[2] += p**4 + q**4 + 4 * p**2 * q**2
        ps[3] += 0
        ps[4] += 2 * p**2 * q + 2 * p * q**2
        ps[5] += p**3 + q**3 + p**2 * q + p * q**2
        ps[6] += 0
        ps[7] += 0
        ps[8] += 1

    print()

    if debug:
        print(ps)
        print(f"calculated in {time.time() - start:.2f} seconds")

    return ps


def get_sample_idx(filename):
    """returns a mapping from sample names to indices"""
    records = VCF(filename)
    return {k: v for v, k in enumerate(records.samples)}


def calc_ibs_segment(i1, i2, ref, start, end):
    """calculates IBS for a segment of a VCF file"""
    ibs = np.zeros((3,), dtype=int)

    global tb
    for record in tb.fetch(
        reference=ref,
        start=start,
        end=end,
        parser=pysam.asVCF(),
        multiple_iterators=True,
    ):
        if "," in record.alt:  # skip snps with multiple alt alleles
            continue

        sorted_gt_0 = sorted(record[i1].split("|"))
        sorted_gt_1 = sorted(record[i2].split("|"))

        if sorted_gt_0[0] == sorted_gt_1[0] and sorted_gt_0[1] == sorted_gt_1[1]:
            index = 2
        elif sorted_gt_0[0] != sorted_gt_1[0] and sorted_gt_0[1] != sorted_gt_1[1]:
            index = 0
        else:
            index = 1

        ibs[index] += 1

    return ibs


def get_ibs(ref, sample_idx, sample1, sample2, debug=False):
    print("Calculating IBS")

    i1 = sample_idx[sample1]
    i2 = sample_idx[sample2]

    start = time.time()

    with mp.Pool() as p:
        result = p.starmap(functools.partial(calc_ibs_segment, i1, i2, ref), regions)

    if debug:
        print(np.sum(result, axis=0))
        print(f"calculated in {time.time() - start:.2f} seconds")

    return np.sum(result, axis=0)


def get_ibd(ps, ibs):
    print("Calculating IBD")

    ibd0 = ibs[0] / ps[0]
    ibd1 = (ibs[1] - ps[1] * ibd0) / ps[4]
    ibd2 = (ibs[2] - ps[2] * ibd0 - ps[5] * ibd1) / ps[8]

    temp = np.array([ibd0, ibd1, ibd2])

    # normalize
    if np.sum(temp > 1) == 1:
        ibd0, ibd1, ibd2 = np.array(temp > 1, dtype=int)
    elif np.sum(temp < 0) == 1:
        value = temp[temp < 0]
        temp[temp < 0] = 0
        temp[temp.argmax()] = temp[temp.argmax()] + value[0]
        ibd0, ibd1, ibd2 = temp

    return ibd0, ibd1, ibd2


def get_samples(filename):
    """returns a list of sample names present in the specified VCF file"""
    records = VCF(filename)
    return records.samples


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="CSE 284 IBD Calculator",
        description="Given a VCF file and two sample names located in the file, calculate the IBD between the two samples",
        epilog="ex) python relative_finding.py vcf/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz HG03750 HG03754",
    )

    parser.add_argument("filename")
    parser.add_argument("sample1", nargs="?", default=None)
    parser.add_argument("sample2", nargs="?", default=None)

    args = parser.parse_args()

    vcf_file = Path(args.filename)
    assert vcf_file.exists(), "Please input a valid VCF file"

    ref, contig_length, _ = get_stats(vcf_file)

    cores = mp.cpu_count()
    segment = contig_length // cores
    regions = [[segment * i, segment * (i + 1) + 1] for i in range(cores)]
    regions[-1][-1] = None

    sample_idx = get_sample_idx(vcf_file)

    if args.sample1 is None and args.sample2 is None:
        with open("sample_list.txt", "w") as file:
            samples = get_samples(vcf_file)
            file.write("\n".join(samples))
    else:
        if (cache := Path(f"cache/{vcf_file.name}.npy")).exists():
            ps = np.load(cache)
        else:
            print("Probability cache not found, performing two passes.")
            ps = get_probs(vcf_file, debug=False)
            np.save(cache, ps)

    global tb  # sorry for use of global, but TabixFile is not picklable
    tb = pysam.TabixFile(str(vcf_file))

    ibs = get_ibs(ref, sample_idx, args.sample1, args.sample2, debug=True)

    p0, p1, p2 = get_ibd(ps, ibs)

    print(f"P(IBD=0) = {p0:.5f}")
    print(f"P(IBD=1) = {p1:.5f}")
    print(f"P(IBD=2) = {p2:.5f}")
