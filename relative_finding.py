import argparse
import time
import numpy as np

from cyvcf2 import VCF
from pathlib import Path


def get_probs(filename, debug=False):
    records = VCF(filename)
    ps = np.zeros((9,))

    start = time.time()

    print("Calculating probs at SNP: ")

    for count, record in enumerate(records):
        print(count, end="\r", flush=True)

        if len(record.ALT) > 1:
            continue

        p = record.aaf
        q = 1 - p

        # TODO: add sample size correction terms
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


def get_ibs(filename, sample1, sample2, debug=False):
    records = VCF(filename)
    records.set_samples([sample1, sample2])

    ibs = np.zeros((3,), dtype=int)

    start = time.time()

    print("Calculating IBS at SNP: ")

    for count, record in enumerate(records):
        print(count, end="\r", flush=True)

        if len(record.ALT) > 1:
            continue

        index = 0

        sorted_gt_0 = sorted(record.genotypes[0])
        sorted_gt_1 = sorted(record.genotypes[1])

        if sorted_gt_0 == sorted_gt_1:
            index = 2
        elif sorted_gt_0[0] in sorted_gt_1 or sorted_gt_0[1] in sorted_gt_1:
            index = 1

        ibs[index] += 1

    print()

    if debug:
        print(ibs)
        print(f"calculated in {time.time() - start:.2f} seconds")

    return ibs


def get_ibd(ps, ibs):
    ibd0 = ibs[0] / ps[0]
    ibd1 = (ibs[1] - ps[1] * ibd0) / ps[4]
    ibd2 = (ibs[2] - ps[2] * ibd0 - ps[5] * ibd1) / ps[8]

    if np.sum(np.array([ibd0, ibd1, ibd2]) > 1) == 1:
        ibd0, ibd1, ibd2 = np.array(np.array([ibd0, ibd1, ibd2]) > 1, dtype=int)

    # normalize
    pi = ibd1 / 2 + ibd2

    if pi**2 <= ibd2:
        ibd0_star = 1 - pi**2
        ibd1_star = 2 * pi * (1 - pi)
        ibd2_star = pi**2
    else:
        ibd0_star = ibd0
        ibd1_star = ibd1
        ibd2_star = ibd2

    return ibd0_star, ibd1_star, ibd2_star


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="CSE 284 IBD Calculator",
        description="Given a VCF file and two sample names located in the file, calculate the IBD between the two samples",
        epilog="ex) python relative_finding.py vcf/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz HG03750 HG03754",
    )

    parser.add_argument("filename")
    parser.add_argument("sample1")
    parser.add_argument("sample2")

    args = parser.parse_args()
    # print(args.filename, args.sample1, args.sample2)

    vcf_file = Path(args.filename)

    assert vcf_file.exists(), "Please input a valid VCF file"

    if (cache := Path(f"cache/{vcf_file.name}.npy")).exists():
        ps = np.load(cache)
    else:
        print("Probability cache not found, performing two passes.")
        ps = get_probs(vcf_file, debug=False)
        np.save(cache, ps)

    ibs = get_ibs(vcf_file, args.sample1, args.sample2, debug=False)

    p0, p1, p2 = get_ibd(ps, ibs)

    print(f"P(IBD=0) = {p0:.5f}")
    print(f"P(IBD=1) = {p1:.5f}")
    print(f"P(IBD=2) = {p2:.5f}")
