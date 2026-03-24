import os
import argparse

## Extract DARLIN sequence (Reads2)
## -j threads
## -A Adapter | 3' adapter to be removed from second read in a pair [R2].
## -G Adapter | 5' adapter to be removed from second read in a pair [R2].
## -q cutoff | Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. 
##           | Set to 0 means do not perform this trimming.
## -m cutoff | Minimum length after trimming.
## --discard-untrimmed | [Important!!] Discard reads that do not contain an adapter.
## partial match
# cmd1 = f"cutadapt -j {cores} -A {prime3} -q {base_quality} -m 10 --discard-untrimmed -o {tmp_fq1} -p {tmp_fq2} {in_fq1} {in_fq2}"
# cmd2 = f"cutadapt -j {cores} -G {prime5} -q {base_quality} -m 10 --discard-untrimmed -o {out_fq1} -p {out_fq2} {tmp_fq1} {tmp_fq2}"

PRIMERS = {
    "CA" : {
        "5prime": "AGCTGTACAAGTAAGCGGC",
        "3prime": "AGAATTCTAACTAGAGCTCGCTGATCAGCCTCGACTGTGCCTTCT"
    },
    "RA" : {
        "5prime":"GTACAAGTAAAGCGGCC",
        "3prime":"GTCTGCTGTGTGCCTTCTAGTT"
    },
    "TA" : {
        "5prime": "TCGGTACCTCGCGAA",
        "3prime": "GTCTTGTCGGTGCCTTCTAGTT"
    }
}

primer_len = 15

def cut(templete, cores, reads1, reads2, output, base_quality, sample):

    prime3 = PRIMERS[templete]['3prime']
    prime5 = PRIMERS[templete]['5prime']

    prime3 = prime3[:primer_len] if len(prime3) > primer_len else prime3
    prime5 = prime5[-primer_len:] if len(prime5) > primer_len else prime5
    len_p3 = len(prime3)
    len_p5 = len(prime5)

    in_fq1 = reads1
    in_fq2 = reads2
    tmp_fq1 = f"{output}/{sample}_R1.tmp.fq.gz"
    tmp_fq2 = f"{output}/{sample}_R2.tmp.fq.gz"
    out_fq1 = f"{output}/{sample}_R1.trimmed.fq.gz"
    out_fq2 = f"{output}/{sample}_R2.trimmed.fq.gz"

    cmd1 = f"cutadapt -j {cores} -A {prime3} --overlap {len_p3} --error-rate 0 -q {base_quality} -m 10 --discard-untrimmed -o {tmp_fq1} -p {tmp_fq2} {in_fq1} {in_fq2}"
    cmd2 = f"cutadapt -j {cores} -G {prime5} --overlap {len_p5} --error-rate 0 -q {base_quality} -m 10 --discard-untrimmed -o {out_fq1} -p {out_fq2} {tmp_fq1} {tmp_fq2}"
    os.system(cmd1)
    if os.system(cmd2) == 0:  # Check if cmd2 completed successfully
        os.system(f"rm {tmp_fq1} {tmp_fq2}")

    return out_fq1, out_fq2

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Cut adaptor sequences from reads')

    parser.add_argument('-l', '--locus', type=str, help='CA, RA, TA')
    parser.add_argument('-t', '--threads', type=int, default = 8, help='Number of threads to use (default: 8)')
    parser.add_argument('-r1', '--reads1', type=str, help='Input directory of raw fastq files')
    parser.add_argument('-r2', '--reads2', type=str, help='Input directory of raw fastq files')
    parser.add_argument('-o', '--output', type=str, help='Output directory of trimmed fastq files')
    parser.add_argument('-q', '--base_quality', type=int, default=10, help='Base quality cutoff (default: 10)')
    parser.add_argument('-s', '--sample', type=str, help='Sample name, prefix of output files')
    
    args = parser.parse_args()

    templete = args.locus
    cores = int(args.threads)
    reads1 = args.reads1
    reads2 = args.reads2
    output = args.output
    base_quality = args.base_quality if args.base_quality else 0
    sample = args.sample

    out_fq1, out_fq2 = cut(templete, cores, reads1, reads2, output, base_quality, sample)