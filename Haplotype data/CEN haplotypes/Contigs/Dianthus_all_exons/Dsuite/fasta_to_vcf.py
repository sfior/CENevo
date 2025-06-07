#!/usr/bin/env python3

# source: https://github.com/pomo-dev/PoMo/blob/master/scripts/FastaToVCF.py


import argparse
from Bio import SeqIO

def is_valid_base(base):
    return base in "ACGT"

def write_vcf_header(handle, sample_names):
    handle.write("##fileformat=VCFv4.2\n")
    handle.write("##source=fasta_to_vcf.py\n")
    handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names) + "\n")

def convert_fasta_to_vcf(fasta_file, vcf_file, ref_name=None):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        raise ValueError("No sequences found in the FASTA file.")

    # Choose the reference sequence
    if ref_name:
        ref_record = next((r for r in records if r.id == ref_name), None)
        if not ref_record:
            raise ValueError(f"Reference sequence '{ref_name}' not found in the FASTA file.")
    else:
        ref_record = records[0]

    ref_seq = str(ref_record.seq).upper()
    sample_records = records  # include all individuals, including reference
    sample_names = [r.id for r in sample_records]

    with open(vcf_file, "w") as vcf_out:
        write_vcf_header(vcf_out, sample_names)
        for i in range(len(ref_seq)):
            pos = i + 1
            bases_at_pos = [str(r.seq[i]).upper() for r in sample_records]

            # Skip positions with any gaps or non-ACGT
            if any(b not in "ACGT" for b in bases_at_pos):
                continue

            ref_base = str(ref_seq[i])
            alt_alleles = sorted(set(b for b in bases_at_pos if b != ref_base))
            if not alt_alleles:
                continue  # monomorphic

            # Assign numeric genotype codes
            alt_map = {allele: str(idx + 1) for idx, allele in enumerate(alt_alleles)}
            genotypes = []
            for base in bases_at_pos:
                genotypes.append(alt_map.get(base, "0"))  # 0 if matches ref

            vcf_out.write(f"chr1\t{pos}\t.\t{ref_base}\t{','.join(alt_alleles)}\t.\tPASS\t.\tGT\t" +
                          "\t".join(genotypes) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert aligned FASTA to VCF, skipping gap positions.")
    parser.add_argument("fasta", help="Input FASTA file (alignment)")
    parser.add_argument("vcf", help="Output VCF file")
    parser.add_argument("-r", "--reference", help="Reference sequence ID (default: first sequence)")
    args = parser.parse_args()

    convert_fasta_to_vcf(args.fasta, args.vcf, args.reference)