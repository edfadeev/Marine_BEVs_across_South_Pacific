#!/usr/bin/env python3
import argparse
import gzip
import re
import sys

def tryptic_cuts(seq):
    # Cut after K or R unless followed by P
    idx = [0]
    L = len(seq)
    for i, aa in enumerate(seq):
        if aa in "KR" and (i + 1 == L or seq[i + 1] != "P"):
            idx.append(i + 1)
    if idx[-1] != L:
        idx.append(L)
    return idx

def tryptic_peptides(seq, missed=2, min_len=7, max_len=35):
    cuts = tryptic_cuts(seq)
    peps = []
    n = len(cuts) - 1
    for i in range(n):
        for m in range(missed + 1):
            j = i + 1 + m
            if j < len(cuts):
                pep = seq[cuts[i]:cuts[j]]
                if min_len <= len(pep) <= max_len:
                    peps.append(pep)
    return peps

def sanitize_header(h):
    # Keep first token; replace spaces and unusual chars with underscore
    h = h.strip()
    tok = h.split()[0]
    tok = re.sub(r"\s+", "_", tok)
    tok = re.sub(r"[^A-Za-z0-9|_.:-]", "_", tok)
    return tok

def normalize_seq(s):
    # Remove stops, uppercase, replace non-IUPAC AAs with X
    s = s.replace("*", "").upper()
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "X", s)

def write_record(fh, header, seq, width=60):
    fh.write(f">{header}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i+width] + "\n")

def fasta_reader(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fin:
        header, seq_chunks = None, []
        for line in fin:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)

def parse_args():
    ap = argparse.ArgumentParser(
        description="Filter proteins by in-silico tryptic peptide content (no explicit protein-length filter)."
    )
    ap.add_argument("inp", help="Input FASTA (optionally .gz)")
    ap.add_argument("out", help="Output FASTA")
    ap.add_argument("--min_peptides", type=int, default=1,
                    help="Minimum number of qualifying tryptic peptides to keep a protein. Default: 1")
    ap.add_argument("--pep_min", type=int, default=7,
                    help="Minimum tryptic peptide length. Default: 7")
    ap.add_argument("--pep_max", type=int, default=35,
                    help="Maximum tryptic peptide length. Default: 35")
    ap.add_argument("--missed", type=int, default=2,
                    help="Missed cleavages allowed. Default: 2")
    return ap.parse_args()

def main():
    args = parse_args()
    total, kept = 0, 0
    seen_ids = set()

    with open(args.out, "w") as fout:
        for raw_header, raw_seq in fasta_reader(args.inp):
            total += 1
            s = normalize_seq(raw_seq)
            peps = tryptic_peptides(s, missed=args.missed,
                                    min_len=args.pep_min, max_len=args.pep_max)
            if len(peps) < args.min_peptides:
                continue
            hid = sanitize_header(raw_header)
            # ensure unique IDs
            base = hid
            if hid in seen_ids:
                n = 2
                while f"{base}|dup{n}" in seen_ids:
                    n += 1
                hid = f"{base}|dup{n}"
            seen_ids.add(hid)
            write_record(fout, hid, s)
            kept += 1

    print(f"Tryptic filter: processed {total} records → kept {kept}", file=sys.stderr)

if __name__ == "__main__":
    main()
