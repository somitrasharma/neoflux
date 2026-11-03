#!/usr/bin/env python3
"""Convert MAF to NeoFlux candidate TSV for scoring."""
import csv
import sys
import json
import uuid
from pathlib import Path

def maf_to_candidates(maf_path: str, sample_ids: list[str], output_path: str):
    rows = []
    with open(maf_path) as f:
        for line in f:
            if not line.startswith("#"):
                break
        header = line.rstrip("\n").split("\t")
        reader = csv.DictReader(f, fieldnames=header, delimiter="\t")
        for row in reader:
            sid = row.get("Tumor_Sample_Barcode", "")
            if sid not in sample_ids:
                continue
            vc = row.get("Variant_Classification", "")
            if vc not in ("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins",
                          "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del",
                          "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation",
                          "De_novo_Start_OutOfFrame", "Start_Codon_SNP", "Stop_Codon_SNP"):
                continue

            chrom = row.get("Chromosome", "").replace("chr", "")
            pos = row.get("Start_Position", "")
            ref = row.get("Reference_Allele", "")
            alt = row.get("Tumor_Seq_Allele2", "") or row.get("Tumor_Seq_Allele1", "")
            gene = row.get("Hugo_Symbol", "")
            protein = row.get("Protein_Change", "") or row.get("HGVSp_Short", "")

            var_type = "SNV"
            if vc in ("Frame_Shift_Ins", "In_Frame_Ins"):
                var_type = "INS"
            elif vc in ("Frame_Shift_Del", "In_Frame_Del"):
                var_type = "DEL"

            if "-" in chrom or chrom in ("M", "MT"):
                continue
            if not ref or not alt or not pos:
                continue
            if ref == alt:
                continue
            if len(ref) > 50 or len(alt) > 50:
                continue

            cid = f"NF-{uuid.uuid4().hex[:12]}"
            rows.append({
                "candidate_id": cid,
                "sample_id": sid,
                "gene": gene,
                "chromosome": chrom,
                "position": pos,
                "ref": ref,
                "alt": alt,
                "variant_type": var_type,
                "consequence": vc,
                "protein_change": protein,
            })

    with open(output_path, "w") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} candidates for {len(sample_ids)} samples to {output_path}")

if __name__ == "__main__":
    maf_path = sys.argv[1]
    output_path = sys.argv[2]
    samples = sys.argv[3:] if len(sys.argv) > 3 else ["ME001", "ME002", "ME009"]
    maf_to_candidates(maf_path, samples, output_path)
