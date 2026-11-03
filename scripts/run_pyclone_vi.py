#!/usr/bin/env python3
"""PyClone-VI wrapper for NeoFlux.

Installs via conda:
  conda install -c bioconda pyclone-vi

Usage:
  python run_pyclone_vi.py \
    --vaf-table <path> \
    --sample-id <id> \
    --output <path>

Input TSV columns: variant_id, vaf, depth
Output TSV columns: variant_id, clonal_status, cancer_cell_fraction

PyClone-VI requires:
  - VAF table with columns: mutation_id, sample_id, ref_counts, alt_counts, ... 1, 1, normal_cn, tumour_cn, tumour_content
  - Copy number and purity inputs

This wrapper converts the NeoFlux VAF table into PyClone-VI input format.
"""

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path


def build_parser():
    p = argparse.ArgumentParser()
    p.add_argument("--vaf-table", required=True)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--cnv-file", default="", help="Optional copy number file")
    p.add_argument("--purity", type=float, default=None, help="Tumor purity (0-1)")
    return p


def vaf_to_pyclone_input(vaf_rows: list[dict], sample_id: str) -> list[dict]:
    """Convert NeoFlux VAF table to PyClone-VI input rows.

    PyClone-VI expects per row:
      mutation_id, sample_id, ref_counts, alt_counts, ... 1, 1, normal_cn, tumour_cn, tumour_content
    """
    out = []
    for row in vaf_rows:
        vaf = float(row.get("vaf") or 0)
        depth = int(row.get("depth") or 0)
        alt = round(depth * vaf)
        ref = depth - alt
        out.append({
            "mutation_id": row.get("variant_id", ""),
            "sample_id": sample_id,
            "ref_counts": str(max(ref, 0)),
            "alt_counts": str(max(alt, 0)),
            "": "1",  # normal_cn used as positional
            "": "1",  # tumour_cn
            "": "1",  # tumour_content placeholder
        })
    return out


def main():
    args = build_parser().parse_args()

    with open(args.vaf_table) as f:
        reader = csv.DictReader(f, delimiter="\t")
        vaf_rows = list(reader)

    pyclone_input = vaf_to_pyclone_input(vaf_rows, args.sample_id)

    # Write temp input file
    tmp_input = Path(tempfile.mktemp(suffix=".tsv"))
    with open(tmp_input, "w") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "mutation_id", "sample_id", "ref_counts", "alt_counts",
        ], delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(pyclone_input)

    # Run PyClone-VI
    cmd = [
        "pyclone-vi", "fit",
        "-i", str(tmp_input),
        "-o", str(Path(args.output).with_suffix(".h5")),
        "-c", "40",
        "-d", "beta-binomial",
        "-r", "10",
    ]
    if args.purity:
        cmd.extend(["-p", str(args.purity)])

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        print("pyclone-vi not found. Install: conda install -c bioconda pyclone-vi")
        print("Falling back to heuristic clonality...")
        # Fallback: call pyclone-stub
        from neoflux.cli import command_pyclone_stub
        import argparse
        command_pyclone_stub(argparse.Namespace(
            vaf_table=args.vaf_table,
            output=args.output,
        ))
        return
    except subprocess.CalledProcessError as e:
        print(f"PyClone-VI failed: {e.stderr}")
        raise

    # Convert PyClone-VI output to NeoFlux format
    # PyClone-VI outputs .h5 file; convert via build_table subcommand
    build_cmd = [
        "pyclone-vi", "build_table",
        "-i", str(Path(args.output).with_suffix(".h5")),
        "-o", args.output,
    ]
    subprocess.run(build_cmd, check=True, capture_output=True, text=True)

    tmp_input.unlink()


if __name__ == "__main__":
    main()
