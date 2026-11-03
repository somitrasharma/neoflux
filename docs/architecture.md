# NeoFlux Architecture Notes

NeoFlux is organized as a Nextflow DSL2 workflow with small process wrappers and
a Python package for data models, feature engineering, evidence integration, and
ranking.

## Module Boundaries

- Nextflow owns orchestration, provenance, containers, and file-level contracts.
- Python owns biological data models, normalization, peptide representation,
  feature engineering, scoring, and cohort summaries.
- External tools stay replaceable. For example, MHCflurry can be used as the
  default class I predictor while NetMHCpan remains an optional licensed plug-in.

## Main Contracts

- HLA alleles: TSV with at least `allele` and optional `hla_class`.
- Candidate neoantigens: JSONL, one `NeoantigenCandidate` per line.
- Peptides: FASTA plus TSV with `candidate_id`, `sample_id`, `peptide`.
- Binding predictions: TSV keyed by `candidate_id` and `hla_allele`.
- MS evidence: TSV keyed by `candidate_id` or peptide sequence.
- Feature matrix: TSV with `candidate_id`, `sample_id`, and ML feature columns.
- Ranked output: TSV from `NeoantigenCandidate.to_output_row()`.

## Realistic First Demo Scope

For a small academic setup, the strongest first milestone is:

1. Tumor/normal WES BAM plus tumor RNA-seq FASTQ for 5 to 20 public samples.
2. HLA allele list from OptiType or a public metadata table.
3. Mutect2 plus Strelka2 SNV/indel calls.
4. STAR-Fusion and rMATS enabled for RNA evidence.
5. MHCflurry as the open default predictor.
6. Optional FragPipe/PepQuery2 if mzML is available.
7. Baseline ranking compared against a trained model once labels are assembled
   from TESLA, IEDB, and HLA Ligand Atlas.

This scope is ambitious enough for a SITC abstract while avoiding a dependency
on proprietary software for the default path.
