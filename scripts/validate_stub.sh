#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

if ! command -v nextflow >/dev/null 2>&1; then
  echo "nextflow is not installed. Install with: brew install nextflow"
  exit 1
fi

make install
make test
make train-demo

nextflow run main.nf -profile test,stub -stub-run

echo "NeoFlux stub validation completed successfully."
