#!/usr/bin/env bash

set -euo pipefail

# OrthoFinder v2.5.5
# Usage:
# bash scripts/01_run_orthofinder.sh /path/to/ensemblProteinSequence

PROTEIN_DIR=${1:-"ensemblProteinSequence"}

orthofinder -f "${PROTEIN_DIR}"