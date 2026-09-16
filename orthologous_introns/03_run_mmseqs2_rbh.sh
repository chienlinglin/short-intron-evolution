#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Configuration
# ============================================================

ROOT=${1:-"."}

REFERENCE_DIR="${ROOT}/extractedJunctionPeptides/Human"
PAIRED_ID_DIR="${ROOT}/extractedJunctionPeptides/PairedIDTable"

OUTPUT_ROOT="${ROOT}/MMseqs2_pairwise_output"
TMP_ROOT="${OUTPUT_ROOT}/tmp"

MMSEQS="mmseqs"

THREADS=4
PARALLEL_JOBS=5

mkdir -p "${OUTPUT_ROOT}" "${TMP_ROOT}"


# ============================================================
# Function for one human-species comparison
# ============================================================

run_rbh() {

    SPECIES=$1
    CSV_NAME=$2

    QUERY_DIR="${ROOT}/extractedJunctionPeptides/${SPECIES}"
    RESULT_DIR="${OUTPUT_ROOT}/Human_${SPECIES}"
    CSV_FILE="${PAIRED_ID_DIR}/${CSV_NAME}"

    mkdir -p "${RESULT_DIR}"
    mkdir -p "${TMP_ROOT}/${SPECIES}"

    echo "Running Human vs ${SPECIES}"

    tail -n +2 "${CSV_FILE}" |
    while IFS=',' read -r HUMAN_ID _ _ _ QUERY_ID _; do

        HUMAN_ID=$(echo "${HUMAN_ID}" | tr -d '"')
        QUERY_ID=$(echo "${QUERY_ID}" | tr -d '"')

        REFERENCE_FILE="${REFERENCE_DIR}/${HUMAN_ID}.fasta"
        QUERY_FILE="${QUERY_DIR}/${QUERY_ID}.fasta"

        RESULT_FILE="${RESULT_DIR}/${HUMAN_ID}_${QUERY_ID}_RBH.tsv"

        TMP_DIR="${TMP_ROOT}/${SPECIES}/${HUMAN_ID}_${QUERY_ID}"

        if [[ -f "${RESULT_FILE}" ]]; then
            continue
        fi

        if [[ ! -f "${REFERENCE_FILE}" ]]; then
            echo "Missing reference: ${REFERENCE_FILE}" >&2
            continue
        fi

        if [[ ! -f "${QUERY_FILE}" ]]; then
            echo "Missing query: ${QUERY_FILE}" >&2
            continue
        fi

        echo \
        "${MMSEQS} easy-rbh \
        ${REFERENCE_FILE} \
        ${QUERY_FILE} \
        ${RESULT_FILE} \
        ${TMP_DIR} \
        --min-seq-id 0.4 \
        --threads ${THREADS}"

    done |
    parallel -j "${PARALLEL_JOBS}"
}


# ============================================================
# Human vs each vertebrate species
# ============================================================

run_rbh \
    "Mouse" \
    "Mouse_vs_human.csv"

run_rbh \
    "Chicken" \
    "Chicken_vs_human.csv"

run_rbh \
    "Turtle" \
    "Western_painted_turtle_vs_human.csv"

run_rbh \
    "Frog" \
    "Western_clawed_frog_vs_human.csv"

run_rbh \
    "Zebrafish" \
    "Zebrafish_vs_human.csv"


echo "All MMseqs2 RBH comparisons completed."