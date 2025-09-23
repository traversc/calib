#!/usr/bin/env bash
set -euo pipefail

TMP_DIR="$(mktemp -d 2>/dev/null || mktemp -d -t calib-cons-smoke)"
cleanup() {
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

CLUSTER_FILE="${TMP_DIR}/input.cluster"
R1_FASTQ="${TMP_DIR}/reads_R1.fastq"
R2_FASTQ="${TMP_DIR}/reads_R2.fastq"
OUT_PREFIX_R1="${TMP_DIR}/consensus.R1."
OUT_PREFIX_R2="${TMP_DIR}/consensus.R2."
OUT_R1_FASTQ="${OUT_PREFIX_R1}fastq"
OUT_R2_FASTQ="${OUT_PREFIX_R2}fastq"

# Generate cluster definitions: 3 clusters with sizes 1, 1001, and 10.
{
  printf '0\t0\t0\n'
  for rid in $(seq 1 1001); do
    printf '1\t%d\t%d\n' "$rid" "$rid"
  done
  for rid in $(seq 1002 1011); do
    printf '2\t%d\t%d\n' "$rid" "$rid"
  done
} > "${CLUSTER_FILE}"

emit_fastq_record() {
  local file=$1
  local read_id=$2
  local mate=$3
  local sequence=$4
  local quality=$5
  printf '@read%04d/%d\n%s\n+\n%s\n' "${read_id}" "${mate}" "${sequence}" "${quality}" >> "${file}"
}

# Prepare FASTQ inputs with deterministic content.
: > "${R1_FASTQ}"
: > "${R2_FASTQ}"
for rid in $(seq 0 1011); do
  if [ "${rid}" -eq 0 ]; then
    r1_seq="ACGT"
    r1_qual="!!!!"
    r2_seq="TGCA"
    r2_qual="####"
  elif [ "${rid}" -le 1001 ]; then
    if [ "${rid}" -eq 1 ]; then
      r1_seq="AAAA"
      r2_seq="TTTT"
    else
      r1_seq="CCCC"
      r2_seq="GGGG"
    fi
    r1_qual="IIII"
    r2_qual="JJJJ"
  else
    if [ "${rid}" -eq 1002 ]; then
      r1_seq="AAAA"
      r2_seq="TTTT"
    else
      r1_seq="CCCC"
      r2_seq="GGGG"
    fi
    r1_qual="KKKK"
    r2_qual="LLLL"
  fi
  emit_fastq_record "${R1_FASTQ}" "${rid}" 1 "${r1_seq}" "${r1_qual}"
  emit_fastq_record "${R2_FASTQ}" "${rid}" 2 "${r2_seq}" "${r2_qual}"
done

# Run calib_cons with explicit format controls and cluster limits.
./consensus/calib_cons \
  --input-format plain \
  --cluster-input-format plain \
  --output-format plain \
  --min-reads-per-cluster 1 \
  --max-reads-per-cluster 2000 \
  -t 1 \
  -c "${CLUSTER_FILE}" \
  -q "${R1_FASTQ}" "${R2_FASTQ}" \
  -o "${OUT_PREFIX_R1}" "${OUT_PREFIX_R2}"

# Basic sanity checks on output structure.
if [ ! -f "${OUT_R1_FASTQ}" ] || [ ! -f "${OUT_R2_FASTQ}" ]; then
  echo "Expected consensus FASTQ outputs were not generated" >&2
  exit 1
fi

r1_records=$(grep -c '^@' "${OUT_R1_FASTQ}")
r2_records=$(grep -c '^@' "${OUT_R2_FASTQ}")
if [ "${r1_records}" -ne 3 ] || [ "${r2_records}" -ne 3 ]; then
  echo "Unexpected number of consensus records: R1=${r1_records}, R2=${r2_records}" >&2
  exit 1
fi

IFS=$'\n' r1_sequences=($(awk 'NR % 4 == 2 { print }' "${OUT_R1_FASTQ}"))
IFS=$'\n' r2_sequences=($(awk 'NR % 4 == 2 { print }' "${OUT_R2_FASTQ}"))
IFS=$' \t\n'

if [ "${r1_sequences[0]}" != "ACGT" ]; then
  echo "Cluster 0 R1 consensus sequence mismatch: ${r1_sequences[0]}" >&2
  exit 1
fi
if [ "${r2_sequences[0]}" != "TGCA" ]; then
  echo "Cluster 0 R2 consensus sequence mismatch: ${r2_sequences[0]}" >&2
  exit 1
fi
if [ "${r1_sequences[1]}" != "CCCC" ] || [ "${r1_sequences[2]}" != "CCCC" ]; then
  echo "Unexpected R1 consensus sequences for multi-read clusters: ${r1_sequences[*]}" >&2
  exit 1
fi
if [ "${r2_sequences[1]}" != "GGGG" ] || [ "${r2_sequences[2]}" != "GGGG" ]; then
  echo "Unexpected R2 consensus sequences for multi-read clusters: ${r2_sequences[*]}" >&2
  exit 1
fi

# Validate that headers reflect capped and full cluster membership counts.
cluster1_header=$(sed -n '5p' "${OUT_R1_FASTQ}")
cluster2_header=$(sed -n '9p' "${OUT_R1_FASTQ}")
cluster1_members=$(grep -o ';' <<< "${cluster1_header}" | wc -l | tr -d ' \t')
cluster2_members=$(grep -o ';' <<< "${cluster2_header}" | wc -l | tr -d ' \t')
if [ "${cluster1_members}" -ne 1000 ]; then
  echo "Cluster 1 header should list 1000 members, found ${cluster1_members}" >&2
  exit 1
fi
if [ "${cluster2_members}" -ne 10 ]; then
  echo "Cluster 2 header should list 10 members, found ${cluster2_members}" >&2
  exit 1
fi

# Display outputs for manual inspection.
truncate_fastq_headers() {
  local file=$1
  while IFS= read -r header && IFS= read -r seq && IFS= read -r plus && IFS= read -r qual; do
    if [ -z "${header}" ]; then
      continue
    fi
    local short_header
    if [ "${#header}" -gt 80 ]; then
      short_header="$(printf '%.80s' "${header}")..."
    else
      short_header="${header}"
    fi
    printf '%s\n%s\n%s\n%s\n' "${short_header}" "${seq}" "${plus}" "${qual}"
  done < "${file}"
}

echo "=== calib_cons R1 consensus FASTQ (headers truncated for display) ==="
truncate_fastq_headers "${OUT_R1_FASTQ}"
echo "=== calib_cons R2 consensus FASTQ (headers truncated for display) ==="
truncate_fastq_headers "${OUT_R2_FASTQ}"
