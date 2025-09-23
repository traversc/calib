#!/usr/bin/env bash
set -xeuo pipefail

if ! command -v gzip >/dev/null 2>&1; then
  echo "gzip is required for this test"
  exit 1
fi
if ! command -v zstd >/dev/null 2>&1; then
  echo "zstd is required for this test"
  exit 1
fi

TMP_DIR="$(mktemp -d 2>/dev/null || mktemp -d -t calib-test)"
cleanup() {
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

R1_PLAIN="${TMP_DIR}/reads_R1.fastq"
R2_PLAIN="${TMP_DIR}/reads_R2.fastq"
cat <<'DATA' > "${R1_PLAIN}"
@read1/1
AACCGGTT
+
FFFFFFFF
@read2/1
AACCGGTA
+
FFFFFFFF
DATA
cat <<'DATA' > "${R2_PLAIN}"
@read1/2
TTGGCCAA
+
FFFFFFFF
@read2/2
TTGGCCAT
+
FFFFFFFF
DATA

gzip -c "${R1_PLAIN}" > "${TMP_DIR}/reads_R1.fastq.gz"
gzip -c "${R2_PLAIN}" > "${TMP_DIR}/reads_R2.fastq.gz"
zstd --quiet -T1 -f -o "${TMP_DIR}/reads_R1.fastq.zst" "${R1_PLAIN}"
zstd --quiet -T1 -f -o "${TMP_DIR}/reads_R2.fastq.zst" "${R2_PLAIN}"

run_calib() {
  local name="$1"; shift
  local in1="$1"; local in2="$2"; shift 2
  local outcomp="$1"; shift
  ./calib \
    --input-forward "${in1}" \
    --input-reverse "${in2}" \
    --output-prefix "${TMP_DIR}/${name}." \
    --output-format "${outcomp}" \
    --barcode-length 1 \
    --ignored-sequence-prefix-length 0 \
    --minimizer-count 1 \
    --kmer-size 4 \
    --error-tolerance 1 \
    --minimizer-threshold 1 \
    --threads 1 \
    --silent "$@"
}

run_calib plain_plain "${R1_PLAIN}" "${R2_PLAIN}" plain
run_calib plain_force "${R1_PLAIN}" "${R2_PLAIN}" plain --input-format plain
run_calib gz_plain "${TMP_DIR}/reads_R1.fastq.gz" "${TMP_DIR}/reads_R2.fastq.gz" plain --input-format gzip
run_calib zst_plain "${TMP_DIR}/reads_R1.fastq.zst" "${TMP_DIR}/reads_R2.fastq.zst" plain --input-format zstd
cmp -s "${TMP_DIR}/plain_plain.cluster" "${TMP_DIR}/plain_force.cluster"
cmp -s "${TMP_DIR}/plain_plain.cluster" "${TMP_DIR}/gz_plain.cluster"
cmp -s "${TMP_DIR}/plain_plain.cluster" "${TMP_DIR}/zst_plain.cluster"
run_calib plain_gzip "${R1_PLAIN}" "${R2_PLAIN}" gzip
run_calib plain_zstd "${R1_PLAIN}" "${R2_PLAIN}" zstd
gzip -dc "${TMP_DIR}/plain_gzip.cluster.gz" > "${TMP_DIR}/plain_gzip.cluster.txt"
zstd --quiet -T1 -dc "${TMP_DIR}/plain_zstd.cluster.zst" > "${TMP_DIR}/plain_zstd.cluster.txt"
cmp -s "${TMP_DIR}/plain_plain.cluster" "${TMP_DIR}/plain_gzip.cluster.txt"
cmp -s "${TMP_DIR}/plain_plain.cluster" "${TMP_DIR}/plain_zstd.cluster.txt"

run_cons() {
  local name="$1"; shift
  local in1="$1"; local in2="$2"; shift 2
  local cluster_file="$1"; shift
  local cluster_comp="$1"; shift
  local out_comp="$1"; shift
  ./consensus/calib_cons \
    -c "${cluster_file}" \
    --cluster-input-format "${cluster_comp}" \
    --output-format "${out_comp}" \
    -q "${in1}" "${in2}" \
    -o "${TMP_DIR}/${name}.R1." "${TMP_DIR}/${name}.R2." \
    -t 1 "$@"
}

run_cons plain_plain "${R1_PLAIN}" "${R2_PLAIN}" "${TMP_DIR}/plain_plain.cluster" auto plain
run_cons plain_force "${R1_PLAIN}" "${R2_PLAIN}" "${TMP_DIR}/plain_plain.cluster" plain plain --input-format plain
run_cons gz_plain_in "${TMP_DIR}/reads_R1.fastq.gz" "${TMP_DIR}/reads_R2.fastq.gz" "${TMP_DIR}/plain_plain.cluster" auto plain --input-format gzip
run_cons zst_plain_in "${TMP_DIR}/reads_R1.fastq.zst" "${TMP_DIR}/reads_R2.fastq.zst" "${TMP_DIR}/plain_plain.cluster" auto plain --input-format zstd
cmp -s "${TMP_DIR}/plain_plain.R1.fastq" "${TMP_DIR}/plain_force.R1.fastq"
cmp -s "${TMP_DIR}/plain_plain.R2.fastq" "${TMP_DIR}/plain_force.R2.fastq"
cmp -s "${TMP_DIR}/plain_plain.R1.fastq" "${TMP_DIR}/gz_plain_in.R1.fastq"
cmp -s "${TMP_DIR}/plain_plain.R2.fastq" "${TMP_DIR}/gz_plain_in.R2.fastq"
cmp -s "${TMP_DIR}/plain_plain.R1.fastq" "${TMP_DIR}/zst_plain_in.R1.fastq"
cmp -s "${TMP_DIR}/plain_plain.R2.fastq" "${TMP_DIR}/zst_plain_in.R2.fastq"
run_cons plain_gzip_out "${R1_PLAIN}" "${R2_PLAIN}" "${TMP_DIR}/plain_plain.cluster" plain gzip
run_cons plain_zstd_out "${R1_PLAIN}" "${R2_PLAIN}" "${TMP_DIR}/plain_plain.cluster" plain zstd
gzip -dc "${TMP_DIR}/plain_gzip_out.R1.fastq.gz" > "${TMP_DIR}/plain_gzip_out.R1.fastq.txt"
gzip -dc "${TMP_DIR}/plain_gzip_out.R2.fastq.gz" > "${TMP_DIR}/plain_gzip_out.R2.fastq.txt"
zstd --quiet -T1 -dc "${TMP_DIR}/plain_zstd_out.R1.fastq.zst" > "${TMP_DIR}/plain_zstd_out.R1.fastq.txt"
zstd --quiet -T1 -dc "${TMP_DIR}/plain_zstd_out.R2.fastq.zst" > "${TMP_DIR}/plain_zstd_out.R2.fastq.txt"
cmp -s "${TMP_DIR}/plain_plain.R1.fastq" "${TMP_DIR}/plain_gzip_out.R1.fastq.txt"
cmp -s "${TMP_DIR}/plain_plain.R2.fastq" "${TMP_DIR}/plain_gzip_out.R2.fastq.txt"
cmp -s "${TMP_DIR}/plain_plain.R1.fastq" "${TMP_DIR}/plain_zstd_out.R1.fastq.txt"
cmp -s "${TMP_DIR}/plain_plain.R2.fastq" "${TMP_DIR}/plain_zstd_out.R2.fastq.txt"

echo "All compression format checks passed"
