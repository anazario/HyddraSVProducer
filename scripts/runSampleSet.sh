#!/bin/bash
#
# runSampleSet.sh - Interactively run parallelRun.sh over a set of 5 TTALP samples
#
# Processes all 5 ctau samples for a chosen mass group (2 GeV or 60 GeV) sequentially.
# Between each sample the temporary parallel_output directory is removed to avoid
# mixing files across samples. Merged outputs are saved with the naming convention:
#
#   {FILENAME}_{ANALYZER}_{TRACKCOLLECTION}_{OPTIONALFLAG}.root
#
#   FILENAME           - sample .txt name with .txt removed
#   ANALYZER           - config file with leading 'test' and trailing '_cfg.py' removed
#   TRACKCOLLECTION    - value passed via --track-collection
#   OPTIONALFLAG       - optional label (omitted from name if left blank)
#
# Usage:
#   bash scripts/runSampleSet.sh
#

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
FILES_DIR="$PROJECT_ROOT/files"
TEST_DIR="$PROJECT_ROOT/test"
PARALLEL_SCRIPT="$SCRIPT_DIR/parallelRun.sh"
TEMP_OUTPUT="$PROJECT_ROOT/parallel_output"

# ── Colors ────────────────────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

# ── Sample file groups ────────────────────────────────────────────────────────
FILES_MASS2=(
    "TTALPto2Mu_MALP-2_ctau-1e-5mm_AODSIM.txt"
    "TTALPto2Mu_MALP-2_ctau-1e0mm_AODSIM.txt"
    "TTALPto2Mu_MALP-2_ctau-1e1mm_AODSIM.txt"
    "TTALPto2Mu_MALP-2_ctau-1e2mm_AODSIM.txt"
    "TTALPto2Mu_MALP-2_ctau-1e3mm_AODSIM.txt"
)

FILES_MASS60=(
    "TTALPto2Mu_MALP-60_ctau-1e-5mm_AODSIM.txt"
    "TTALPto2Mu_MALP-60_ctau-1e0mm_AODSIM.txt"
    "TTALPto2Mu_MALP-60_ctau-1e1mm_AODSIM.txt"
    "TTALPto2Mu_MALP-60_ctau-1e2mm_AODSIM.txt"
    "TTALPto2Mu_MALP-60_ctau-1e3mm_AODSIM.txt"
)

# ── Helpers ───────────────────────────────────────────────────────────────────
print_header() {
    echo ""
    echo -e "${BLUE}${BOLD}================================================${NC}"
    echo -e "${BLUE}${BOLD}       TTALP Sample Set Processor${NC}"
    echo -e "${BLUE}${BOLD}================================================${NC}"
    echo ""
}

# prompt_yn <question> <default y|n>  →  echoes "y" or "n"
prompt_yn() {
    local question="$1"
    local default="${2:-n}"
    local answer
    while true; do
        if [[ "$default" == "y" ]]; then
            read -rp "$(echo -e "${CYAN}${question} [Y/n]: ${NC}")" answer
            answer="${answer:-y}"
        else
            read -rp "$(echo -e "${CYAN}${question} [y/N]: ${NC}")" answer
            answer="${answer:-n}"
        fi
        case "$answer" in
            [yY]*) echo "y"; return ;;
            [nN]*) echo "n"; return ;;
            *)     echo -e "${RED}  Please enter y or n.${NC}" >&2 ;;
        esac
    done
}

# ── Sanity checks ─────────────────────────────────────────────────────────────
if [[ ! -f "$PARALLEL_SCRIPT" ]]; then
    echo -e "${RED}Error: parallelRun.sh not found at $PARALLEL_SCRIPT${NC}"
    exit 1
fi

if [[ ! -d "$FILES_DIR" ]]; then
    echo -e "${RED}Error: files/ directory not found at $FILES_DIR${NC}"
    exit 1
fi

# ═════════════════════════════════════════════════════════════════════════════
# Interactive prompts
# ═════════════════════════════════════════════════════════════════════════════
print_header

# ── Step 1: Mass group ────────────────────────────────────────────────────────
echo -e "${BOLD}Step 1: Select mass group${NC}"
echo "  1) MALP-2   (2 GeV  — 5 samples)"
echo "  2) MALP-60  (60 GeV — 5 samples)"
echo "  3) Both     (10 samples total)"
echo ""

MASS_CHOICE=""
while [[ -z "$MASS_CHOICE" ]]; do
    read -rp "$(echo -e "${CYAN}Mass group [1/2/3]: ${NC}")" MASS_CHOICE
    case "$MASS_CHOICE" in
        1)
            SELECTED_FILES=("${FILES_MASS2[@]}")
            MASS_LABEL="MALP-2 (2 GeV)"
            ;;
        2)
            SELECTED_FILES=("${FILES_MASS60[@]}")
            MASS_LABEL="MALP-60 (60 GeV)"
            ;;
        3)
            SELECTED_FILES=("${FILES_MASS2[@]}" "${FILES_MASS60[@]}")
            MASS_LABEL="Both (2 GeV + 60 GeV)"
            ;;
        *)
            echo -e "${RED}  Please enter 1, 2, or 3.${NC}"
            MASS_CHOICE=""
            ;;
    esac
done

# ── Step 2: Config file ───────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}Step 2: Select config file${NC}"

mapfile -t CONFIGS < <(find "$TEST_DIR" -maxdepth 1 -name "test*_cfg.py" \
                        -exec basename {} \; 2>/dev/null | sort)

if [[ ${#CONFIGS[@]} -eq 0 ]]; then
    echo -e "${RED}No config files found in $TEST_DIR${NC}"
    exit 1
fi

for i in "${!CONFIGS[@]}"; do
    printf "  %2d) %s\n" "$((i+1))" "${CONFIGS[$i]}"
done
echo ""

CONFIG_CHOICE=""
while [[ -z "$CONFIG_CHOICE" ]]; do
    read -rp "$(echo -e "${CYAN}Config [1-${#CONFIGS[@]}]: ${NC}")" CONFIG_CHOICE
    if [[ "$CONFIG_CHOICE" =~ ^[0-9]+$ ]] \
        && [[ "$CONFIG_CHOICE" -ge 1 ]] \
        && [[ "$CONFIG_CHOICE" -le "${#CONFIGS[@]}" ]]; then
        CONFIG="${CONFIGS[$((CONFIG_CHOICE-1))]}"
        CONFIG_PATH="$TEST_DIR/$CONFIG"
    else
        echo -e "${RED}  Please enter a number between 1 and ${#CONFIGS[@]}.${NC}"
        CONFIG_CHOICE=""
    fi
done

# Derive analyzer name: strip leading 'test' and trailing '_cfg.py'
CONFIG_NO_CFG="${CONFIG%_cfg.py}"
ANALYZER_NAME="${CONFIG_NO_CFG#test}"

# ── Step 3: Track collection ──────────────────────────────────────────────────
echo ""
echo -e "${BOLD}Step 3: Track collection${NC}"
echo "  Common values: general, generalWithDSA, dsa"
echo ""

TRACK_COLLECTION=""
while [[ -z "$TRACK_COLLECTION" ]]; do
    read -rp "$(echo -e "${CYAN}Track collection: ${NC}")" TRACK_COLLECTION
    [[ -z "$TRACK_COLLECTION" ]] && echo -e "${RED}  Track collection cannot be empty.${NC}"
done

# ── Step 4: Output directory ──────────────────────────────────────────────────
echo ""
echo -e "${BOLD}Step 4: Output directory for merged .root files${NC}"
echo ""
read -rp "$(echo -e "${CYAN}Output directory [output]: ${NC}")" OUTPUT_DIR
OUTPUT_DIR="${PROJECT_ROOT}/${OUTPUT_DIR:-output}"

# ── Step 5: Optional flag ─────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}Step 5: Optional label appended to output filename (leave blank to omit)${NC}"
echo "  Examples: leptonic, hadronic, noCuts, applyCuts, noGen"
echo ""
read -rp "$(echo -e "${CYAN}Optional flag [none]: ${NC}")" OPTIONAL_FLAG

# ── Step 6: Parallel jobs ─────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}Step 6: parallelRun options${NC}"
echo ""
read -rp "$(echo -e "${CYAN}Number of parallel jobs [8]: ${NC}")" N_JOBS
N_JOBS="${N_JOBS:-8}"

# ── Step 7: Additional options ────────────────────────────────────────────────
echo ""
echo -e "${BOLD}Step 7: Additional options${NC}"
echo ""

EXTRA_ARGS=""

# Process mode
echo "  Process mode: both (default), leptonic, hadronic"
read -rp "$(echo -e "${CYAN}  Process mode [both]: ${NC}")" PROCESS_MODE
if [[ -n "$PROCESS_MODE" && "$PROCESS_MODE" != "both" ]]; then
    EXTRA_ARGS="$EXTRA_ARGS --process-mode $PROCESS_MODE"
fi

echo ""

# No-gen
NO_GEN=$(prompt_yn "  Disable gen info (--no-gen)?" "n")
[[ "$NO_GEN" == "y" ]] && EXTRA_ARGS="$EXTRA_ARGS --no-gen"

# Apply cuts
APPLY_CUTS=$(prompt_yn "  Enable track quality cuts (--apply-cuts)?" "n")
if [[ "$APPLY_CUTS" == "y" ]]; then
    EXTRA_ARGS="$EXTRA_ARGS --apply-cuts"
    read -rp "$(echo -e "${CYAN}    Min pT in GeV [1.0]: ${NC}")" MIN_PT
    [[ -n "$MIN_PT" ]] && EXTRA_ARGS="$EXTRA_ARGS --min-pt $MIN_PT"
    read -rp "$(echo -e "${CYAN}    Min |sip2D| [4.0]: ${NC}")" MIN_SIP2D
    [[ -n "$MIN_SIP2D" ]] && EXTRA_ARGS="$EXTRA_ARGS --min-abs-sip2d $MIN_SIP2D"
    read -rp "$(echo -e "${CYAN}    Max normalized chi2 [5.0]: ${NC}")" MAX_CHI2
    [[ -n "$MAX_CHI2" ]] && EXTRA_ARGS="$EXTRA_ARGS --max-chi2 $MAX_CHI2"
fi

echo ""

# Custom extra args (pass-through to parallelRun.sh)
read -rp "$(echo -e "${CYAN}  Any additional parallelRun.sh flags [none]: ${NC}")" CUSTOM_EXTRA
[[ -n "$CUSTOM_EXTRA" ]] && EXTRA_ARGS="$EXTRA_ARGS $CUSTOM_EXTRA"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo -e "${GREEN}${BOLD}================================================${NC}"
echo -e "${GREEN}${BOLD}  Configuration Summary${NC}"
echo -e "${GREEN}${BOLD}================================================${NC}"
echo ""
printf "  %-20s %s\n" "Mass group:"      "$MASS_LABEL"
printf "  %-20s %s\n" "Samples:"         "${#SELECTED_FILES[@]}"
for f in "${SELECTED_FILES[@]}"; do printf "    - %s\n" "$f"; done
printf "  %-20s %s\n" "Config:"          "$CONFIG"
printf "  %-20s %s\n" "Analyzer name:"   "$ANALYZER_NAME"
printf "  %-20s %s\n" "Track collection:" "$TRACK_COLLECTION"
printf "  %-20s %s\n" "Output dir:"      "$OUTPUT_DIR"
printf "  %-20s %s\n" "Optional flag:"   "${OPTIONAL_FLAG:-<none>}"
printf "  %-20s %s\n" "Parallel jobs:"   "$N_JOBS"
printf "  %-20s %s\n" "Extra args:"      "${EXTRA_ARGS:-<none>}"
echo ""
echo "  Output naming example:"
EXAMPLE_BASE="${SELECTED_FILES[0]%.txt}"
if [[ -n "$OPTIONAL_FLAG" ]]; then
    EXAMPLE_OUT="${EXAMPLE_BASE}_${ANALYZER_NAME}_${TRACK_COLLECTION}_${OPTIONAL_FLAG}.root"
else
    EXAMPLE_OUT="${EXAMPLE_BASE}_${ANALYZER_NAME}_${TRACK_COLLECTION}.root"
fi
echo "    $EXAMPLE_OUT"
echo ""

CONFIRM=$(prompt_yn "Proceed with these settings?" "y")
if [[ "$CONFIRM" != "y" ]]; then
    echo -e "${YELLOW}Aborted.${NC}"
    exit 0
fi

# ═════════════════════════════════════════════════════════════════════════════
# Process each sample
# ═════════════════════════════════════════════════════════════════════════════
mkdir -p "$OUTPUT_DIR"

N_TOTAL="${#SELECTED_FILES[@]}"
N_SUCCESS=0
N_FAILED=0
FAILED_SAMPLES=()

for i in "${!SELECTED_FILES[@]}"; do
    SAMPLE="${SELECTED_FILES[$i]}"
    SAMPLE_NUM=$((i+1))

    FILENAME_NOEXT="${SAMPLE%.txt}"
    if [[ -n "$OPTIONAL_FLAG" ]]; then
        OUTPUT_NAME="${FILENAME_NOEXT}_${ANALYZER_NAME}_${TRACK_COLLECTION}_${OPTIONAL_FLAG}.root"
    else
        OUTPUT_NAME="${FILENAME_NOEXT}_${ANALYZER_NAME}_${TRACK_COLLECTION}.root"
    fi

    echo ""
    echo -e "${BLUE}${BOLD}------------------------------------------------${NC}"
    echo -e "${BLUE}${BOLD}  Sample $SAMPLE_NUM / $N_TOTAL${NC}"
    echo -e "${BLUE}${BOLD}  $SAMPLE${NC}"
    echo -e "${BLUE}${BOLD}------------------------------------------------${NC}"
    echo "  Output: $OUTPUT_DIR/$OUTPUT_NAME"
    echo ""

    # Verify the input file exists
    if [[ ! -f "$FILES_DIR/$SAMPLE" ]]; then
        echo -e "${RED}  ERROR: Input file not found: $FILES_DIR/$SAMPLE${NC}"
        echo -e "${RED}  Skipping.${NC}"
        N_FAILED=$((N_FAILED+1))
        FAILED_SAMPLES+=("$SAMPLE")
        continue
    fi

    # Clean up any leftover temp output from a previous run
    [[ -d "$TEMP_OUTPUT" ]] && rm -rf "$TEMP_OUTPUT"

    # Run parallelRun.sh; capture exit code without letting set -e propagate
    set +e
    bash "$PARALLEL_SCRIPT" \
        "$FILES_DIR/$SAMPLE" \
        -c "$CONFIG_PATH" \
        --track-collection "$TRACK_COLLECTION" \
        -j "$N_JOBS" \
        -o "$TEMP_OUTPUT" \
        -n "$OUTPUT_NAME" \
        $EXTRA_ARGS
    PARALLEL_EXIT=$?
    set -e

    if [[ $PARALLEL_EXIT -ne 0 ]]; then
        echo -e "${RED}  parallelRun.sh exited with status $PARALLEL_EXIT for $SAMPLE${NC}"
        N_FAILED=$((N_FAILED+1))
        FAILED_SAMPLES+=("$SAMPLE")
    elif [[ -f "$TEMP_OUTPUT/$OUTPUT_NAME" ]]; then
        mv "$TEMP_OUTPUT/$OUTPUT_NAME" "$OUTPUT_DIR/$OUTPUT_NAME"
        echo -e "${GREEN}  Saved: $OUTPUT_DIR/$OUTPUT_NAME${NC}"
        N_SUCCESS=$((N_SUCCESS+1))
    else
        echo -e "${RED}  WARNING: parallelRun.sh succeeded but merged output not found: $TEMP_OUTPUT/$OUTPUT_NAME${NC}"
        N_FAILED=$((N_FAILED+1))
        FAILED_SAMPLES+=("$SAMPLE")
    fi

    # Always discard the temp directory between samples
    if [[ -d "$TEMP_OUTPUT" ]]; then
        echo "  Removing $TEMP_OUTPUT..."
        rm -rf "$TEMP_OUTPUT"
    fi
done

# ═════════════════════════════════════════════════════════════════════════════
# Final summary
# ═════════════════════════════════════════════════════════════════════════════
echo ""
echo -e "${GREEN}${BOLD}================================================${NC}"
echo -e "${GREEN}${BOLD}  Final Summary${NC}"
echo -e "${GREEN}${BOLD}================================================${NC}"
echo ""
printf "  %-20s %s\n" "Total samples:"   "$N_TOTAL"
echo -e "  ${GREEN}$(printf "%-20s" "Successful:")      $N_SUCCESS${NC}"
if [[ $N_FAILED -gt 0 ]]; then
    echo -e "  ${RED}$(printf "%-20s" "Failed:")          $N_FAILED${NC}"
    for f in "${FAILED_SAMPLES[@]}"; do
        echo -e "    ${RED}- $f${NC}"
    done
fi
echo ""
echo "  Outputs in: $OUTPUT_DIR/"
if ls "$OUTPUT_DIR"/*.root &>/dev/null; then
    ls -lh "$OUTPUT_DIR"/*.root | awk '{printf "    %-60s %s\n", $NF, $5}'
else
    echo "    (no .root files found)"
fi
echo ""
echo -e "${GREEN}${BOLD}Done!${NC}"
