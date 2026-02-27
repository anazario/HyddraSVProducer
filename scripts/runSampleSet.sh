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
#   TRACKCOLLECTION    - track/vertex collection value
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

# ── Track collection menus per format ─────────────────────────────────────────
TRACK_COLLECTIONS_MINIAOD=(
    "sip2DMuonEnhanced"
    "sip2DMuonEnhancedWithEle"
    "sip2DSlimmedDisplacedMuonEnhancedWithEle"
    "promptMuonExtracted"
    "displacedMuonExtracted"
    "pf"
    "lost"
    "eleLost"
    "merged"
    "mergedWithEle"
    "mergedAll"
)

TRACK_COLLECTIONS_FASTSIM=(
    "fastSimSip2DMuonEnhanced"
    "fastSimMuonEnhanced"
    "fastSimSip2D"
    "fastSimSelected"
    "general"
    "generalFiltered"
)

TRACK_COLLECTIONS_AOD=(
    "sip2DMuonEnhanced"
    "promptMuonBestTrack"
    "displacedMuonBestTrack"
    "promptMuonExtracted"
    "displacedMuonExtracted"
    "promptMuonPriority"
    "displacedMuonPriority"
    "muonEnhanced"
    "sip2D"
    "general"
    "generalFiltered"
    "selected"
    "displacedGlobalMuon"
    "displacedStandAloneMuon"
)

VERTEX_COLLECTIONS_LLPNANO=(
    "PatMuonVertex"
    "PatDSAMuonVertex"
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

# Detect format type from config filename
if [[ "$CONFIG" == *"_miniAOD_"* ]]; then
    FORMAT_TYPE="miniAOD"
elif [[ "$CONFIG" == *"_fastSim_"* ]]; then
    FORMAT_TYPE="fastSim"
else
    FORMAT_TYPE="AOD"
fi

# Detect which optional features this config supports
HAS_TRACK_COLLECTION=false
HAS_VERTEX_COLLECTION=false
HAS_PROCESS_MODE=false
HAS_APPLY_CUTS=false

grep -q "options.register('trackCollection'" "$CONFIG_PATH" 2>/dev/null && HAS_TRACK_COLLECTION=true
grep -q "options.register('collection'"      "$CONFIG_PATH" 2>/dev/null && HAS_VERTEX_COLLECTION=true
grep -q "options.register('processMode'"     "$CONFIG_PATH" 2>/dev/null && HAS_PROCESS_MODE=true
grep -q "options.register('applyCuts'"       "$CONFIG_PATH" 2>/dev/null && HAS_APPLY_CUTS=true

echo ""
echo -e "  Detected format: ${YELLOW}${FORMAT_TYPE}${NC}"
$HAS_PROCESS_MODE && echo -e "  processMode:     ${GREEN}supported${NC}"
$HAS_APPLY_CUTS   && echo -e "  applyCuts:       ${GREEN}supported${NC}"
echo ""

# ── Step 3: Track / vertex collection ────────────────────────────────────────
echo -e "${BOLD}Step 3: Select track collection${NC}"
echo ""

TRACK_COLLECTION=""
COLLECTION_ARG=""   # the parallelRun.sh flag to use (--track-collection or --collection)

if $HAS_TRACK_COLLECTION; then
    COLLECTION_ARG="--track-collection"

    case "$FORMAT_TYPE" in
        miniAOD) MENU_COLLECTIONS=("${TRACK_COLLECTIONS_MINIAOD[@]}") ;;
        fastSim) MENU_COLLECTIONS=("${TRACK_COLLECTIONS_FASTSIM[@]}") ;;
        *)       MENU_COLLECTIONS=("${TRACK_COLLECTIONS_AOD[@]}") ;;
    esac

    for i in "${!MENU_COLLECTIONS[@]}"; do
        printf "  %2d) %s\n" "$((i+1))" "${MENU_COLLECTIONS[$i]}"
    done
    N_COLL="${#MENU_COLLECTIONS[@]}"
    printf "  %2d) %s\n" "$((N_COLL+1))" "custom (enter manually)"
    echo ""

    COLL_CHOICE=""
    while [[ -z "$COLL_CHOICE" ]]; do
        read -rp "$(echo -e "${CYAN}Track collection [1-$((N_COLL+1))]: ${NC}")" COLL_CHOICE
        if [[ "$COLL_CHOICE" =~ ^[0-9]+$ ]] && [[ "$COLL_CHOICE" -ge 1 ]] && [[ "$COLL_CHOICE" -le "$((N_COLL+1))" ]]; then
            if [[ "$COLL_CHOICE" -le "$N_COLL" ]]; then
                TRACK_COLLECTION="${MENU_COLLECTIONS[$((COLL_CHOICE-1))]}"
            else
                read -rp "$(echo -e "${CYAN}  Enter collection name: ${NC}")" TRACK_COLLECTION
                [[ -z "$TRACK_COLLECTION" ]] && { echo -e "${RED}  Cannot be empty.${NC}"; COLL_CHOICE=""; }
            fi
        else
            echo -e "${RED}  Please enter a number between 1 and $((N_COLL+1)).${NC}"
            COLL_CHOICE=""
        fi
    done

elif $HAS_VERTEX_COLLECTION; then
    COLLECTION_ARG="--collection"

    echo "  (This config uses a vertex collection rather than a track collection)"
    echo ""
    for i in "${!VERTEX_COLLECTIONS_LLPNANO[@]}"; do
        printf "  %2d) %s\n" "$((i+1))" "${VERTEX_COLLECTIONS_LLPNANO[$i]}"
    done
    echo ""

    COLL_CHOICE=""
    while [[ -z "$COLL_CHOICE" ]]; do
        read -rp "$(echo -e "${CYAN}Vertex collection [1-${#VERTEX_COLLECTIONS_LLPNANO[@]}]: ${NC}")" COLL_CHOICE
        if [[ "$COLL_CHOICE" =~ ^[0-9]+$ ]] \
            && [[ "$COLL_CHOICE" -ge 1 ]] \
            && [[ "$COLL_CHOICE" -le "${#VERTEX_COLLECTIONS_LLPNANO[@]}" ]]; then
            TRACK_COLLECTION="${VERTEX_COLLECTIONS_LLPNANO[$((COLL_CHOICE-1))]}"
        else
            echo -e "${RED}  Please enter 1 or 2.${NC}"
            COLL_CHOICE=""
        fi
    done

else
    echo "  (This config does not register a track/vertex collection — skipping.)"
    TRACK_COLLECTION="default"
fi

# ── Step 4: Process mode (menu, only if config supports it) ───────────────────
EXTRA_ARGS=""

if $HAS_PROCESS_MODE; then
    echo ""
    echo -e "${BOLD}Step 4: Process mode${NC}"
    echo ""
    echo "  1) both       (default — run leptonic and hadronic)"
    echo "  2) leptonic   (leptonic SVs only)"
    echo "  3) hadronic   (hadronic SVs only)"
    echo ""

    PM_CHOICE=""
    while [[ -z "$PM_CHOICE" ]]; do
        read -rp "$(echo -e "${CYAN}Process mode [1/2/3, default 1]: ${NC}")" PM_CHOICE
        PM_CHOICE="${PM_CHOICE:-1}"
        case "$PM_CHOICE" in
            1) ;;   # 'both' is the default; no flag needed
            2) EXTRA_ARGS="$EXTRA_ARGS --process-mode leptonic" ;;
            3) EXTRA_ARGS="$EXTRA_ARGS --process-mode hadronic" ;;
            *)
                echo -e "${RED}  Please enter 1, 2, or 3.${NC}"
                PM_CHOICE=""
                ;;
        esac
    done
fi

# ── Step 5: Apply cuts (only if config supports it) ───────────────────────────
if $HAS_APPLY_CUTS; then
    echo ""
    echo -e "${BOLD}Step 5: Track quality cuts${NC}"
    echo ""
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
fi

# ── Step 6: Output directory ──────────────────────────────────────────────────
echo ""
echo -e "${BOLD}Step 6: Output directory for merged .root files${NC}"
echo ""
read -rp "$(echo -e "${CYAN}Output directory [output]: ${NC}")" OUTPUT_DIR_INPUT
OUTPUT_DIR="${PROJECT_ROOT}/${OUTPUT_DIR_INPUT:-output}"

# ── Step 7: Optional flag ─────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}Step 7: Optional label appended to output filename (leave blank to omit)${NC}"
echo "  Examples: leptonic, hadronic, noCuts, applyCuts, noGen"
echo ""
read -rp "$(echo -e "${CYAN}Optional flag [none]: ${NC}")" OPTIONAL_FLAG

# ── Step 8: Extra parallelRun flags ──────────────────────────────────────────
echo ""
echo -e "${BOLD}Step 8: Additional parallelRun.sh flags${NC}"
echo ""
echo "  Available options (not yet specified above):"
echo "    --no-gen                  Disable gen info (for data)"
echo "    --collection X            Vertex collection (e.g. PatMuonVertex)"
echo "    --files-per-job N         Input files per job (default: 1)"
echo "    --continue                Skip files completed in a previous run"
echo "    --no-merge                Skip the hadd merge step"
if $HAS_APPLY_CUTS && [[ "$APPLY_CUTS" != "y" ]]; then
    echo "    --apply-cuts              Enable track quality cuts"
    echo "    --min-pt VALUE            Minimum track pT in GeV (default: 1.0)"
    echo "    --min-abs-sip2d VALUE     Minimum |sip2D| (default: 4.0)"
    echo "    --max-chi2 VALUE          Maximum normalized chi2 (default: 5.0)"
fi
echo ""
read -rp "$(echo -e "${CYAN}Extra flags [none]: ${NC}")" CUSTOM_EXTRA
[[ -n "$CUSTOM_EXTRA" ]] && EXTRA_ARGS="$EXTRA_ARGS $CUSTOM_EXTRA"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo -e "${GREEN}${BOLD}================================================${NC}"
echo -e "${GREEN}${BOLD}  Configuration Summary${NC}"
echo -e "${GREEN}${BOLD}================================================${NC}"
echo ""
printf "  %-22s %s\n" "Mass group:"        "$MASS_LABEL"
printf "  %-22s %s\n" "Samples:"           "${#SELECTED_FILES[@]}"
for f in "${SELECTED_FILES[@]}"; do printf "    - %s\n" "$f"; done
printf "  %-22s %s\n" "Config:"            "$CONFIG"
printf "  %-22s %s\n" "Format:"            "$FORMAT_TYPE"
printf "  %-22s %s\n" "Analyzer name:"     "$ANALYZER_NAME"
printf "  %-22s %s\n" "Track collection:"  "$TRACK_COLLECTION"
printf "  %-22s %s\n" "Output dir:"        "$OUTPUT_DIR"
printf "  %-22s %s\n" "Optional flag:"     "${OPTIONAL_FLAG:-<none>}"
printf "  %-22s %s\n" "Parallel jobs:"     "8"
printf "  %-22s %s\n" "Extra args:"        "${EXTRA_ARGS:-<none>}"
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

    # Run parallelRun.sh; capture exit code without propagating set -e
    set +e
    bash "$PARALLEL_SCRIPT" \
        "$FILES_DIR/$SAMPLE" \
        -c "$CONFIG_PATH" \
        $COLLECTION_ARG "$TRACK_COLLECTION" \
        -j 8 \
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
