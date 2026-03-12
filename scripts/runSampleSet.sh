#!/bin/bash
#
# runSampleSet.sh - Interactively run parallelRun.sh over a set of 5 TTALP samples
#
# Processes all 5 ctau samples for a chosen mass group (2 GeV or 60 GeV) sequentially.
# State is saved in the temporary workspace so a run can be resumed after interruption.
# Merged outputs are saved with the naming convention:
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
#   bash scripts/runSampleSet.sh --continue WORKSPACEDIR   # Resume using a specific workspace
#   bash scripts/runSampleSet.sh --continue                # Resume, searching for workspaces
#

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
FILES_DIR="$PROJECT_ROOT/files"
TEST_DIR="$PROJECT_ROOT/test"
PARALLEL_SCRIPT="$SCRIPT_DIR/parallelRun.sh"
TEMP_OUTPUT=""  # set interactively below

# ── Continue mode flag ────────────────────────────────────────────────────────
IS_CONTINUE=false
CONTINUE_WORKSPACE=""
if [[ "$1" == "--continue" ]]; then
    IS_CONTINUE=true
    if [[ -n "$2" ]]; then
        # Resolve workspace path (relative → under project root)
        if [[ "$2" == /* ]]; then
            CONTINUE_WORKSPACE="$2"
        else
            CONTINUE_WORKSPACE="$PROJECT_ROOT/$2"
        fi
    fi
fi

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
    "promptMuonLLPNano"
    "displacedMuonLLPNano"
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

# save_state — write (or overwrite) $TEMP_OUTPUT/.runstate with all run
# settings and the current COMPLETED_SAMPLES progress list.
save_state() {
    mkdir -p "$TEMP_OUTPUT"
    {
        echo "# runSampleSet state — updated $(date)"
        declare -p SELECTED_FILES
        declare -p COMPLETED_SAMPLES
        for _sv in CONFIG_PATH CONFIG ANALYZER_NAME FORMAT_TYPE \
                   TRACK_COLLECTION COLLECTION_ARG EXTRA_ARGS OPTIONAL_FLAG \
                   TEMP_OUTPUT OUTPUT_DIR MASS_LABEL HYDDRA_PRESET; do
            declare -p "$_sv" 2>/dev/null || true
        done
    } > "$TEMP_OUTPUT/.runstate"
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
# Interactive prompts  -OR-  state loading (continue mode)
# ═════════════════════════════════════════════════════════════════════════════
print_header

if $IS_CONTINUE; then
    # ── Continue mode: load saved state from workspace ────────────────────────
    echo -e "${BOLD}Continue mode — loading saved run state...${NC}"
    echo ""

    if [[ -n "$CONTINUE_WORKSPACE" ]]; then
        # Workspace specified directly on the command line
        CHOSEN_STATE="$CONTINUE_WORKSPACE/.runstate"
        if [[ ! -f "$CHOSEN_STATE" ]]; then
            echo -e "${RED}No state file found at $CHOSEN_STATE${NC}"
            echo "  Check the workspace path or start a fresh run."
            exit 1
        fi
    else
        # No workspace given — search for state files under the project root
        mapfile -t STATE_FILES < <(find "$PROJECT_ROOT" -maxdepth 2 -name ".runstate" 2>/dev/null | sort)

        if [[ ${#STATE_FILES[@]} -eq 0 ]]; then
            echo -e "${RED}No saved state files (.runstate) found under $PROJECT_ROOT${NC}"
            echo "  Start a fresh run (without --continue) first, or specify a workspace:"
            echo "    bash scripts/runSampleSet.sh --continue WORKSPACEDIR"
            exit 1
        fi

        if [[ ${#STATE_FILES[@]} -eq 1 ]]; then
            CHOSEN_STATE="${STATE_FILES[0]}"
            echo "  Auto-selected: $CHOSEN_STATE"
        else
            echo "  Found ${#STATE_FILES[@]} saved runs:"
            echo ""
            for i in "${!STATE_FILES[@]}"; do
                printf "  %2d) %s\n" "$((i+1))" "${STATE_FILES[$i]}"
                _od=$(grep '^declare -- OUTPUT_DIR=' "${STATE_FILES[$i]}" 2>/dev/null \
                      | sed 's/declare -- OUTPUT_DIR="\(.*\)"/\1/')
                _cfg=$(grep '^declare -- CONFIG=' "${STATE_FILES[$i]}" 2>/dev/null \
                       | sed 's/declare -- CONFIG="\(.*\)"/\1/')
                [[ -n "$_od" ]]  && printf "       Output dir: %s\n" "$_od"
                [[ -n "$_cfg" ]] && printf "       Config:     %s\n" "$_cfg"
            done
            echo ""

            STATE_CHOICE=""
            while [[ -z "$STATE_CHOICE" ]]; do
                read -rp "$(echo -e "${CYAN}Select run to continue [1-${#STATE_FILES[@]}]: ${NC}")" STATE_CHOICE
                if [[ "$STATE_CHOICE" =~ ^[0-9]+$ ]] \
                    && [[ "$STATE_CHOICE" -ge 1 ]] \
                    && [[ "$STATE_CHOICE" -le "${#STATE_FILES[@]}" ]]; then
                    CHOSEN_STATE="${STATE_FILES[$((STATE_CHOICE-1))]}"
                else
                    echo -e "${RED}  Please enter a number between 1 and ${#STATE_FILES[@]}.${NC}"
                    STATE_CHOICE=""
                fi
            done
        fi
    fi

    # Load state; provide empty defaults for arrays in case of old state files
    COMPLETED_SAMPLES=()
    # shellcheck source=/dev/null
    source "$CHOSEN_STATE"

    # If a workspace was specified, ensure TEMP_OUTPUT points to it
    # (guards against path changes between machines/sessions)
    [[ -n "$CONTINUE_WORKSPACE" ]] && TEMP_OUTPUT="$CONTINUE_WORKSPACE"

    echo -e "  ${GREEN}Loaded: $CHOSEN_STATE${NC}"
    echo ""

else
    # ── Fresh run: interactive prompts ────────────────────────────────────────
    COMPLETED_SAMPLES=()

    # ── Step 1: Mass group ────────────────────────────────────────────────────
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

    # ── Step 2: Config file ───────────────────────────────────────────────────
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
    HAS_HYDDRA_PRESET=false

    grep -q "options.register('trackCollection'" "$CONFIG_PATH" 2>/dev/null && HAS_TRACK_COLLECTION=true
    grep -q "options.register('collection'"      "$CONFIG_PATH" 2>/dev/null && HAS_VERTEX_COLLECTION=true
    grep -q "options.register('processMode'"     "$CONFIG_PATH" 2>/dev/null && HAS_PROCESS_MODE=true
    grep -q "options.register('applyCuts'"       "$CONFIG_PATH" 2>/dev/null && HAS_APPLY_CUTS=true
    grep -q "options.register('hyddraPreset'"    "$CONFIG_PATH" 2>/dev/null && HAS_HYDDRA_PRESET=true

    echo ""
    echo -e "  Detected format: ${YELLOW}${FORMAT_TYPE}${NC}"
    $HAS_PROCESS_MODE && echo -e "  processMode:     ${GREEN}supported${NC}"
    $HAS_APPLY_CUTS   && echo -e "  applyCuts:       ${GREEN}supported${NC}"
    echo ""

    # ── Step 3: Track / vertex collection ────────────────────────────────────
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

    # ── Step 4: Process mode (menu, only if config supports it) ───────────────
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

    # ── Step 5: Apply cuts (only if config supports it) ───────────────────────
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

    # ── Step 6: Output directory ──────────────────────────────────────────────
    echo ""
    echo -e "${BOLD}Step 6: Output directory for merged .root files${NC}"
    echo ""
    read -rp "$(echo -e "${CYAN}Output directory [output]: ${NC}")" OUTPUT_DIR_INPUT
    OUTPUT_DIR="${PROJECT_ROOT}/${OUTPUT_DIR_INPUT:-output}"

    # ── Step 7: Optional flag ─────────────────────────────────────────────────
    echo ""
    echo -e "${BOLD}Step 7: Optional label appended to output filename (leave blank to omit)${NC}"
    echo "  Examples: leptonic, hadronic, noCuts, applyCuts, noGen"
    echo ""
    read -rp "$(echo -e "${CYAN}Optional flag [none]: ${NC}")" OPTIONAL_FLAG

    # ── Step 8: Extra parallelRun flags ──────────────────────────────────────
    echo ""
    echo -e "${BOLD}Step 8: Additional parallelRun.sh flags${NC}"
    echo ""
    echo "  Available options (not yet specified above):"
    echo "    --no-gen                  Disable gen info (for data)"
    echo "    --collection X            Vertex collection (e.g. PatMuonVertex)"
    echo "    --files-per-job N         Input files per job (default: 1)"
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

    # ── Step 9: HYDDRA leptonic preset ───────────────────────────────────────
    HYDDRA_PRESET="default"
    if $HAS_HYDDRA_PRESET; then
        echo ""
        echo -e "${BOLD}Step 9: HYDDRA leptonic algorithm preset${NC}"
        echo ""
        echo "  1) default   (use config as-is)"
        echo "  2) NonIso    (merging OFF, cleaning OFF, smoothing ON, filter: charge neutrality only)"
        echo "  3) TightIso  (merging ON,  cleaning OFF, smoothing ON, filter: charge neutrality only)"
        echo ""

        HYDDRA_PRESET_CHOICE=""
        while [[ -z "$HYDDRA_PRESET_CHOICE" ]]; do
            read -rp "$(echo -e "${CYAN}Preset [1/2/3, default 1]: ${NC}")" HYDDRA_PRESET_CHOICE
            HYDDRA_PRESET_CHOICE="${HYDDRA_PRESET_CHOICE:-1}"
            case "$HYDDRA_PRESET_CHOICE" in
                1) HYDDRA_PRESET="default" ;;
                2) HYDDRA_PRESET="NonIso"
                   EXTRA_ARGS="$EXTRA_ARGS --hyddra-preset NonIso" ;;
                3) HYDDRA_PRESET="TightIso"
                   EXTRA_ARGS="$EXTRA_ARGS --hyddra-preset TightIso" ;;
                *)
                    echo -e "${RED}  Please enter 1, 2, or 3.${NC}"
                    HYDDRA_PRESET_CHOICE=""
                    ;;
            esac
        done
    fi

    # ── Step 10: Temporary working directory ─────────────────────────────────
    echo ""
    echo -e "${BOLD}Step 10: Temporary working directory name${NC}"
    echo "  (persists for the full run; state is saved here for --continue)"
    echo ""
    read -rp "$(echo -e "${CYAN}Temp dir name [parallel_output]: ${NC}")" TEMP_DIR_INPUT
    TEMP_OUTPUT="$PROJECT_ROOT/${TEMP_DIR_INPUT:-parallel_output}"

fi  # end if $IS_CONTINUE / else

# ═════════════════════════════════════════════════════════════════════════════
# Summary (shared between fresh and continue modes)
# ═════════════════════════════════════════════════════════════════════════════
N_TOTAL="${#SELECTED_FILES[@]}"

echo ""
echo -e "${GREEN}${BOLD}================================================${NC}"
echo -e "${GREEN}${BOLD}  Configuration Summary${NC}"
echo -e "${GREEN}${BOLD}================================================${NC}"
echo ""

if $IS_CONTINUE; then
    echo -e "  ${YELLOW}${BOLD}Mode: CONTINUE (completed samples will be skipped)${NC}"
    echo ""
fi

printf "  %-22s %s\n" "Mass group:"        "${MASS_LABEL:-<from state>}"
printf "  %-22s %s\n" "Samples:"           "$N_TOTAL"

for f in "${SELECTED_FILES[@]}"; do
    # Check against COMPLETED_SAMPLES list from state
    _done=false
    for _cs in "${COMPLETED_SAMPLES[@]}"; do
        [[ "$_cs" == "$f" ]] && { _done=true; break; }
    done
    if $_done; then
        printf "    ${GREEN}[done]${NC} %s\n" "$f"
    else
        printf "    - %s\n" "$f"
    fi
done

printf "  %-22s %s\n" "Config:"            "${CONFIG:-<from state>}"
printf "  %-22s %s\n" "Format:"            "${FORMAT_TYPE:-<from state>}"
printf "  %-22s %s\n" "Analyzer name:"     "$ANALYZER_NAME"
printf "  %-22s %s\n" "Track collection:"  "$TRACK_COLLECTION"
printf "  %-22s %s\n" "HYDDRA preset:"     "${HYDDRA_PRESET:-<from state>}"
printf "  %-22s %s\n" "Output dir:"        "$OUTPUT_DIR"
printf "  %-22s %s\n" "Optional flag:"     "${OPTIONAL_FLAG:-<none>}"
printf "  %-22s %s\n" "Workspace:"         "$TEMP_OUTPUT"
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

# Save initial state (fresh run) or re-save to confirm the loaded path (continue)
mkdir -p "$OUTPUT_DIR"
mkdir -p "$TEMP_OUTPUT"
save_state
if ! $IS_CONTINUE; then
    echo -e "  ${CYAN}State saved to $TEMP_OUTPUT/.runstate${NC}"
    echo ""
fi

# ═════════════════════════════════════════════════════════════════════════════
# Process each sample
# ═════════════════════════════════════════════════════════════════════════════

N_SUCCESS=${#COMPLETED_SAMPLES[@]}   # already-done samples count toward success
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

    # Skip samples recorded as completed in state
    _done=false
    for _cs in "${COMPLETED_SAMPLES[@]}"; do
        [[ "$_cs" == "$SAMPLE" ]] && { _done=true; break; }
    done
    if $_done; then
        echo -e "  ${GREEN}Already completed (recorded in state) — skipping.${NC}"
        continue
    fi

    # Verify the input file exists
    if [[ ! -f "$FILES_DIR/$SAMPLE" ]]; then
        echo -e "${RED}  ERROR: Input file not found: $FILES_DIR/$SAMPLE${NC}"
        echo -e "${RED}  Skipping.${NC}"
        N_FAILED=$((N_FAILED+1))
        FAILED_SAMPLES+=("$SAMPLE")
        continue
    fi

    # Each sample gets its own subdir inside the workspace.
    # If the subdir already exists the sample was interrupted — pass --continue
    # so parallelRun re-checks per-file logs and skips successful files.
    SAMPLE_WORK_DIR="$TEMP_OUTPUT/sample_${i}"
    SAMPLE_CONTINUE_FLAG=""
    if [[ -d "$SAMPLE_WORK_DIR" ]]; then
        echo -e "  ${YELLOW}Existing work dir found — passing --continue to parallelRun.sh${NC}"
        SAMPLE_CONTINUE_FLAG="--continue"
    fi

    # Run parallelRun.sh; capture exit code without propagating set -e
    set +e
    bash "$PARALLEL_SCRIPT" \
        "$FILES_DIR/$SAMPLE" \
        -c "$CONFIG_PATH" \
        $COLLECTION_ARG "$TRACK_COLLECTION" \
        -j 8 \
        -o "$SAMPLE_WORK_DIR" \
        -n "$OUTPUT_NAME" \
        $SAMPLE_CONTINUE_FLAG \
        $EXTRA_ARGS
    PARALLEL_EXIT=$?
    set -e

    SAMPLE_FAILED=false
    if [[ $PARALLEL_EXIT -ne 0 ]]; then
        echo -e "${RED}  parallelRun.sh exited with status $PARALLEL_EXIT for $SAMPLE${NC}"
        N_FAILED=$((N_FAILED+1))
        FAILED_SAMPLES+=("$SAMPLE")
        SAMPLE_FAILED=true
    elif [[ -f "$SAMPLE_WORK_DIR/$OUTPUT_NAME" ]]; then
        mv "$SAMPLE_WORK_DIR/$OUTPUT_NAME" "$OUTPUT_DIR/$OUTPUT_NAME"
        echo -e "${GREEN}  Saved: $OUTPUT_DIR/$OUTPUT_NAME${NC}"
        N_SUCCESS=$((N_SUCCESS+1))
        # Record completion in state before cleaning up the work dir
        COMPLETED_SAMPLES+=("$SAMPLE")
        save_state
    else
        echo -e "${RED}  WARNING: parallelRun.sh succeeded but merged output not found: $SAMPLE_WORK_DIR/$OUTPUT_NAME${NC}"
        N_FAILED=$((N_FAILED+1))
        FAILED_SAMPLES+=("$SAMPLE")
        SAMPLE_FAILED=true
    fi

    # On success, clean the per-sample work dir (output already moved to OUTPUT_DIR).
    # On failure, preserve it so logs and partial ntuples can be inspected.
    if [[ -d "$SAMPLE_WORK_DIR" ]]; then
        if $SAMPLE_FAILED; then
            echo -e "${YELLOW}  Preserving $SAMPLE_WORK_DIR for log inspection.${NC}"
            echo -e "${YELLOW}  Resume with: bash scripts/runSampleSet.sh --continue $(basename "$TEMP_OUTPUT")${NC}"
        else
            echo "  Removing $SAMPLE_WORK_DIR..."
            rm -rf "$SAMPLE_WORK_DIR"
        fi
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
