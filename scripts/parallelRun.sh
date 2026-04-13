#!/bin/bash
#
# parallelRun.sh - Run cmsRun in parallel on multiple input files
#
# Processes files from a text file list in groups, merges outputs at the end.
#
# Usage:
#   ./parallelRun.sh <file_list.txt> [options]
#
# Examples:
#   ./parallelRun.sh myfiles.txt
#   ./parallelRun.sh myfiles.txt -j 8 -c testHyddraSVAnalyzer_cfg.py
#   ./parallelRun.sh myfiles.txt -j 4 -o my_output --track-collection general
#   ./parallelRun.sh myfiles.txt -c testTrackAnalyzer_miniAOD_cfg.py --apply-cuts --min-pt 2.0
#

set -e

# Default values
N_JOBS=8
FILES_PER_JOB=1
CONFIG="testHyddraSVAnalyzer_cfg.py"
OUTPUT_DIR="parallel_output"
OUTPUT_NAME="merged_ntuple.root"
EXTRA_ARGS=""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# ── CMSSW environment check ───────────────────────────────────────────────────
check_cmsenv() {
    if [[ -n "$CMSSW_BASE" ]]; then
        return 0
    fi

    # Walk up the directory tree from the script location to find a .SCRAM dir
    local dir
    dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    local cmssw_base=""
    while [[ "$dir" != "/" ]]; do
        if [[ -d "$dir/.SCRAM" ]]; then
            cmssw_base="$dir"
            break
        fi
        dir="$(dirname "$dir")"
    done

    if [[ -z "$cmssw_base" ]]; then
        echo -e "${RED}Error: cmsenv has not been run and no CMSSW release area was found.${NC}"
        echo "  cd to your CMSSW release directory and run: cmsenv"
        exit 1
    fi

    if ! command -v scram &>/dev/null; then
        echo -e "${RED}Error: cmsenv has not been run and 'scram' is not in PATH.${NC}"
        echo "  Source the CMS software environment, then: cd $cmssw_base && cmsenv"
        exit 1
    fi

    echo -e "${YELLOW}[cmsenv] Auto-detected CMSSW area: $cmssw_base${NC}"
    eval "$(cd "$cmssw_base" && scram runtime -sh)"
    echo -e "${GREEN}[cmsenv] Environment ready (CMSSW_BASE=$CMSSW_BASE)${NC}"
    echo ""
}
check_cmsenv

print_usage() {
    echo "Usage: $0 <file_list.txt> [options]"
    echo ""
    echo "Arguments:"
    echo "  file_list.txt    Text file containing input ROOT files (one per line)"
    echo ""
    echo "Options:"
    echo "  -j, --jobs N          Number of parallel jobs (default: 8)
  --files-per-job N     Input files to process per job (default: 1)"
    echo "  -c, --config FILE     cmsRun config file (default: testHyddraSVAnalyzer_cfg.py)"
    echo "  -o, --output-dir DIR  Output directory (default: parallel_output)"
    echo "  -n, --output-name     Merged output filename (default: merged_ntuple.root)"
    echo "  --continue            Skip files that completed successfully in a previous run"
    echo "  --track-collection X  Track collection option to pass to config"
    echo "  --mother-pdg-id X     PDG ID of the signal mother particle (e.g. 1000023 for iDM)"
    echo "  --hyddra-preset X     HYDDRA leptonic preset: default, NonIso, TightIso"
    echo "  --process-mode MODE   SV type: both (default), leptonic, or hadronic"
    echo "  --collection X        Vertex collection (e.g. PatMuonVertex, PatDSAMuonVertex)"
    echo "  --no-gen              Disable gen info (for data)"
    echo "  --no-merge            Skip hadd merge step"
    echo "  -h, --help            Show this help message"
    echo ""
    echo "Track cut options (for MiniAOD configs):"
    echo "  --apply-cuts          Enable track quality cuts"
    echo "  --min-pt VALUE        Minimum track pT in GeV (default: 1.0)"
    echo "  --min-abs-sip2d VALUE Minimum |sip2D| for displaced selection (default: 4.0)"
    echo "  --max-chi2 VALUE      Maximum normalized chi2 (default: 5.0)"
    echo ""
    echo "Examples:"
    echo "  $0 files.txt"
    echo "  $0 files.txt -j 4 -o results"
    echo "  $0 files.txt --track-collection general --no-gen"
    echo "  $0 files.txt -c testTrackAnalyzer_miniAOD_cfg.py --apply-cuts --min-pt 2.0"
}

# Check for GNU parallel
check_parallel() {
    if ! command -v parallel &> /dev/null; then
        echo -e "${RED}Error: GNU parallel is not installed.${NC}"
        echo "Install with:"
        echo "  macOS:  brew install parallel"
        echo "  Ubuntu: sudo apt-get install parallel"
        exit 1
    fi
}

# Parse arguments
if [[ $# -lt 1 ]]; then
    print_usage
    exit 1
fi

INPUT_LIST="$1"
shift

DO_MERGE=true
CONTINUE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -j|--jobs)
            N_JOBS="$2"
            shift 2
            ;;
        --files-per-job)
            FILES_PER_JOB="$2"
            shift 2
            ;;
        -c|--config)
            CONFIG="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -n|--output-name)
            OUTPUT_NAME="$2"
            shift 2
            ;;
        --track-collection)
            EXTRA_ARGS="$EXTRA_ARGS trackCollection=$2"
            shift 2
            ;;
        --mother-pdg-id)
            EXTRA_ARGS="$EXTRA_ARGS motherPdgId=$2"
            shift 2
            ;;
        --gen-dr-cut)
            EXTRA_ARGS="$EXTRA_ARGS genDRCut=$2"
            shift 2
            ;;
        --max-norm-chi2)
            EXTRA_ARGS="$EXTRA_ARGS maxNormChi2=$2"
            shift 2
            ;;
        --hyddra-preset)
            EXTRA_ARGS="$EXTRA_ARGS hyddraPreset=$2"
            shift 2
            ;;
        --collection)
            EXTRA_ARGS="$EXTRA_ARGS collection=$2"
            shift 2
            ;;
        --process-mode)
            EXTRA_ARGS="$EXTRA_ARGS processMode=$2"
            shift 2
            ;;
        --no-gen)
            EXTRA_ARGS="$EXTRA_ARGS hasGenInfo=False"
            shift
            ;;
        --no-merge)
            DO_MERGE=false
            shift
            ;;
        --continue)
            CONTINUE=true
            shift
            ;;
        --apply-cuts)
            EXTRA_ARGS="$EXTRA_ARGS applyCuts=True"
            shift
            ;;
        --min-pt)
            EXTRA_ARGS="$EXTRA_ARGS minPt=$2"
            shift 2
            ;;
        --min-abs-sip2d)
            EXTRA_ARGS="$EXTRA_ARGS minAbsSip2D=$2"
            shift 2
            ;;
        --max-chi2)
            EXTRA_ARGS="$EXTRA_ARGS maxNormalizedChi2=$2"
            shift 2
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            print_usage
            exit 1
            ;;
    esac
done

# Validate inputs
if [[ ! -f "$INPUT_LIST" ]]; then
    echo -e "${RED}Error: Input file list not found: $INPUT_LIST${NC}"
    exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
    # Try in test/ directory
    if [[ -f "test/$CONFIG" ]]; then
        CONFIG="test/$CONFIG"
    else
        echo -e "${RED}Error: Config file not found: $CONFIG${NC}"
        exit 1
    fi
fi

check_parallel

# Count files (excluding comments and empty lines)
N_FILES=$(grep -v '^#' "$INPUT_LIST" | grep -v '^$' | wc -l | tr -d ' ')

if [[ $N_FILES -eq 0 ]]; then
    echo -e "${RED}Error: No files found in $INPUT_LIST${NC}"
    exit 1
fi

# Setup output directory
mkdir -p "$OUTPUT_DIR"
LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$LOG_DIR"

# Build the list of files to process
FILES_TO_PROCESS=$(mktemp)
grep -v '^#' "$INPUT_LIST" | grep -v '^$' > "$FILES_TO_PROCESS"

N_SKIPPED=0
if [[ "$CONTINUE" == true ]]; then
    FILTERED=$(mktemp)
    while IFS= read -r INPUT_FILE; do
        BASENAME=$(basename "$INPUT_FILE" .root)
        BASENAME=${BASENAME#file:}
        LOG_FILE="$LOG_DIR/${BASENAME}.log"

        if [[ -f "$LOG_FILE" ]] && grep -q "CMSRUN_EXIT_SUCCESS" "$LOG_FILE"; then
            N_SKIPPED=$((N_SKIPPED + 1))
        else
            echo "$INPUT_FILE" >> "$FILTERED"
        fi
    done < "$FILES_TO_PROCESS"
    mv "$FILTERED" "$FILES_TO_PROCESS"
fi

N_TO_PROCESS=$(wc -l < "$FILES_TO_PROCESS" | tr -d ' ')

echo -e "${GREEN}============================================${NC}"
echo -e "${GREEN}  Parallel cmsRun Processing${NC}"
echo -e "${GREEN}============================================${NC}"
echo ""
echo "  Input file list: $INPUT_LIST"
echo "  Total files:     $N_FILES"
if [[ "$CONTINUE" == true ]]; then
echo "  Already done:    $N_SKIPPED"
fi
echo "  Files to run:    $N_TO_PROCESS"
echo "  Files per job:   $FILES_PER_JOB"
echo "  Parallel jobs:   $N_JOBS"
echo "  Config file:     $CONFIG"
echo "  Output dir:      $OUTPUT_DIR"
echo "  Extra args:      $EXTRA_ARGS"
echo ""

if [[ $N_TO_PROCESS -eq 0 ]]; then
    echo -e "${GREEN}All files already processed. Nothing to do.${NC}"
    rm -f "$FILES_TO_PROCESS"
    # Still merge if requested
    N_OUTPUTS=$(ls -1 "$OUTPUT_DIR"/*_ntuple.root 2>/dev/null | wc -l | tr -d ' ')
    if [[ "$DO_MERGE" == true ]] && [[ $N_OUTPUTS -gt 0 ]]; then
        echo ""
        echo -e "${YELLOW}Merging outputs with hadd...${NC}"
        MERGED_FILE="$OUTPUT_DIR/$OUTPUT_NAME"
        if hadd -f "$MERGED_FILE" "$OUTPUT_DIR"/*_ntuple.root; then
            echo -e "${GREEN}Merged output: $MERGED_FILE${NC}"
            MERGED_SIZE=$(ls -lh "$MERGED_FILE" | awk '{print $5}')
            echo "  Size: $MERGED_SIZE"
        else
            echo -e "${RED}hadd merge failed${NC}"
            exit 1
        fi
    fi
    echo ""
    echo -e "${GREEN}Done!${NC}"
    exit 0
fi

# Split FILES_TO_PROCESS into chunks of FILES_PER_JOB lines each.
# Each chunk is a temp file whose path is fed to parallel as one job.
CHUNK_DIR=$(mktemp -d)
CHUNK_LIST=$(mktemp)
chunk_idx=0
line_count=0
current_chunk=""

while IFS= read -r f; do
    if [[ $line_count -eq 0 ]]; then
        current_chunk="$CHUNK_DIR/chunk_$(printf '%04d' $chunk_idx)"
        echo "$current_chunk" >> "$CHUNK_LIST"
        chunk_idx=$((chunk_idx + 1))
    fi
    echo "$f" >> "$current_chunk"
    line_count=$((line_count + 1))
    if [[ $line_count -ge $FILES_PER_JOB ]]; then
        line_count=0
        current_chunk=""
    fi
done < "$FILES_TO_PROCESS"

N_CHUNKS=$(wc -l < "$CHUNK_LIST" | tr -d ' ')

# Create a temporary script that processes all files in a chunk sequentially
TEMP_SCRIPT=$(mktemp)
cat > "$TEMP_SCRIPT" << 'SCRIPT_EOF'
#!/bin/bash
CHUNK_FILE="$1"
CONFIG="$2"
OUTPUT_DIR="$3"
EXTRA_ARGS="$4"

CHUNK_FAILED=0
while IFS= read -r INPUT_FILE; do
    [[ -z "$INPUT_FILE" ]] && continue
    BASENAME=$(basename "$INPUT_FILE" .root)
    BASENAME=${BASENAME#file:}
    OUTPUT_FILE="$OUTPUT_DIR/${BASENAME}_ntuple.root"
    LOG_FILE="$OUTPUT_DIR/logs/${BASENAME}.log"

    echo "[$(date '+%H:%M:%S')] Starting: $BASENAME"

    if cmsRun "$CONFIG" inputFiles="$INPUT_FILE" outputFile="$OUTPUT_FILE" $EXTRA_ARGS > "$LOG_FILE" 2>&1; then
        echo "CMSRUN_EXIT_SUCCESS" >> "$LOG_FILE"
        echo "[$(date '+%H:%M:%S')] Completed: $BASENAME"
    else
        echo "[$(date '+%H:%M:%S')] FAILED: $BASENAME (see $LOG_FILE)"
        CHUNK_FAILED=1
    fi
done < "$CHUNK_FILE"
exit $CHUNK_FAILED
SCRIPT_EOF
chmod +x "$TEMP_SCRIPT"

# Run parallel processing
echo -e "${YELLOW}Starting parallel processing ($N_CHUNKS chunks, $FILES_PER_JOB file(s)/job)...${NC}"
echo ""

START_TIME=$(date +%s)

cat "$CHUNK_LIST" | \
    parallel --bar -j "$N_JOBS" "$TEMP_SCRIPT" {} "$CONFIG" "$OUTPUT_DIR" "'$EXTRA_ARGS'"

PARALLEL_EXIT=$?

# Cleanup temp files
rm -f "$TEMP_SCRIPT" "$FILES_TO_PROCESS" "$CHUNK_LIST"
rm -rf "$CHUNK_DIR"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

# Count successful outputs — more reliable than parallel's exit code, which may
# return the last job's status rather than a non-zero aggregate.
N_OUTPUTS=$(ls -1 "$OUTPUT_DIR"/*_ntuple.root 2>/dev/null | wc -l | tr -d ' ')
N_EXPECTED=$((N_TO_PROCESS + N_SKIPPED))
if [[ $N_OUTPUTS -lt $N_EXPECTED ]]; then
    PARALLEL_EXIT=1
fi

echo ""
if [[ $PARALLEL_EXIT -ne 0 ]]; then
    echo -e "${RED}Some jobs failed. Check logs in $LOG_DIR${NC}"
fi

echo -e "${GREEN}Processing completed in ${ELAPSED}s${NC}"
echo "  Successful outputs: $N_OUTPUTS / $N_FILES (${N_SKIPPED} from previous run)"

if [[ "$DO_MERGE" == true ]] && [[ $N_OUTPUTS -eq 0 ]]; then
    echo -e "${RED}No output files produced — all jobs failed. Check logs in $LOG_DIR${NC}"
    exit 1
fi

# Merge outputs
if [[ "$DO_MERGE" == true ]] && [[ $N_OUTPUTS -gt 0 ]]; then
    echo ""
    echo -e "${YELLOW}Merging outputs with hadd...${NC}"

    MERGED_FILE="$OUTPUT_DIR/$OUTPUT_NAME"

    if hadd -f "$MERGED_FILE" "$OUTPUT_DIR"/*_ntuple.root; then
        echo -e "${GREEN}Merged output: $MERGED_FILE${NC}"

        # Show file size
        MERGED_SIZE=$(ls -lh "$MERGED_FILE" | awk '{print $5}')
        echo "  Size: $MERGED_SIZE"
    else
        echo -e "${RED}hadd merge failed${NC}"
        exit 1
    fi
fi

echo ""
echo -e "${GREEN}Done!${NC}"
exit $PARALLEL_EXIT
