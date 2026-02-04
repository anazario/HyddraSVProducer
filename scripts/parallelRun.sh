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
#

set -e

# Default values
N_JOBS=8
CONFIG="testHyddraSVAnalyzer_cfg.py"
OUTPUT_DIR="parallel_output"
OUTPUT_NAME="merged_ntuple.root"
EXTRA_ARGS=""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

print_usage() {
    echo "Usage: $0 <file_list.txt> [options]"
    echo ""
    echo "Arguments:"
    echo "  file_list.txt    Text file containing input ROOT files (one per line)"
    echo ""
    echo "Options:"
    echo "  -j, --jobs N          Number of parallel jobs (default: 8)"
    echo "  -c, --config FILE     cmsRun config file (default: testHyddraSVAnalyzer_cfg.py)"
    echo "  -o, --output-dir DIR  Output directory (default: parallel_output)"
    echo "  -n, --output-name     Merged output filename (default: merged_ntuple.root)"
    echo "  --track-collection X  Track collection option to pass to config"
    echo "  --no-gen              Disable gen info (for data)"
    echo "  --no-merge            Skip hadd merge step"
    echo "  -h, --help            Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 files.txt"
    echo "  $0 files.txt -j 4 -o results"
    echo "  $0 files.txt --track-collection general --no-gen"
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

while [[ $# -gt 0 ]]; do
    case $1 in
        -j|--jobs)
            N_JOBS="$2"
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
        --no-gen)
            EXTRA_ARGS="$EXTRA_ARGS hasGenInfo=False"
            shift
            ;;
        --no-merge)
            DO_MERGE=false
            shift
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

echo -e "${GREEN}============================================${NC}"
echo -e "${GREEN}  Parallel cmsRun Processing${NC}"
echo -e "${GREEN}============================================${NC}"
echo ""
echo "  Input file list: $INPUT_LIST"
echo "  Number of files: $N_FILES"
echo "  Parallel jobs:   $N_JOBS"
echo "  Config file:     $CONFIG"
echo "  Output dir:      $OUTPUT_DIR"
echo "  Extra args:      $EXTRA_ARGS"
echo ""

# Create a temporary script for parallel to run
# This handles the file naming properly
TEMP_SCRIPT=$(mktemp)
cat > "$TEMP_SCRIPT" << 'SCRIPT_EOF'
#!/bin/bash
INPUT_FILE="$1"
CONFIG="$2"
OUTPUT_DIR="$3"
EXTRA_ARGS="$4"

# Extract base name for output
BASENAME=$(basename "$INPUT_FILE" .root)
# Handle file: prefix
BASENAME=${BASENAME#file:}
OUTPUT_FILE="$OUTPUT_DIR/${BASENAME}_ntuple.root"
LOG_FILE="$OUTPUT_DIR/logs/${BASENAME}.log"

echo "[$(date '+%H:%M:%S')] Starting: $BASENAME"

# Run cmsRun
if cmsRun "$CONFIG" inputFiles="$INPUT_FILE" outputFile="$OUTPUT_FILE" $EXTRA_ARGS > "$LOG_FILE" 2>&1; then
    echo "[$(date '+%H:%M:%S')] Completed: $BASENAME"
    exit 0
else
    echo "[$(date '+%H:%M:%S')] FAILED: $BASENAME (see $LOG_FILE)"
    exit 1
fi
SCRIPT_EOF
chmod +x "$TEMP_SCRIPT"

# Run parallel processing
echo -e "${YELLOW}Starting parallel processing...${NC}"
echo ""

START_TIME=$(date +%s)

# Use parallel to process files
# --bar shows progress, --halt soon,fail=1 stops on first failure (optional)
grep -v '^#' "$INPUT_LIST" | grep -v '^$' | \
    parallel --bar -j "$N_JOBS" "$TEMP_SCRIPT" {} "$CONFIG" "$OUTPUT_DIR" "'$EXTRA_ARGS'"

PARALLEL_EXIT=$?

# Cleanup temp script
rm -f "$TEMP_SCRIPT"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
if [[ $PARALLEL_EXIT -ne 0 ]]; then
    echo -e "${RED}Some jobs failed. Check logs in $LOG_DIR${NC}"
fi

echo -e "${GREEN}Processing completed in ${ELAPSED}s${NC}"

# Count successful outputs
N_OUTPUTS=$(ls -1 "$OUTPUT_DIR"/*_ntuple.root 2>/dev/null | wc -l | tr -d ' ')
echo "  Successful outputs: $N_OUTPUTS / $N_FILES"

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
