#!/bin/bash
#
# fetchDASFileLists.sh - Query DAS for iDM dielectron MiniAOD signal samples
#                        and write file lists to files/
#
# Requires a valid grid proxy:  voms-proxy-init --voms cms --valid 168:00
#
# Usage:
#   bash scripts/fetchDASFileLists.sh
#   bash scripts/fetchDASFileLists.sh --redirector root://cmsxrootd.fnal.gov/
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FILES_DIR="$(cd "$SCRIPT_DIR/../files" && pwd)"

REDIRECTOR="root://cmsxrootd.fnal.gov/"
if [[ "$1" == "--redirector" && -n "$2" ]]; then
    REDIRECTOR="$2"
fi

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Check dependencies
if ! command -v dasgoclient &> /dev/null; then
    echo -e "${RED}Error: dasgoclient not found. Source a CMSSW environment first.${NC}"
    exit 1
fi

if ! voms-proxy-info --exists &> /dev/null; then
    echo -e "${RED}Error: No valid grid proxy. Run: voms-proxy-init --voms cms --valid 168:00${NC}"
    exit 1
fi

# ---------------------------------------------------------------------------
# Sample definitions: output filename → DAS dataset path
# ---------------------------------------------------------------------------
declare -A SAMPLES

SAMPLES["iDM_DarkPhotonToEE_Mchi-5p5_dMchi-1p0_ctau-1_MiniAODSIM.txt"]="/iDM_DarkPhotonToEE_Mchi-5p5_dMchi-1p0_ctau-1_mA-15p0_HT80_TuneCP5_13TeV_madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3/MINIAODSIM"
SAMPLES["iDM_DarkPhotonToEE_Mchi-5p5_dMchi-1p0_ctau-10_MiniAODSIM.txt"]="/iDM_DarkPhotonToEE_Mchi-5p5_dMchi-1p0_ctau-10_mA-15p0_HT80_TuneCP5_13TeV_madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3/MINIAODSIM"
SAMPLES["iDM_DarkPhotonToEE_Mchi-5p5_dMchi-1p0_ctau-100_MiniAODSIM.txt"]="/iDM_DarkPhotonToEE_Mchi-5p5_dMchi-1p0_ctau-100_mA-15p0_HT80_TuneCP5_13TeV_madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3/MINIAODSIM"
SAMPLES["iDM_DarkPhotonToEE_Mchi-55p0_dMchi-10p0_ctau-1_MiniAODSIM.txt"]="/iDM_DarkPhotonToEE_Mchi-55p0_dMchi-10p0_ctau-1_mA-150p0_HT80_TuneCP5_13TeV_madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3/MINIAODSIM"
SAMPLES["iDM_DarkPhotonToEE_Mchi-55p0_dMchi-10p0_ctau-10_MiniAODSIM.txt"]="/iDM_DarkPhotonToEE_Mchi-55p0_dMchi-10p0_ctau-10_mA-150p0_HT80_TuneCP5_13TeV_madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3/MINIAODSIM"
SAMPLES["iDM_DarkPhotonToEE_Mchi-55p0_dMchi-10p0_ctau-100_MiniAODSIM.txt"]="/iDM_DarkPhotonToEE_Mchi-55p0_dMchi-10p0_ctau-100_mA-150p0_HT80_TuneCP5_13TeV_madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3/MINIAODSIM"

# ---------------------------------------------------------------------------
echo ""
echo "Fetching iDM MiniAOD file lists from DAS"
echo "  Redirector: $REDIRECTOR"
echo "  Output dir: $FILES_DIR"
echo ""

N_OK=0
N_FAIL=0

for OUTFILE in "${!SAMPLES[@]}"; do
    DATASET="${SAMPLES[$OUTFILE]}"
    OUTPATH="$FILES_DIR/$OUTFILE"

    printf "  %-70s " "$OUTFILE"

    if dasgoclient --query="file dataset=$DATASET" 2>/dev/null \
        | sort \
        | sed "s|^|${REDIRECTOR}|" \
        > "$OUTPATH"; then

        N_LINES=$(wc -l < "$OUTPATH" | tr -d ' ')
        if [[ $N_LINES -eq 0 ]]; then
            echo -e "${YELLOW}WARNING: 0 files returned${NC}"
            rm -f "$OUTPATH"
            N_FAIL=$((N_FAIL + 1))
        else
            echo -e "${GREEN}${N_LINES} files${NC}"
            N_OK=$((N_OK + 1))
        fi
    else
        echo -e "${RED}FAILED${NC}"
        rm -f "$OUTPATH"
        N_FAIL=$((N_FAIL + 1))
    fi
done

echo ""
echo -e "${GREEN}Done: ${N_OK} / $((N_OK + N_FAIL)) file lists written to $FILES_DIR/${NC}"
if [[ $N_FAIL -gt 0 ]]; then
    echo -e "${RED}  ${N_FAIL} failed — check your proxy and dataset names.${NC}"
fi
echo ""
