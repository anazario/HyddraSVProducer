#!/usr/bin/env python3
"""
MultiCrab submission script for HyddraSVAnalyzer

Usage:
    # Submit with a filelist (request name derived from filename: "my_sample")
    python crab_hyddraSVAnalyzer.py --filelist my_sample.txt

    # Submit with a dataset (request name derived from primary dataset name)
    python crab_hyddraSVAnalyzer.py --dataset /SMS-GlGl.../AODSIM

    # Override the auto-derived name
    python crab_hyddraSVAnalyzer.py --filelist files.txt --name custom_name

    # Append a trial tag to differentiate runs
    python crab_hyddraSVAnalyzer.py --filelist my_sample.txt --trial v2

    # Process only leptonic vertices
    python crab_hyddraSVAnalyzer.py --filelist my_sample.txt --mode leptonic

    # Dry run (don't actually submit)
    python crab_hyddraSVAnalyzer.py --filelist my_sample.txt --dry-run

    # Use a specific track collection and vertex collection
    python crab_hyddraSVAnalyzer.py --filelist my_sample.txt --track-collection general --collection PatDSAMuonVertex

    # Apply track quality cuts
    python crab_hyddraSVAnalyzer.py --filelist my_sample.txt --apply-cuts --min-pt 2.0 --max-chi2 3.0
"""

import os
import sys
import argparse
from datetime import datetime

def get_crab_config(args):
    """Generate CRAB configuration based on arguments."""

    from CRABClient.UserUtilities import config

    cfg = config()

    # General settings
    cfg.General.requestName = args.name
    cfg.General.workArea = args.workarea
    cfg.General.transferOutputs = True
    cfg.General.transferLogs = True

    # Job type settings
    cfg.JobType.pluginName = 'Analysis'
    cfg.JobType.psetName = os.path.join(os.path.dirname(__file__), 'testHyddraSVAnalyzer_cfg.py')
    cfg.JobType.maxMemoryMB = args.memory
    cfg.JobType.numCores = args.cores

    # Build pyCfgParams based on mode
    pycfg_params = []
    pycfg_params.append(f'hasGenInfo={args.gen_info}')
    pycfg_params.append(f'processMode={args.mode}')

    if args.max_events > 0:
        pycfg_params.append(f'maxEvents={args.max_events}')

    if args.track_collection:
        pycfg_params.append(f'trackCollection={args.track_collection}')
    if args.collection:
        pycfg_params.append(f'collection={args.collection}')
    if args.apply_cuts:
        pycfg_params.append('applyCuts=True')
    if args.min_pt is not None:
        pycfg_params.append(f'minPt={args.min_pt}')
    if args.min_abs_sip2d is not None:
        pycfg_params.append(f'minAbsSip2D={args.min_abs_sip2d}')
    if args.max_chi2 is not None:
        pycfg_params.append(f'maxNormalizedChi2={args.max_chi2}')

    cfg.JobType.pyCfgParams = pycfg_params

    # Data settings
    if args.dataset:
        cfg.Data.inputDataset = args.dataset
        cfg.Data.inputDBS = args.dbs
    elif args.filelist:
        cfg.Data.userInputFiles = read_filelist(args.filelist)
        cfg.Data.splitting = 'FileBased'
        cfg.Data.unitsPerJob = args.files_per_job

    if not args.filelist:
        cfg.Data.splitting = 'FileBased'
        cfg.Data.unitsPerJob = args.files_per_job

    cfg.Data.outLFNDirBase = args.output_dir
    cfg.Data.publication = False
    cfg.Data.outputDatasetTag = args.name

    # Site settings
    cfg.Site.storageSite = args.storage_site

    if args.whitelist:
        cfg.Site.whitelist = args.whitelist.split(',')
    if args.blacklist:
        cfg.Site.blacklist = args.blacklist.split(',')

    return cfg


def read_filelist(filepath):
    """Read a txt file containing one root file path per line."""
    files = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                files.append(line)
    return files


def create_modified_pset(args):
    """Create a modified pset that handles the processing mode."""

    pset_template = '''
# Auto-generated CRAB pset wrapper for HyddraSVAnalyzer
# Processing mode: {mode}

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# Setup command line options
options = VarParsing.VarParsing('analysis')
options.register('hasGenInfo',
                 {has_gen_info},
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Process gen-level information (set False for data)")
options.register('processMode',
                 '{mode}',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Processing mode: both, leptonic, or hadronic")
options.setDefault('maxEvents', -1)
options.setDefault('outputFile', 'hyddraSV_ntuple.root')
options.parseArguments()

process = cms.Process("HYDDRAANALYSIS")

# Message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        options.inputFiles if options.inputFiles else ['file:input.root']
    )
)

# Load geometry and magnetic field
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# Remove CASTOR if present (not available in all geometries)
if hasattr(process, 'CaloGeometryBuilder') and 'CASTOR' in process.CaloGeometryBuilder.SelectedCalos:
    process.CaloGeometryBuilder.SelectedCalos.remove('CASTOR')

process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

# GlobalTag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '{global_tag}', '')

# TFileService for output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

# Load producers
process.load("KUCMSNtupleizer.KUCMSNtupleizer.MuonEnhancedTracks_cfi")
process.load("KUCMSNtupleizer.KUCMSNtupleizer.ECALTracks_cfi")
process.load("KUCMSNtupleizer.HyddraSVProducer.hyddra_cfi")
process.load("KUCMSNtupleizer.HyddraSVProducer.hyddraSVAnalyzer_cfi")

# Configure analyzer
process.hyddraSVAnalyzer.hasGenInfo = cms.bool(options.hasGenInfo)
process.hyddraSVAnalyzer.mergedSCs = cms.InputTag("ecalTracks", "displacedElectronSCs")

# Configure based on processing mode
if options.processMode == 'leptonic':
    # Only process leptonic vertices - set hadronic to empty
    process.hyddraSVAnalyzer.hadronicVertices = cms.InputTag("")
elif options.processMode == 'hadronic':
    # Only process hadronic vertices - set leptonic to empty
    process.hyddraSVAnalyzer.leptonicVertices = cms.InputTag("")
# else: process both (default)

# Path
process.p = cms.Path(
    process.muonEnhancedTracks +
    process.ecalTracks +
    process.hyddraSVs +
    process.hyddraSVAnalyzer
)

process.schedule = cms.Schedule(process.p)
'''

    return pset_template.format(
        mode=args.mode,
        has_gen_info=str(args.gen_info),
        global_tag=args.global_tag
    )


def main():
    parser = argparse.ArgumentParser(
        description='MultiCrab submission for HyddraSVAnalyzer',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--dataset', '-d',
                            help='Published dataset name (e.g., /SMS-GlGl.../AODSIM)')
    input_group.add_argument('--filelist', '-f',
                            help='Text file containing list of root files')

    # Naming options
    parser.add_argument('--name', '-n', default=None,
                       help='Request name (default: derived from filelist basename)')
    parser.add_argument('--trial', '-t', default=None,
                       help='Trial name appended to request name (e.g. v1, noCuts)')

    # Processing mode
    parser.add_argument('--mode', '-m',
                       choices=['both', 'leptonic', 'hadronic'],
                       default='both',
                       help='Processing mode: both (default), leptonic, or hadronic')

    # Gen info
    parser.add_argument('--gen-info',
                       type=lambda x: x.lower() == 'true',
                       default=True,
                       help='Process gen-level information (default: True)')
    parser.add_argument('--data', action='store_true',
                       help='Shortcut for --gen-info=False')

    # Job configuration
    parser.add_argument('--files-per-job', type=int, default=1,
                       help='Number of files per job (default: 1)')
    parser.add_argument('--max-events', type=int, default=-1,
                       help='Max events to process per job (default: -1 = all)')
    parser.add_argument('--memory', type=int, default=2500,
                       help='Memory per job in MB (default: 2500)')
    parser.add_argument('--cores', type=int, default=1,
                       help='Number of cores per job (default: 1)')

    # Output configuration
    parser.add_argument('--output-dir', '-o',
                       default='/store/group/lpcsusylep/anazario/HyddraSVs',
                       help='Output LFN directory base')
    parser.add_argument('--workarea', '-w',
                       default='crab_projects_{date}'.format(
                           date=datetime.now().strftime('%Y%m%d')),
                       help='CRAB work area directory')

    # Site configuration
    parser.add_argument('--storage-site', '-s',
                       default='T3_US_FNALLPC',
                       help='Storage site (default: T3_US_FNALLPC)')
    parser.add_argument('--whitelist',
                       default='T3_US_FNALLPC',
                       help='Comma-separated list of whitelisted sites (default: T3_US_FNALLPC)')
    parser.add_argument('--blacklist',
                       help='Comma-separated list of blacklisted sites')

    # Database
    parser.add_argument('--dbs', default='global',
                       help='DBS instance (default: global)')

    # Track and vertex options
    parser.add_argument('--track-collection',
                       help='Track collection option to pass to config')
    parser.add_argument('--collection',
                       help='Vertex collection (e.g. PatMuonVertex, PatDSAMuonVertex)')

    # Track cut options
    parser.add_argument('--apply-cuts', action='store_true',
                       help='Enable track quality cuts')
    parser.add_argument('--min-pt', type=float, default=None,
                       help='Minimum track pT in GeV (default: 1.0)')
    parser.add_argument('--min-abs-sip2d', type=float, default=None,
                       help='Minimum |sip2D| for displaced selection (default: 4.0)')
    parser.add_argument('--max-chi2', type=float, default=None,
                       help='Maximum normalized chi2 (default: 5.0)')

    # GlobalTag
    parser.add_argument('--global-tag', '-g',
                       default='124X_mcRun3_2022_realistic_postEE_v1',
                       help='GlobalTag to use')

    # Execution options
    parser.add_argument('--dry-run', action='store_true',
                       help='Print configuration without submitting')
    parser.add_argument('--generate-pset', action='store_true',
                       help='Generate modified pset file for debugging')

    args = parser.parse_args()

    # Handle --data shortcut
    if args.data:
        args.gen_info = False

    # Derive request name from filelist if not provided
    if args.name is None:
        if args.filelist:
            args.name = os.path.splitext(os.path.basename(args.filelist))[0]
        elif args.dataset:
            # Use the primary dataset name (first component after leading /)
            args.name = args.dataset.strip('/').split('/')[0]

    # Append trial name if provided
    if args.trial:
        args.name = f'{args.name}_{args.trial}'

    # Validate filelist if provided
    if args.filelist and not os.path.exists(args.filelist):
        print(f"Error: File list not found: {args.filelist}")
        sys.exit(1)

    # Generate modified pset if requested
    if args.generate_pset:
        pset_content = create_modified_pset(args)
        pset_filename = f'pset_{args.name}_{args.mode}.py'
        with open(pset_filename, 'w') as f:
            f.write(pset_content)
        print(f"Generated pset: {pset_filename}")
        if not args.dry_run:
            return

    # Print configuration summary
    print("=" * 60)
    print("HyddraSVAnalyzer CRAB Submission")
    print("=" * 60)
    print(f"Request name:     {args.name}")
    print(f"Processing mode:  {args.mode}")
    print(f"Has gen info:     {args.gen_info}")
    print(f"Files per job:    {args.files_per_job}")
    if args.dataset:
        print(f"Dataset:          {args.dataset}")
    else:
        files = read_filelist(args.filelist)
        print(f"File list:        {args.filelist} ({len(files)} files)")
    print(f"Storage site:     {args.storage_site}")
    print(f"Output dir:       {args.output_dir}")
    print(f"GlobalTag:        {args.global_tag}")
    if args.track_collection:
        print(f"Track collection: {args.track_collection}")
    if args.collection:
        print(f"Vertex collection:{args.collection}")
    if args.apply_cuts:
        print(f"Apply cuts:       True")
        if args.min_pt is not None:
            print(f"  Min pT:         {args.min_pt}")
        if args.min_abs_sip2d is not None:
            print(f"  Min |sip2D|:    {args.min_abs_sip2d}")
        if args.max_chi2 is not None:
            print(f"  Max chi2:       {args.max_chi2}")
    print("=" * 60)

    if args.dry_run:
        print("\n[DRY RUN] Would submit with above configuration")
        print("\nTo generate the pset for inspection, add --generate-pset")
        return

    # Import CRAB and submit
    try:
        from CRABAPI.RawCommand import crabCommand
        from CRABClient.ClientExceptions import ClientException
        from http.client import HTTPException
    except ImportError:
        print("Error: CRAB environment not set up. Run:")
        print("  source /cvmfs/cms.cern.ch/common/crab-setup.sh")
        sys.exit(1)

    # Get configuration
    cfg = get_crab_config(args)

    # Submit
    try:
        print("\nSubmitting to CRAB...")
        crabCommand('submit', config=cfg)
        print("\nSubmission successful!")
    except HTTPException as e:
        print(f"\nHTTP Error: {e}")
        sys.exit(1)
    except ClientException as e:
        print(f"\nCRAB Client Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
