#!/bin/bash

#########################################################################
#                                                                       #
#  Bash script for running TS reoptimization from existing NEB         #
#                                                                       #
#########################################################################

# Default values
NPROCS=102
ORCA_PATH="/home/rudenko/orca/orca"
WORK_DIR=$(pwd)
SPECIFIC_DIR=""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Help function
show_help() {
    cat << EOF
Usage: ./run_reoptimize.sh [OPTIONS]

This script reoptimizes TS from existing NEB calculations and recalculates frequencies.

OPTIONS:
    -h, --help              Show this help message
    -p, --path PATH         Path to working directory (default: current directory)
    -o, --orca PATH         Path to ORCA executable (default: /home/rudenko/orca/orca)
    -n, --nprocs NUM        Number of processors (default: 102)
    -d, --dir NAME          Process only specific directory (default: all directories)
    
EXAMPLES:
    # Process all directories with default settings
    ./run_reoptimize.sh
    
    # Process all directories with custom ORCA path and 64 processors
    ./run_reoptimize.sh -o /opt/orca/orca -n 64
    
    # Process only specific directory
    ./run_reoptimize.sh -d CH3NO2_H2O_near
    
    # Process from different working directory
    ./run_reoptimize.sh -p /path/to/calculations -n 48

EOF
    exit 0
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            ;;
        -p|--path)
            WORK_DIR="$2"
            shift 2
            ;;
        -o|--orca)
            ORCA_PATH="$2"
            shift 2
            ;;
        -n|--nprocs)
            NPROCS="$2"
            shift 2
            ;;
        -d|--dir)
            SPECIFIC_DIR="$2"
            shift 2
            ;;
        *)
            echo -e "${RED}Error: Unknown option $1${NC}"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Check if working directory exists
if [ ! -d "$WORK_DIR" ]; then
    echo -e "${RED}Error: Working directory $WORK_DIR does not exist!${NC}"
    exit 1
fi

# Check if for_orca directory exists
if [ ! -d "$WORK_DIR/for_orca" ]; then
    echo -e "${RED}Error: for_orca directory not found in $WORK_DIR${NC}"
    exit 1
fi

# Check if ORCA executable exists
if [ ! -f "$ORCA_PATH" ]; then
    echo -e "${RED}Error: ORCA executable not found at $ORCA_PATH${NC}"
    exit 1
fi

# Check if Python script exists
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="$SCRIPT_DIR/reoptimize_ts.py"

if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo -e "${RED}Error: Python script not found at $PYTHON_SCRIPT${NC}"
    exit 1
fi

# Print configuration
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}TS Reoptimization Script${NC}"
echo -e "${GREEN}========================================${NC}"
echo -e "Working directory: ${YELLOW}$WORK_DIR${NC}"
echo -e "ORCA path: ${YELLOW}$ORCA_PATH${NC}"
echo -e "Number of processors: ${YELLOW}$NPROCS${NC}"
if [ -n "$SPECIFIC_DIR" ]; then
    echo -e "Processing directory: ${YELLOW}$SPECIFIC_DIR${NC}"
else
    echo -e "Processing: ${YELLOW}ALL directories${NC}"
fi
echo -e "${GREEN}========================================${NC}"
echo ""

# Ask for confirmation
read -p "Continue? (y/n) " -n 1 -r
echo
if [[ ! $REPL =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 0
fi

# Build Python command
CMD="python3 $PYTHON_SCRIPT -p $WORK_DIR -o $ORCA_PATH -n $NPROCS"

if [ -n "$SPECIFIC_DIR" ]; then
    CMD="$CMD -d $SPECIFIC_DIR"
fi

# Run the script
echo -e "${GREEN}Starting reoptimization...${NC}"
echo "Command: $CMD"
echo ""

$CMD

EXIT_CODE=$?

# Check exit status
if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}Reoptimization completed successfully!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo -e "Check ${YELLOW}$WORK_DIR/reoptimization_times.txt${NC} for details"
else
    echo ""
    echo -e "${RED}========================================${NC}"
    echo -e "${RED}Reoptimization failed with exit code $EXIT_CODE${NC}"
    echo -e "${RED}========================================${NC}"
    echo -e "Check ${YELLOW}$WORK_DIR/reoptimization_times.txt${NC} for details"
fi

exit $EXIT_CODE
