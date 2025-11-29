#!/usr/bin/env python3
#########################################################################
#                                                                       #
#  Script to check status of TS reoptimization calculations            #
#                                                                       #
#########################################################################

import os
import sys
import argparse
from pathlib import Path


def check_file_exists(filepath):
    """Check if file exists"""
    return os.path.exists(filepath)


def check_normally_terminated(filepath):
    """Check if ORCA calculation finished normally"""
    if not os.path.exists(filepath):
        return False
    try:
        with open(filepath, 'r', encoding='UTF-8') as f:
            content = f.read()
            return "NORMALLY" in content
    except:
        return False


def check_imaginary_modes(filepath):
    """Count imaginary modes in frequency calculation"""
    if not os.path.exists(filepath):
        return None
    try:
        count = 0
        with open(filepath, 'r', encoding='UTF-8') as f:
            for line in f:
                if "***imaginary mode***" in line:
                    count += 1
        return count
    except:
        return None


def check_directory_status(dir_path):
    """
    Check status of calculations in directory
    Returns dict with status information
    """
    status = {
        'neb_exists': False,
        'ts_opt_exists': False,
        'ts_opt_done': False,
        'freq_H_exists': False,
        'freq_H_done': False,
        'freq_H_imag': None,
        'freq_D_exists': False,
        'freq_D_done': False,
        'freq_D_imag': None,
        'start_freq_H_exists': False,
        'start_freq_H_done': False,
        'start_freq_H_imag': None,
        'start_freq_D_exists': False,
        'start_freq_D_done': False,
        'start_freq_D_imag': None,
        'overall_status': 'NOT_STARTED'
    }
    
    orca_path = os.path.join(dir_path, 'orca')
    
    # Check NEB
    neb_xyz = os.path.join(orca_path, 'neb', 'neb.xyz')
    status['neb_exists'] = check_file_exists(neb_xyz)
    
    # Check TS optimization
    ts_opt_path = os.path.join(orca_path, 'ts_opt')
    ts_opt_out = os.path.join(ts_opt_path, 'ts_opt.out')
    status['ts_opt_exists'] = os.path.exists(ts_opt_path)
    status['ts_opt_done'] = check_normally_terminated(ts_opt_out)
    
    # Check TS frequencies
    freq_H_out = os.path.join(orca_path, 'freq_H', 'f_H.out')
    freq_D_out = os.path.join(orca_path, 'freq_D', 'f_D.out')
    status['freq_H_exists'] = os.path.exists(os.path.join(orca_path, 'freq_H'))
    status['freq_H_done'] = check_normally_terminated(freq_H_out)
    status['freq_H_imag'] = check_imaginary_modes(freq_H_out)
    status['freq_D_exists'] = os.path.exists(os.path.join(orca_path, 'freq_D'))
    status['freq_D_done'] = check_normally_terminated(freq_D_out)
    status['freq_D_imag'] = check_imaginary_modes(freq_D_out)
    
    # Check start frequencies
    start_freq_H_out = os.path.join(orca_path, 'start_freq_H', 'start_f_H.out')
    start_freq_D_out = os.path.join(orca_path, 'start_freq_D', 'start_f_D.out')
    status['start_freq_H_exists'] = os.path.exists(os.path.join(orca_path, 'start_freq_H'))
    status['start_freq_H_done'] = check_normally_terminated(start_freq_H_out)
    status['start_freq_H_imag'] = check_imaginary_modes(start_freq_H_out)
    status['start_freq_D_exists'] = os.path.exists(os.path.join(orca_path, 'start_freq_D'))
    status['start_freq_D_done'] = check_normally_terminated(start_freq_D_out)
    status['start_freq_D_imag'] = check_imaginary_modes(start_freq_D_out)
    
    # Determine overall status
    if not status['neb_exists']:
        status['overall_status'] = 'NO_NEB'
    elif not status['ts_opt_exists']:
        status['overall_status'] = 'NOT_STARTED'
    elif not status['ts_opt_done']:
        status['overall_status'] = 'TS_OPT_RUNNING'
    elif status['ts_opt_done'] and not status['freq_H_done']:
        status['overall_status'] = 'FREQ_RUNNING'
    elif (status['freq_H_done'] and status['freq_D_done'] and 
          status['start_freq_H_done'] and status['start_freq_D_done']):
        # Check imaginary modes
        if status['freq_H_imag'] == 1 and status['freq_D_imag'] == 1:
            if status['start_freq_H_imag'] == 0 and status['start_freq_D_imag'] == 0:
                status['overall_status'] = 'COMPLETED_OK'
            else:
                status['overall_status'] = 'COMPLETED_START_IMAG'
        else:
            status['overall_status'] = 'COMPLETED_TS_IMAG_ERROR'
    else:
        status['overall_status'] = 'IN_PROGRESS'
    
    return status


def print_status_line(name, status, verbose=False):
    """Print status line for directory"""
    status_str = status['overall_status']
    
    # Color codes
    colors = {
        'COMPLETED_OK': '\033[0;32m',      # Green
        'NOT_STARTED': '\033[0;37m',       # White
        'NO_NEB': '\033[0;31m',            # Red
        'TS_OPT_RUNNING': '\033[1;33m',    # Yellow
        'FREQ_RUNNING': '\033[1;33m',      # Yellow
        'IN_PROGRESS': '\033[1;33m',       # Yellow
        'COMPLETED_START_IMAG': '\033[0;35m',  # Magenta
        'COMPLETED_TS_IMAG_ERROR': '\033[0;31m',  # Red
    }
    RESET = '\033[0m'
    
    color = colors.get(status_str, '\033[0;37m')
    
    if not verbose:
        print(f"{color}{name:40s} {status_str:30s}{RESET}")
    else:
        print(f"\n{color}{'='*70}")
        print(f"Directory: {name}")
        print(f"Status: {status_str}")
        print(f"{'='*70}{RESET}")
        print(f"  NEB exists: {status['neb_exists']}")
        print(f"  TS optimization:")
        print(f"    - Folder exists: {status['ts_opt_exists']}")
        print(f"    - Completed: {status['ts_opt_done']}")
        print(f"  TS frequencies:")
        print(f"    - H: Done={status['freq_H_done']}, Imag={status['freq_H_imag']}")
        print(f"    - D: Done={status['freq_D_done']}, Imag={status['freq_D_imag']}")
        print(f"  Start frequencies:")
        print(f"    - H: Done={status['start_freq_H_done']}, Imag={status['start_freq_H_imag']}")
        print(f"    - D: Done={status['start_freq_D_done']}, Imag={status['start_freq_D_imag']}")


def main():
    parser = argparse.ArgumentParser(
        description='Check status of TS reoptimization calculations')
    parser.add_argument('-p', '--path', type=str, required=True,
                       help='Path to directory with for_orca folder')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output with detailed information')
    parser.add_argument('-d', '--dir', type=str, default=None,
                       help='Check specific directory only')
    parser.add_argument('--filter', type=str, default=None,
                       choices=['completed', 'not_started', 'in_progress', 'errors'],
                       help='Filter by status')
    
    args = parser.parse_args()
    
    path_to_orca_dirs = os.path.join(args.path, "for_orca")
    
    if not os.path.exists(path_to_orca_dirs):
        print(f"ERROR: Directory {path_to_orca_dirs} does not exist!")
        sys.exit(1)
    
    # Get directories to check
    if args.dir:
        orca_dirs = [args.dir]
    else:
        orca_dirs = sorted([d for d in os.listdir(path_to_orca_dirs) 
                           if os.path.isdir(os.path.join(path_to_orca_dirs, d))])
    
    # Statistics
    stats = {
        'total': 0,
        'completed_ok': 0,
        'not_started': 0,
        'in_progress': 0,
        'errors': 0,
        'no_neb': 0
    }
    
    # Print header
    print("\n" + "="*70)
    print("TS Reoptimization Status Check")
    print("="*70)
    if not args.verbose:
        print(f"{'Directory':40s} {'Status':30s}")
        print("-"*70)
    
    # Check each directory
    for dir_name in orca_dirs:
        dir_path = os.path.join(path_to_orca_dirs, dir_name)
        
        if not os.path.exists(dir_path):
            continue
        
        status = check_directory_status(dir_path)
        stats['total'] += 1
        
        # Update statistics
        if status['overall_status'] == 'COMPLETED_OK':
            stats['completed_ok'] += 1
        elif status['overall_status'] == 'NOT_STARTED':
            stats['not_started'] += 1
        elif status['overall_status'] == 'NO_NEB':
            stats['no_neb'] += 1
        elif 'ERROR' in status['overall_status'] or 'IMAG' in status['overall_status']:
            stats['errors'] += 1
        else:
            stats['in_progress'] += 1
        
        # Apply filter
        if args.filter:
            if args.filter == 'completed' and status['overall_status'] != 'COMPLETED_OK':
                continue
            elif args.filter == 'not_started' and status['overall_status'] != 'NOT_STARTED':
                continue
            elif args.filter == 'in_progress' and status['overall_status'] not in ['TS_OPT_RUNNING', 'FREQ_RUNNING', 'IN_PROGRESS']:
                continue
            elif args.filter == 'errors' and 'ERROR' not in status['overall_status'] and 'IMAG' not in status['overall_status']:
                continue
        
        print_status_line(dir_name, status, args.verbose)
    
    # Print summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Total directories:    {stats['total']}")
    print(f"Completed OK:         {stats['completed_ok']} ({100*stats['completed_ok']/max(stats['total'],1):.1f}%)")
    print(f"Not started:          {stats['not_started']} ({100*stats['not_started']/max(stats['total'],1):.1f}%)")
    print(f"In progress:          {stats['in_progress']} ({100*stats['in_progress']/max(stats['total'],1):.1f}%)")
    print(f"Errors/Issues:        {stats['errors']} ({100*stats['errors']/max(stats['total'],1):.1f}%)")
    print(f"No NEB data:          {stats['no_neb']} ({100*stats['no_neb']/max(stats['total'],1):.1f}%)")
    print("="*70)
    
    # Print legend
    print("\nSTATUS LEGEND:")
    print("  \033[0;32mCOMPLETED_OK\033[0m          - All calculations completed successfully")
    print("  \033[0;37mNOT_STARTED\033[0m           - Reoptimization not started yet")
    print("  \033[1;33mTS_OPT_RUNNING\033[0m        - TS optimization in progress")
    print("  \033[1;33mFREQ_RUNNING\033[0m          - Frequency calculations in progress")
    print("  \033[1;33mIN_PROGRESS\033[0m           - Other calculations in progress")
    print("  \033[0;35mCOMPLETED_START_IMAG\033[0m  - Completed but start has imaginary modes")
    print("  \033[0;31mCOMPLETED_TS_IMAG_ERROR\033[0m - Completed but TS has wrong number of imaginary modes")
    print("  \033[0;31mNO_NEB\033[0m                - NEB data not found")
    print()


if __name__ == '__main__':
    main()
