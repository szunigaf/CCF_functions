import sys
sys.path.append("./MODULES")

import numpy as np
from calc_rv_new_chapman import master_make_ccf_vsini
import glob
import argparse
import os
import configparser

def parse_arguments():
    """Parse command line arguments for CCF analysis parameters"""
    parser = argparse.ArgumentParser(description='CCF Analysis for High-Resolution Spectra')
    
    # Configuration file option
    parser.add_argument('--config', type=str, 
                       help='Configuration file (INI format)')
    
    # Basic parameters
    parser.add_argument('--data_dir', type=str, default='./INPUT/', 
                       help='Directory containing FITS files (default: ./INPUT/)')
    parser.add_argument('--mask', type=str, default='K0', 
                       choices=['K0', 'K5', 'M0', 'M2', 'M4', 'M5', 'R50G2', 'G2', 'F0','M2_IR','M_5_IR'],
                       help='Mask type for CCF analysis (default: K0)')
    parser.add_argument('--instrument', type=str, default='HARPS',
                       choices=['HARPS', 'NIRPS', 'UVES', 'UVES_combine', 'FEROS', 'XSHOO', 'ELODIE'],
                       help='Instrument type (default: HARPS)')
    parser.add_argument('--output_dir', type=str, default='./OUTPUTS/',
                       help='Output directory for results (default: ./OUTPUTS/)')
    
    # UVES_combine specific arguments
    parser.add_argument('--file1', type=str,
                       help='First FITS file for UVES_combine (e.g., blue arm)')
    parser.add_argument('--file2', type=str,
                       help='Second FITS file for UVES_combine (e.g., red arm)')
    parser.add_argument('--pairs_file', type=str,
                       help='Text file containing pairs of files for batch UVES_combine processing')
    
    # Advanced CCF parameters
    parser.add_argument('--vel_step', type=float, default=1.0,
                       help='Velocity step size in km/s (default: 1.0)')
    parser.add_argument('--vel_min', type=float, default=-150.0,
                       help='Minimum velocity range in km/s (default: -150)')
    parser.add_argument('--vel_max', type=float, default=150.0,
                       help='Maximum velocity range in km/s (default: 150)')
    parser.add_argument('--vsini_step', type=float, default=1.0,
                       help='VSini step for profile fitting in km/s (default: 1.0)')
    parser.add_argument('--wlen_step', type=float, default=0.02,
                       help='Wavelength interpolation step in Angstroms (default: 0.02)')
    parser.add_argument('--ncpu', type=int, default=8,
                       help='Number of CPU cores to use (default: 8)')
    
    # Processing options
    parser.add_argument('--file_pattern', type=str, default='*fits',
                       help='File pattern to search for (default: *fits)')
    parser.add_argument('--resume', action='store_true',
                       help='Resume processing (skip already processed files)')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose output')
    parser.add_argument('--plot_flag', action='store_true',
                       help='Display intermediate plots during processing')
    # ... rest of arguments ...
    
    args = parser.parse_args()
    
    # Validation for UVES_combine
    if args.instrument == 'UVES_combine':
        if not args.pairs_file and not (args.file1 and args.file2):
            parser.error("UVES_combine requires either --pairs_file or both --file1 and --file2")
        if args.file1 and args.file2 and args.pairs_file:
            parser.error("Use either --pairs_file OR individual files (--file1, --file2), not both")
    
    # Load configuration file if specified
    if args.config:
        load_config(args)
    
    return args

def load_config(args):
    """Load parameters from configuration file"""
    if not os.path.exists(args.config):
        raise ValueError(f"Configuration file not found: {args.config}")
    
    config = configparser.ConfigParser()
    config.read(args.config)
    
    # Update args with config file values (command line args take precedence)
    if 'Basic' in config:
        for key, value in config['Basic'].items():
            if not hasattr(args, key) or getattr(args, key) == argparse.ArgumentParser().get_default(key):
                if key == 'resume' or key == 'verbose':
                    setattr(args, key, config.getboolean('Basic', key))
                else:
                    setattr(args, key, value)
    
    if 'CCF_Parameters' in config:
        for key, value in config['CCF_Parameters'].items():
            if not hasattr(args, key) or getattr(args, key) == argparse.ArgumentParser().get_default(key):
                if key in ['vel_step', 'vel_min', 'vel_max', 'vsini_step', 'wlen_step']:
                    setattr(args, key, config.getfloat('CCF_Parameters', key))
                elif key == 'ncpu':
                    setattr(args, key, config.getint('CCF_Parameters', key))
    
    if 'Processing' in config:
        for key, value in config['Processing'].items():
            if not hasattr(args, key) or getattr(args, key) == argparse.ArgumentParser().get_default(key):
                if key == 'resume' or key == 'verbose':
                    setattr(args, key, config.getboolean('Processing', key))
                else:
                    setattr(args, key, value)

def validate_parameters(args):
    """Validate input parameters and show warnings"""
    warnings = []
    
    # Check data directory
    if not os.path.exists(args.data_dir):
        raise ValueError(f"Data directory does not exist: {args.data_dir}")
    
    # Check output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Validate velocity range
    if args.vel_min >= args.vel_max:
        raise ValueError("vel_min must be less than vel_max")
    
    # Performance warnings
    vel_range = args.vel_max - args.vel_min
    n_vel_points = int(vel_range / args.vel_step)
    
    if args.vel_step < 0.5:
        warnings.append(f"Small vel_step ({args.vel_step}) will increase processing time significantly")
    
    if n_vel_points > 500:
        warnings.append(f"Large velocity range ({n_vel_points} points) may slow processing")
    
    if args.wlen_step < 0.01:
        warnings.append(f"Small wlen_step ({args.wlen_step}) will significantly increase processing time")
    
    if args.vsini_step < 1:
        warnings.append(f"Small vsini_step ({args.vsini_step}) will increase processing time")
    
    if args.ncpu > 16:
        warnings.append(f"High CPU count ({args.ncpu}) may not improve performance")
    
    return warnings

def print_configuration(args):
    """Print current configuration"""
    print("\n" + "="*60)
    print("CCF ANALYSIS CONFIGURATION")
    print("="*60)
    print(f"Data Directory:     {args.data_dir}")
    print(f"Output Directory:   {args.output_dir}")
    print(f"Mask Type:          {args.mask}")
    print(f"Instrument:         {args.instrument}")
    print(f"File Pattern:       {args.file_pattern}")
    print(f"Resume Mode:        {args.resume}")
    print("\nCCF Parameters:")
    print(f"  Velocity Range:   {args.vel_min} to {args.vel_max} km/s")
    print(f"  Velocity Step:    {args.vel_step} km/s")
    print(f"  VSini Step:       {args.vsini_step} km/s")
    print(f"  Wavelength Step:  {args.wlen_step} Å")
    print(f"  CPU Cores:        {args.ncpu}")
    
    # Calculate processing points
    n_vel = int((args.vel_max - args.vel_min) / args.vel_step)
    print(f"\nProcessing Points:  {n_vel} velocity steps")
    print("="*60)

def get_files_to_process(args):
    """Get list of files to process - handles manual UVES_combine specification"""
    
    if args.instrument == 'UVES_combine':
        files_to_process = []
        
        if args.file1 and args.file2:
            # Single pair mode
            if not os.path.exists(args.file1):
                raise ValueError(f"File1 not found: {args.file1}")
            if not os.path.exists(args.file2):
                raise ValueError(f"File2 not found: {args.file2}")
            
            files_to_process = [[args.file1, args.file2]]
            all_files = [args.file1, args.file2]
            
            if args.verbose:
                print(f"UVES_combine pair:")
                print(f"  File1: {args.file1}")
                print(f"  File2: {args.file2}")
        
        elif args.pairs_file:
            # Batch mode from pairs file
            if not os.path.exists(args.pairs_file):
                raise ValueError(f"Pairs file not found: {args.pairs_file}")
            
            with open(args.pairs_file, 'r') as f:
                lines = f.readlines()
            
            all_files = []
            for line_num, line in enumerate(lines, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split()
                if len(parts) != 2:
                    print(f"Warning: Invalid format in line {line_num}: {line}")
                    continue
                
                file1, file2 = parts
                
                # Convert relative paths to absolute if needed
                if not os.path.isabs(file1):
                    file1 = os.path.join(args.data_dir, file1)
                if not os.path.isabs(file2):
                    file2 = os.path.join(args.data_dir, file2)
                
                if not os.path.exists(file1):
                    print(f"Warning: File1 not found in line {line_num}: {file1}")
                    continue
                if not os.path.exists(file2):
                    print(f"Warning: File2 not found in line {line_num}: {file2}")
                    continue
                
                files_to_process.append([file1, file2])
                all_files.extend([file1, file2])
                
                if args.verbose:
                    print(f"Pair {len(files_to_process)}: {os.path.basename(file1)} + {os.path.basename(file2)}")
    
    else:
        # Normal single file processing for other instruments
        search_pattern = os.path.join(args.data_dir, args.file_pattern)
        data_files = glob.glob(search_pattern)
        
        if not data_files:
            raise ValueError(f"No files found matching pattern: {search_pattern}")
        
        files_to_process = data_files
        all_files = data_files
        
        if args.verbose:
            print(f"\nFound {len(data_files)} files:")
            for fname in [os.path.basename(f) for f in data_files[:10]]:  # Show first 10
                print(f"  {fname}")
            if len(data_files) > 10:
                print(f"  ... and {len(data_files) - 10} more files")
    
    # Handle resume mode for UVES_combine
    if args.resume and args.instrument == 'UVES_combine':
        files_to_process = handle_uves_resume(files_to_process, args)
    elif args.resume:
        files_to_process = handle_normal_resume(files_to_process, args)
    
    return files_to_process, all_files

def handle_uves_resume(files_to_process, args):
    """Handle resume mode for UVES_combine processing"""
    master_output = os.path.join(args.output_dir, "master_output.txt")
    
    try:
        existing_data = np.loadtxt(master_output, 
                                 dtype={'names':('file','obj_name','ra_deg','dec_deg',
                                               'mjd_obs','rv','width','depth','vsini_calc_val',
                                               'min_chi','snr_median','mask','central_wlen',
                                               'wlen_range','BIS','b_b', 'b_b_err','c_b','ad_sig'), 
                                      'formats':('U100', 'U50', 'f8', 'f8', 'f8', 'f8', 
                                               'f8', 'f8', 'f8', 'f8', 'f8', 'U2', 'f8', 
                                               'f8', 'f8', 'f8', 'f8', 'f8', 'f8')}, 
                                 skiprows=1, unpack=True)
        
        existing_files = set(existing_data[0])
        
        # Filter out already processed pairs
        filtered_pairs = []
        for file_pair in files_to_process:
            file1, file2 = file_pair
            # Create identifier for this pair
            pair_id = f"{os.path.basename(file1)}+{os.path.basename(file2)}"
            
            if pair_id not in existing_files:
                filtered_pairs.append(file_pair)
                
        if args.verbose:
            processed_count = len(files_to_process) - len(filtered_pairs)
            print(f"\nResume mode: {processed_count} pairs already processed")
            
        return filtered_pairs
        
    except (OSError, IOError, ValueError):
        print("No existing output file found - processing all pairs")
        return files_to_process


def handle_normal_resume(files_to_process, args):
    """Handle resume mode for single-file instruments"""
    master_output = os.path.join(args.output_dir, "master_output.txt")
    if not os.path.exists(master_output):
        print("No existing output file found - processing all files")
        return files_to_process
    try:
        with open(master_output, 'r') as f:
            lines = f.readlines()[1:]  # skip header
        existing_files = set()
        for line in lines:
            if not line.strip():
                continue
            first_token = line.split()[0]
            existing_files.add(first_token)
        filtered = []
        for fpath in files_to_process:
            base = os.path.basename(fpath)
            if base not in existing_files:
                filtered.append(fpath)
        if args.verbose:
            processed_count = len(files_to_process) - len(filtered)
            print(f"\nResume mode: {processed_count} files already processed")
        return filtered
    except Exception:
        print("Could not parse existing output file - processing all files")
        return files_to_process

def process_files(files_to_process, args):
    """Process the files with CCF analysis - handles manual UVES_combine pairs"""
    
    master_output = os.path.join(args.output_dir, "master_output.txt")
    
    # Create header if file doesn't exist
    if not os.path.exists(master_output):
        header = "file     obj_name      ra_deg    dec_deg    mjd_obs    rv    width    depth    vsini_calc_val      min_chi    snr_median    mask      central_wlen    wlen_range    BIS    b_b    b_b_err    c_b    ad_sig\n"
        with open(master_output, "w") as f:
            f.write(header)
    
    total_files = len(files_to_process)
    file_type = "pairs" if args.instrument == 'UVES_combine' else "files"
    print(f"\nProcessing {total_files} {file_type}...")
    print("-" * 60)
    
    for i, file_input in enumerate(files_to_process):
        
        if args.instrument == 'UVES_combine':
            # file_input is a list of two files
            file1, file2 = file_input
            filename = f"{os.path.basename(file1)}+{os.path.basename(file2)}"
            
            if args.verbose:
                print(f"[{i+1}/{total_files}] Processing UVES pair:")
                print(f"  File1: {os.path.basename(file1)}")
                print(f"  File2: {os.path.basename(file2)}")
            else:
                print(f"[{i+1}/{total_files}] {os.path.basename(file1)} + {os.path.basename(file2)}")
        else:
            # Single file processing
            file1 = file_input
            filename = os.path.basename(file_input)
            
            if args.verbose:
                print(f"[{i+1}/{total_files}] Processing: {filename}")
            else:
                print(f"[{i+1}/{total_files}] {filename}")
        
        try:
            # Call the analysis function
            result = master_make_ccf_vsini(
                file_input if args.instrument == 'UVES_combine' else file1,
                flag='fits', 
                mask_type=args.mask, 
                instrument=args.instrument,
                vel_step=args.vel_step,
                vel_range_min=args.vel_min,
                vel_range_max=args.vel_max,
                vsini_step=args.vsini_step,
                wlen_step=args.wlen_step,
                ncpu=args.ncpu
            )
            
            # Unpack results
            obj_name, ra_deg, dec_deg, mjd_obs, rv, width, depth, \
            vsini_calc_val, min_chi, snr_median, central_wlen, \
            wlen_range, BIS, b_b, stderr, c_b, ad_sig = result
            
            # Format output
            output_string = f"{filename}     {obj_name}      {ra_deg:.5f}    {dec_deg:.5f}    {mjd_obs:.3f}    {rv:.2f}    {width:.2f}    {depth:.3f}    {vsini_calc_val:d}      {min_chi:.3f}    {snr_median:.2f}    {args.mask}      {central_wlen:.1f}    {wlen_range:.1f}    {BIS:.3f}    {b_b:.3f}    {stderr:.3f}    {c_b:.3f}    {ad_sig:.2f}\n"
            
            # Write to output file
            with open(master_output, "a") as text_file:
                text_file.write(output_string)
                
            if args.verbose:
                print(f"  -> RV: {rv:.2f} km/s, VSini: {vsini_calc_val} km/s, S/N: {snr_median:.1f}")
                
        except Exception as e:
            print(f"  ERROR processing {filename}: {str(e)}")
            continue
    
    print("-" * 60)
    print(f"Processing complete! Results saved to: {master_output}")

def main():
    """Main function"""
    try:
        # Parse arguments
        args = parse_arguments()
        
        # Validate parameters and show warnings
        warnings = validate_parameters(args)
        
        # Print configuration
        print_configuration(args)
        
        # Show warnings if any
        if warnings:
            print("\n⚠️  PERFORMANCE WARNINGS:")
            for warning in warnings:
                print(f"  • {warning}")
            print()
            
        # Get files to process
        files_to_process, all_files = get_files_to_process(args)
        
        if not files_to_process:
            print("No files to process!")
            return
            
        print(f"\n{len(files_to_process)} / {len(all_files)} files need processing")
        
        # Ask for confirmation if processing time might be long
        total_points = int((args.vel_max - args.vel_min) / args.vel_step)
        if total_points > 300 or args.wlen_step < 0.015:
            response = input("\nProcessing may take a long time. Continue? (y/N): ")
            if response.lower() not in ['y', 'yes']:
                print("Processing cancelled.")
                return
        
        # Process files
        process_files(files_to_process, args)
        
    except KeyboardInterrupt:
        print("\n\nProcessing interrupted by user.")
    except Exception as e:
        print(f"\nError: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
