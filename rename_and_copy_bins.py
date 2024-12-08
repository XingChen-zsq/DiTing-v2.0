import os
import sys
import shutil
from glob import glob

def rename_and_copy_bins(refinement_dir, output_dir, basename):
    """
    Rename and copy .fa files from the refinement directory to the output directory.

    Parameters:
    refinement_dir (str): Path to the refinement directory containing metawrap_50_10_bins.
    output_dir (str): Path to the output directory where renamed files will be saved.
    basename (str): Sample basename to prefix file names.
    """
    # Path to the folder containing .fa files
    bins_dir = os.path.join(refinement_dir, "metawrap_50_10_bins")
    
    # Check if the directory exists
    if not os.path.exists(bins_dir):
        print(f"Warning: Directory {bins_dir} does not exist. Skipping sample {basename}.")
        return

    # Get all .fa files in the directory
    fa_files = sorted(glob(os.path.join(bins_dir, "*.fa")))

    # Check if there are any .fa files
    if not fa_files:
        print(f"Warning: No .fa files found in {bins_dir}. Skipping sample {basename}.")
        return

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Rename and copy each .fa file
    for fa_file in fa_files:
        original_name = os.path.basename(fa_file)
        new_name = f"{basename}_{original_name}"
        new_path = os.path.join(output_dir, new_name)
        
        # Copy the file with the new name
        shutil.copy(fa_file, new_path)
        print(f"Copied and renamed: {fa_file} -> {new_path}")

if __name__ == "__main__":
    # Check command-line arguments
    if len(sys.argv) != 4:
        print("Usage: python rename_and_copy_bins.py <refinement_dir> <output_dir> <basename>")
        sys.exit(1)

    # Parse arguments
    refinement_dir = sys.argv[1]
    output_dir = sys.argv[2]
    basename = sys.argv[3]

    # Call the function
    rename_and_copy_bins(refinement_dir, output_dir, basename)

