#!/usr/bin/env python3

import nibabel as nib
import numpy as np
from scipy.ndimage import zoom
import argparse
import os

def transform_to_freesurfer_space(in_file, out_file):
    """
    Transform T1 MPRAGE to FreeSurfer's standard space (256x256x256, 1mm isotropic)
    and save as NIfTI
    
    Parameters:
    -----------
    in_file : str
        Path to input T1 MPRAGE image
    out_file : str
        Path where transformed image will be saved (as NIfTI)
    """
    # Load the input image
    img = nib.load(in_file)
    data = img.get_fdata()
    
    # Get current dimensions and voxel sizes
    current_shape = data.shape
    voxel_sizes = img.header.get_zooms()
    
    # Calculate zoom factors to achieve 1mm isotropic in 256x256x256
    target_shape = (256, 256, 256)
    zoom_factors = [t/c for t, c in zip(target_shape, current_shape)]
    
    # Account for current voxel sizes in zoom factors
    zoom_factors = [z*v for z, v in zip(zoom_factors, voxel_sizes)]
    
    # Resample the image
    resampled_data = zoom(data, zoom_factors, order=3)  # order=3 for cubic interpolation
    
    # Create new affine matrix for 1mm isotropic voxels
    new_affine = np.eye(4)
    new_affine[0:3, 0:3] = np.diag([1.0, 1.0, 1.0])  # 1mm isotropic voxels
    
    # Center the image in the new space
    origin = np.array(target_shape) / 2
    new_affine[0:3, 3] = -origin
    
    # Create new NIfTI image
    new_img = nib.Nifti1Image(resampled_data, new_affine)
    
    # Set intent code to NONE (consistent with FreeSurfer's orig.mgz)
    new_img.header.set_intent(code=0)
    
    # Update header information
    new_img.header['pixdim'][1:4] = [1.0, 1.0, 1.0]  # Set voxel dimensions
    new_img.header['qform_code'] = 1  # Use scanner-based coordinate system
    new_img.header['sform_code'] = 1
    
    # Save the transformed image
    nib.save(new_img, out_file)

def main():
    parser = argparse.ArgumentParser(description='Transform T1 MPRAGE to FreeSurfer standard space')
    parser.add_argument('in_file', help='Path to input T1 MPRAGE image')
    parser.add_argument('out_file', help='Path where transformed image will be saved')
    args = parser.parse_args()
    
    # Ensure output file has .nii.gz extension
    if not args.out_file.endswith(('.nii', '.nii.gz')):
        base, _ = os.path.splitext(args.out_file)
        args.out_file = base + '.nii.gz'
    
    transform_to_freesurfer_space(args.in_file, args.out_file)

if __name__ == '__main__':
    main()