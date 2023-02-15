import sys
from nilearn.input_data import NiftiMasker

# Extract the timeseries from a single session
def extract_time_series(fmri_image, mask):
    """Extract the timeseries from a single brain (fMRI image).
    
    The function returns the timeseries data from each brain voxel.
    
    Args:
        fmri_image (string): Brain image (fMRI).
        mask (string): Mask image.
        
    Returns:
        brain_time_series (ndarray): A timeseries array in the format [n_volumes, n_voxels].
        brain_masker (NiftiMasker object): A NiftiMasker object that could be used to extract timeseries from fMRI images or transform 2D timeseries arrays back into 4D fMRI images.
        
    """
    
    brain_masker = NiftiMasker(mask, memory='nilearn_cache', memory_level=1, verbose=0)
    brain_time_series = brain_masker.fit_transform(fmri_image)
    return brain_time_series, brain_masker

if __name__ == "__main__":
    fmri_image = sys.argv[1]
    mask = sys.argv[2]
    extract_time_series(fmri_image, mask)


