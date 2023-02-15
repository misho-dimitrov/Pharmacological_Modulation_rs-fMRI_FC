import sys
import numpy
from nilearn.input_data import NiftiMasker

# Transform the FC or DC array back into a NIfTI volume and save
def array_to_nifti(brain_masker,array):
    """ Transform the degree centrality array into a NIfTI volume and save.

    The function returns a degree centrality image.
    
    Args:
        brain_masker (NiftiMasker): A NiftiMasker object that could be used to extract timeseries from fMRI images or transform 2D timeseries arrays back into 4D fMRI images.
        dc_array (ndarray): A 1D voxel-wise degree centrality array.
        
    Returns:
        img (nibabel.nifti1.Nifti1Image:): Degree centrality image.
        
    """
    
    img = brain_masker.inverse_transform(array.T)
    return img

if __name__ == "__main__":
    brain_masker = sys.argv[1]
    array = sys.argv[2]
    array_to_nifti(array)


