import sys
import numpy
from nilearn.input_data import NiftiMasker

# Transform the FC or DC array back into a NIfTI volume and save
def array_to_nifti(brain_masker,array):
    img = brain_masker.inverse_transform(array.T)
    return img

if __name__ == "__main__":
    brain_masker = sys.argv[1]
    array = sys.argv[2]
    array_to_nifti(array)


