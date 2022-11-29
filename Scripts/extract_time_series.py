import sys
from nilearn.input_data import NiftiMasker

# Extract the timeseries from a single session
def extract_time_series(fmri_image, mask):
    brain_masker = NiftiMasker(mask, memory='nilearn_cache', memory_level=1, verbose=0)
    brain_time_series = brain_masker.fit_transform(fmri_image)
    return brain_time_series, brain_masker

if __name__ == "__main__":
    fmri_image = sys.argv[1]
    mask = sys.argv[2]
    extract_time_series(fmri_image, mask)


