import sys
import os
import nibabel as nib

# Load my resampled GM mask
def mask_img_load(folder, img):
    """Load a brain mask image.
    
    The function returns a brain mask image as a nibabel image object.
    
    Args:
        folder (string): Mask image directory.
        img (string): Mask image filename.
        
    Returns:
        mask (nibabel.nifti1.Nifti1Image:): Mask image.
        
    """
    
    img_path = os.path.join(folder, img)
    mask = nib.load(img_path)
    return mask

if __name__ == "__main__":
    folder = sys.argv[1]
    img = sys.argv[2]
    mask_img_load(folder, img)



