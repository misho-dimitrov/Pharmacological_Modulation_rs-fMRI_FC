import sys
import os
import nibabel as nib


# Load an image from a single brain:
def img_load(folder, img):
    """Load an image from a single brain.
    
    The function returns a brain image as a nibabel image object.
    
    Args:
        folder (string): Image directory.
        img (string): Image filename.
        
    Returns:
        fmri_image (nibabel.nifti1.Nifti1Image:): Brain image.
        
    """
    
    img_path = os.path.join(folder, img)
    fmri_image = nib.load(img_path)
    print(fmri_image.shape)
    return fmri_image

if __name__ == "__main__":
    folder = sys.argv[1]
    img = sys.argv[2]
    img_load(folder, img)
