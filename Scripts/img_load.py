import sys
import os
import nibabel as nib


# Load an image from a single brain:
def img_load(folder, img):
    img_path = os.path.join(folder, img)
    fmri_image = nib.load(img_path)
    print(fmri_image.shape)
    return fmri_image

if __name__ == "__main__":
    folder = sys.argv[1]
    img = sys.argv[2]
    img_load(folder, img)
