import sys
import os
import nibabel as nib

# Load my resampled GM mask
def mask_img_load(folder, img):
    img_path = os.path.join(folder, img)
    mask = nib.load(img_path)
    return mask

if __name__ == "__main__":
    folder = sys.argv[1]
    img = sys.argv[2]
    mask_img_load(folder, img)



