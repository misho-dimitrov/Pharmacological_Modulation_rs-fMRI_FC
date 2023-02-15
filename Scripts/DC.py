import csv
import sys
from img_load import img_load
from mask_img_load import mask_img_load
from extract_time_series import extract_time_series
from cos_sim_func import cos_sim_func
from thresh_mat import thresh_mat
from calc_dc import calc_dc
#from z_score_dc import z_score_dc
from array_to_nifti import array_to_nifti
import nibabel as nib
import pickle


# To calculate functional connectivity and degree centrality from an adjacency matrix
def calc_dc_auto(fMRIroot, fMRI_txt, MASKroot, GM_mask):
    """ Calculate functional connectivity and degree centrality from an adjacency matrix and save the DC images to disk.
    
    Args:
        fMRIroot (string): fMRI image directory.
        fMRI_txt (string): A txt file containing the fMRI image filenames.
        MASKroot (string): Mask image directory.
        GM_mask (string): Mask image filename.
        
    """
    
    # transform the txt file into a Python list... of lists!
    fmri_list = []
    with open(fMRI_txt, newline='') as inputfile:
        for row in csv.reader(inputfile):
            fmri_list.append(row)
    
    for i in range(len(fmri_list)):
        # load an image from a single brain:
        fmri_image = img_load(fMRIroot, fmri_list[i][0])
        print("Participant " + fmri_list[i][0])
        print("The image dimensions are: " + str(fmri_image.shape))
        
        # load the GM mask
        gm_mask = mask_img_load(MASKroot, GM_mask)
        print("GM mask loaded.")
        
        # create a NiftiMasker object and extract the time series
        brain_time_series, brain_masker = extract_time_series(fmri_image, gm_mask)
        print("Time series extracted.")
        
        # calculate the correlations between each pair of voxels
        cos_sim = cos_sim_func(brain_time_series)
        print("Correlations calculated.")
        
        # threshold the matrix
        adjacency_matrix = thresh_mat(cos_sim)
        print("Adjacency matrix thresholded.")
        
        # calculate the degree centrality (DC)
        dc = calc_dc(adjacency_matrix)
        print("Raw wDC values calcualted.")
        
        # convert the raw DC matrix and save it as a NIfTI image
        dc_img = array_to_nifti(brain_masker,dc)
        dc_img_name = fmri_list[i][0][:6] + "_wDC_raw"
        nib.save(dc_img, dc_img_name)
        print("Raw wDC image of the whole GM of " + fmri_list[i][0][:6] + " saved.")
        
    print("Done!")

if __name__ == "__main__":
    fMRIroot = sys.argv[1]
    fMRI_txt = sys.argv[2]
    MASKroot = sys.argv[3]
    GM_mask = sys.argv[4]

    calc_dc_auto(fMRIroot, fMRI_txt, MASKroot, GM_mask)


