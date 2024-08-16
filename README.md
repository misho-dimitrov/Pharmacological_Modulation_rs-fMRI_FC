# Pharmacological Modulation of Functional Brain Connectivity in Autism

## Estimate resting-state functional brain connectivity (FC) from fMRI images and extract graph features. <br /> Assess the effect of various drugs on FC and graph properties. <br /> Pinpoint potential pharmacological mechanisms of these drugs using the REACT method described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6547164/.  

### Main Scripts and Commands: 

<ins>Pre-processing:</ins> <br />
A) Extract the DICOM images from the TAR archive and convert them to NIfTI (using SGE) - *untar_dcm2nii_final.csh* <br />
B) Generate the script for the AFNI-based pre-processing (https://pubmed.ncbi.nlm.nih.gov/8812068/) - *afni_multi-echo.rtf* (the output file, i.e. proc.SUBJ, is then run using *tcsh -xef proc.SUBJ | \& tee output.proc.SUBJ*) <br />

<ins>Functional Connectivity Estimation:</ins> <br />
C) Compute FC and (static) weighted degree centrality (wDC) - *voxel_wise_FC_calc.ipynb* and *z-score.ipynb* <br />

<ins>Post-processing:</ins> <br />
D) Resample and binarise the Holiga et al. EU-AIMS masks (https://www.science.org/doi/10.1126/scitranslmed.aat9223), and intersect each with the intersection mask generated using the grey matter masks from all participants (all sessions) - *generate_intersection_mask.ipynb*; Estimate mask-averaged values for each mask for each participant session - *generate_mean_DC.ipynb* <br />
E) Compute the REACT (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6547164/) maps (after masking out the reference region in the PET image, if necessary) - using the instructions on https://github.com/ottaviadipasquale/react-fmri <br />

<ins>Analysis and plotting:</ins> <br />
F) Check if the study groups are balanced in terms of age, IQ and in-scanner movement as well as depression and anxiety scores - *balance_check.ipynb* and *HAM-DA_check.R*<br />

G) Mask-averaged/mean wDC analyses - *mean_wDC_LMM.R* <br />
H) Individual mean wDC trajectories - *spaghetti_individual_trajectories.R* <br />

I) Voxel-wise wDC analyses
        i.) Main effects of group and drug as well as their interaction - *voxel-wise_wDC_LMM_Main.R* <br />
        ii.) Within-group drug effects - *voxel-wise_wDC_LMM_Within.R* <br />
J) Visualisation of voxel-wise wDC results - *plot_final_results.ipynb* <br />

K) Voxel-wise REACT analysis - Randomise (FSL - https://pubmed.ncbi.nlm.nih.gov/15501092/) <br />
        i.) Within-group comparisons: *randomise -i study_group_wDC_img.nii.gz -o drug_effect_study_group -m intersection_mask.nii.gz -1 -v 5 -T* (results visualised in FSLeyes, FSL) <br />

<ins>*Several fslmaths (FSL - https://pubmed.ncbi.nlm.nih.gov/15501092/) commands have been used to additionally manipulate images, e.g. to threshold, binarise or resample images.* </ins>

<ins>*N.B. conda environment file: Project_3_TIA.yml*</ins>

<ins>Please check my PhD thesis (***Link to be added..***) for more details.</ins>