# Pharmacological Modulation of Functional Brain Connectivity in Autism

## Calculate resting-state functional brain connectivity (FC) from fMRI images and extract graph features. <br /> Assess the effect of various drugs on FC and graph properties. <br /> Pinpoint potential pharmacological mechanisms of these drugs using the REACT method described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6547164/. <br /> Investigate the ability of baseline FC to predict drug response. 

### Main Scripts and Commands: 

<ins>Pre-processing:</ins> <br />
A) Extract the DICOM images from the TAR archive and convert them to NIfTI (using SGE) - *untar_dcm2nii_final.csh* <br />
B) Generate the script for the AFNI-based pre-processing (https://pubmed.ncbi.nlm.nih.gov/8812068/) - *afni_single-echo.rtf* OR *afni_multi-echo.rtf* (the output file, i.e. proc.SUBJ, is then run using *tcsh -xef proc.SUBJ | \& tee output.proc.SUBJ*) <br />

<ins>Functional Connectivity Estimation:</ins> <br />
C) Calculate the FC and the (static) weighted degree centrality (wDC) - *voxel_wise_weighted_DC.ipynb* and *z-score.ipynb* <br />

<ins>Post-processing:</ins> <br />
D) Check if the study groups are balanced - *possible_covariate_check.ipynb* <br />
E) Calculate the Holiga et al. EU-AIMS (https://www.science.org/doi/10.1126/scitranslmed.aat9223) mask-averaged wDC and intersect these masks with the intersection mask generated using all participant images - *generate_mean_DC.ipynb* and *generate_intersection_mask.ipynb* <br />
F) Calculate the REACT (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6547164/) maps (after masking out the reference region in the PET image, if necessary) - using the instructions on https://github.com/ottaviadipasquale/react-fmri <br />

<ins>Analysis and plotting:</ins> <br />
G) Between-group comparison (Holiga et al. masks) - *mean_dc_group_comparison_two_sample.R* <br />
H) Within-group comparison (Holiga et al. masks) - *mean_dc_group_comparison_paired_sample.R* <br />
I) Three-group comparison (Holiga et al. masks) - *mean_dc_group_comparison_lme.R* <br />
J) Voxel-wise group comparisons (both wDC and REACT) - *generate_whole_group_images.ipynb* + Randomise (FSL - https://pubmed.ncbi.nlm.nih.gov/15501092/) <br />
        i.) Between-group comparison: *randomise -i between_group_wDC_img.nii -o group_comparison -d design.mat -t design.con -m intersection_mask.nii.gz -n 5000 -D -T* <br />
        ii.) Within-group comparison: *randomise -i study_group_wDC_img.nii.gz -o drug_effect_study_group -m intersection_mask.nii.gz -1 -v 5 -T* <br />
        iii.) Plotting: *plot_final_results_lateral_and_medial.ipynb* <br />
K) Individual trajectories - *spaghetti_individual_trajectories.R* OR *spaghetti_individual_trajectories_3_conditions.R* <br />
L) Clustering of drug responders and non-responders - ***GitHub repo link to be added..*** <br />
M) Additional correlations - *correlation.R* <br />

*Several fslmaths (FSL - https://pubmed.ncbi.nlm.nih.gov/15501092/) commands have been used to additionally manipulate images, e.g. to threshold, binarise or resample images. Please check my PhD thesis (***Link to be added..***) for more details.* 

**N.B. conda environment file: Project_3_TIA.yml**