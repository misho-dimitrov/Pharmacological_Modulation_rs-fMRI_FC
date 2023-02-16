#!/bin/csh

# Set up the sge queue required defaults
#$ -M mihail.dimitrov@kcl.ac.uk
#$ -m ae
#$ -o /data/project/TRADA/nii/o_logs/untar_dcm2nii/
#$ -e /data/project/TRADA/nii/e_logs/untar_dcm2nii/
#$ -q global
#$ -N trada

# Load the required environment defaults
module load afni/21.1.07

# Set the working variables 
set study_name=/data/project/TRADA/nii

# list of subject files
set sge_index=${study_name}/sge_index

# Search the file for the SGE_TASK_ID number as a line number
set file="`awk 'FNR==$SGE_TASK_ID' ${sge_index}`"
#cd $file/DICOM/??????/
cd $file
pwd

# Set each archive to the correct 'type'
set fmri=`\zsh -c 'print -r ./*(DoL[-1])'`
set fmri_dcm=`\echo $fmri | awk -F'[.]' '{print $3}'`
set t1=`\zsh -c 'print -r ./*(DoL[-2])'`
set t1_dcm=`\echo $t1 | awk -F'[.]' '{print $3}'`
#set t2=`\zsh -c 'print -r ./*(DoL[-3])'`
#set t2_dcm=`\echo $t2 | awk -F'[.]' '{print $3}'`

#echo $fmri
#echo $fmri_dcm
#echo $t1
#echo $t1_dcm
echo "Variable names set."

# Untar and dcm2niix the files
#fMRI
tar -xjf $fmri
echo "Untaring of fMRI completed."
dcm2niix -b y -f "f_e%e" -o ./ -v y $fmri_dcm
#cp medicom_c ./$fmri_dcm
#echo "$fmri_dcm"
#cd $fmri_dcm
#./medicom_c fmri
echo "Conversion of fMRI completed."
cd ./..
#T1
tar -xjf $t1
echo "Untaring of T1 completed."
dcm2nii -g n -x y $t1_dcm
echo "Conversion of T1 completed."
#T2
#tar -xjf $t2
#echo "Untaring of T2 completed."
#dcm2nii -g n $t2_dcm
#echo "Conversion of T2 completed."

echo "Done."

