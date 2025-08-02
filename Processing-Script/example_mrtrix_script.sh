#-------------------Step 1: Preparing data--------------------#


#-------------------------Dictionary setup------------------------------------------#
SUBJECT = 100206

module load mrtrix3
module load ants
module load fsl

MainDIR=/xxxxx/SC_Processing
PrepareDIR=${MainDIR}/File_for_processing
INPUTDIR=${MainDIR}/Data_ready_for_tractography

diffdir=${INPUTDIR}/${SUBJECT}
regimg=${diffdir}/nodif_brain_mask.nii.gz

SCRATCHDIR=${MainDIR}/scratch/


#-------------------Step 2: Preprocessing before tractography--------------------#
# compute CSD, bias correction, 5 tissue model (by FSL FAST)
###Convert data and compute CSD
mrconvert ${diffdir}/data.nii.gz ${diffdir}/DWI.mif -fslgrad ${diffdir}/bvecs ${diffdir}/bvals -datatype float32 -stride 0,0,0,1 -quiet -force
dwibiascorrect ants ${diffdir}/DWI.mif ${diffdir}/DWI_bc.mif -bias ${diffdir}/DWI_biasfield.mif -mask ${diffdir}/nodif_brain_mask.nii.gz -scratch ${SCRATCHDIR} -force
dwi2response dhollander ${diffdir}/DWI_bc.mif ${diffdir}/RF_wm_dhollander.txt ${diffdir}/RF_gm_dhollander.txt ${diffdir}/RF_csf_dhollander.txt -voxels ${diffdir}/RF_voxels_dhollander.mif -scratch ${SCRATCHDIR} -force
dwi2fod msmt_csd ${diffdir}/DWI_bc.mif ${diffdir}/RF_wm_dhollander.txt ${diffdir}/RF_wm_dhollander.mif ${diffdir}/RF_gm_dhollander.txt ${diffdir}/RF_gm_dhollander.mif ${diffdir}/RF_csf_dhollander.txt ${diffdir}/RF_csf_dhollander.mif -force


#Make 5TT tissue model using FSL FAST, resampled to diffusion volume space
5ttgen fsl ${diffdir}/T1w_acpc_dc_restore.nii.gz ${diffdir}/5TT_fsl_T1space.mif -mask ${diffdir}/T1w_acpc_dc_restore_brain.nii.gz  -scratch ${SCRATCHDIR} -nocrop -force
mrconvert ${diffdir}/5TT_fsl_T1space.mif ${diffdir}/5TT_fsl_T1space.nii.gz -force
#rm -f ${diffdir}/5TT_fsl_T1space.mif
$FSLDIR/bin/applywarp --rel --interp=trilinear -i ${diffdir}/5TT_fsl_T1space.nii.gz -r ${regimg} -o ${diffdir}/5TT_fsl_dwi.nii.gz


#--------------------------------Step 3.1: Probalisitic tractography----------------------------------------#
# parameter setting: iFOD2 with ACT (select 5 Million tracts, maximum length as 300, cutoff as 0.05)
# using SIFT2 tract weight for further connectcome computation

#-------------------Probalisitic tractography output dictionary--------------------#
OutputDIR=${MainDIR}/Tractography/iFOD2/tck_files
Subject_output=${OutputDIR}/${SUBJECT}
mkdir -p ${Subject_output}

# Probabilistic tractography with ACT
algo=iFOD2
algo_str="ifod2act5Mfsl"
tckfile=${SUBJECT}_CSD_${algo_str}_cut05_5M.tck
siftfile=${tckfile/.tck/_sift2.txt}

# Probalisitic tractography
tckgen ${diffdir}/RF_wm_dhollander.mif ${Subject_output}/${tckfile} -algorithm ${algo} -seed_dynamic ${diffdir}/RF_wm_dhollander.mif -select 5M -maxlength 300 -cutoff 0.05 -act ${diffdir}/5TT_fsl_dwi.nii.gz -force

#Compute SIFT2 tract weights
tcksift2 ${Subject_output}/${tckfile} ${diffdir}/RF_wm_dhollander.mif ${Subject_output}/${siftfile} -fd_thresh 0.05 -force


#--------------------------------Step 3.2: Deterministic tractography----------------------------------------#
# parameter setting: SD stream with ACT (select 5 Million tracts, maximum length as 300, cutoff as 0.05)
# using SIFT2 tract weight for further connectcome computation

#-------------------Deterministic tractography output dictionary--------------------#
OutputDIR=${MainDIR}/Tractography/SD_Stream/tck_files/
Subject_output=${OutputDIR}/${SUBJECT}
mkdir -p ${Subject_output}

# SD stream tractography with ACT
algo=SD_Stream
algo_str="SDStreamAct5Mfsl"
tckfile=${SUBJECT}_CSD_${algo_str}_cut05_5M.tck
siftfile=${tckfile/.tck/_sift2.txt}

# Deterministic tractography
tckgen ${diffdir}/RF_wm_dhollander.mif ${Subject_output}/${tckfile} -algorithm ${algo} -seed_dynamic ${diffdir}/RF_wm_dhollander.mif -select 5M -maxlength 300 -cutoff 0.05 -act ${diffdir}/5TT_fsl_dwi.nii.gz -force

#Compute SIFT2 tract weights
tcksift2 ${Subject_output}/${tckfile} ${diffdir}/RF_wm_dhollander.mif ${Subject_output}/${siftfile} -fd_thresh 0.05 -force


#--------------------------------Step 4.1: Connectcome generation (Probalisitic tractography)----------------------------------------#

#-------------------Tractography dictionary--------------------#
algo=iFOD2
algo_str="ifod2act5Mfsl"
tckfile=${SUBJECT}_CSD_${algo_str}_cut05_5M.tck
TRACT_DIR=${MainDIR}/Tractography/iFOD2/tck_files
TRACTO_FILE=${TRACT_DIR}/${SUBJECT}/${tckfile}
SIFT_FILE=${TRACT_DIR}/${SUBJECT}/${tckfile/.tck/_sift2.txt}

#--------------------Atlas wrapping: Schaefer atlas-----------------------#
#Warp Schaefer atlas from MNI space to diffusion volume space
atlasfile="Schaefer2018_200Parcels_Kong2022_17Networks_order_FSLMNI152_1mm.nii.gz"
roiname="Schaefer200"
roifile="Schaefer200_dwi.nii.gz"
$FSLDIR/bin/applywarp -i ${PrepareDIR}/${atlasfile} -r ${regimg} -w ${diffdir}/standard2acpc_dc.nii.gz -o ${diffdir}/${roifile} --interp=nn

#--------------------FA calculate-----------------------#
# Estimation of the diffusion tensor
dwi2tensor ${diffdir}/DWI_bc.mif ${diffdir}/DTI_tensor.mif -force

# Compute FA from the tensor
tensor2metric -fa ${diffdir}/FA_map.mif ${diffdir}/DTI_tensor.mif -force
mrconvert ${diffdir}/FA_map.mif ${diffdir}/FA_map.nii.gz -force

#-----------------Output dictionary---------------------#
RESULTSDIR=${MainDIR}/Connectome/iFOD2/${roiname}/${SUBJECT}
mkdir -p ${RESULTSDIR}

#-----------------Output --------------------------#
# Calculate mean FA and length per streamline
tcksample ${TRACTO_FILE} ${diffdir}/FA_map.mif ${RESULTSDIR}/${SUBJECT}_mean_FA_per_streamline.csv -stat_tck mean -force
tckstats ${TRACTO_FILE} -dump ${RESULTSDIR}/${SUBJECT}_length_per_streamline.csv -force

# Calculate the connectivity matrix for mean FA, output FA
tck2connectome -symmetric -zero_diagonal ${TRACTO_FILE} ${diffdir}/${roifile} ${RESULTSDIR}/${SUBJECT}_mean_FA_connectome.csv -scale_file ${RESULTSDIR}/${SUBJECT}_mean_FA_per_streamline.csv -stat_edge mean -tck_weights_in ${SIFT_FILE} -force

# output the counts
tck2connectome -symmetric -zero_diagonal ${TRACTO_FILE} ${diffdir}/${roifile} ${RESULTSDIR}/${SUBJECT}_count_connectome.csv -tck_weights_in ${SIFT_FILE} -out_assignments ${RESULTSDIR}/${SUBJECT}_assignments.txt -force


#--------------------------------Step 4.2: Connectcome generation (Deterministic tractography)----------------------------------------#

#-------------------Tractography dictionary--------------------#
algo=SD_Stream
algo_str="SDStreamAct5Mfsl"
tckfile=${SUBJECT}_CSD_${algo_str}_cut05_5M.tck
TRACT_DIR=${MainDIR}/Tractography/SD_Stream/tck_files
TRACTO_FILE=${TRACT_DIR}/${SUBJECT}/${tckfile}
SIFT_FILE=${TRACT_DIR}/${SUBJECT}/${tckfile/.tck/_sift2.txt}

#--------------------FA calculate-----------------------#
# Don't need to calculate FA again. 

#-----------------Output dictionary--------------------------#

RESULTSDIR=${MainDIR}/Connectome/SD_Stream/${roiname}/${SUBJECT}
mkdir -p ${RESULTSDIR}

#-----------------Output --------------------------#

# Calculate mean FA and length per streamline
tcksample ${TRACTO_FILE} ${diffdir}/FA_map.mif ${RESULTSDIR}/${SUBJECT}_mean_FA_per_streamline.csv -stat_tck mean -force
tckstats ${TRACTO_FILE} -dump ${RESULTSDIR}/${SUBJECT}_length_per_streamline.csv -force

# Calculate the connectivity matrix for mean FA, output FA
tck2connectome -symmetric -zero_diagonal ${TRACTO_FILE} ${diffdir}/${roifile} ${RESULTSDIR}/${SUBJECT}_mean_FA_connectome.csv -scale_file ${RESULTSDIR}/${SUBJECT}_mean_FA_per_streamline.csv -stat_edge mean -tck_weights_in ${SIFT_FILE} -force

# output the counts
tck2connectome -symmetric -zero_diagonal ${TRACTO_FILE} ${diffdir}/${roifile} ${RESULTSDIR}/${SUBJECT}_count_connectome.csv -tck_weights_in ${SIFT_FILE} -out_assignments ${RESULTSDIR}/${SUBJECT}_assignments.txt -force



