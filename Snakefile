configfile: 'config/config.yml'
rule all:
    input: 
        # Preprocessed diffusion-weighted image 
        expand( 
            "derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_desc-preproc_dwi.nii.gz", 
            subject=config['subjects']
        ),
        # Create data subsets   
        expand(
            "derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_LTE_inds.txt",
            subject=config['subjects']
        ),        
        expand(
            "derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_STE_inds.txt",
            subject=config['subjects']
        ),             
        # Brain mask
        expand(
            "derivatives/preproc/sub-{subject}/sub-{subject}_b0mean.nii.gz",
            subject=config['subjects']
        ),
        expand(
            "derivatives/preproc/sub-{subject}/sub-{subject}_brain_mask_acq-mde_dwi.nii.gz",
            subject=config['subjects']
        ),
        # STE eddy 
        expand( 
            "derivatives/preproc/sub-{subject}/sub-{subject}_STE_eddy.nii.gz", 
            subject=config['subjects']    
        ),                        
        # Eddy-corrected bvec and bval files for LTE
        expand(
            "derivatives/preproc/sub-{subject}/sub-{subject}_diff_preproc_LTE.bvec",
            subject=config['subjects']
        ),     
        expand(
            "derivatives/preproc/sub-{subject}/sub-{subject}_diff_preproc_LTE.bval",
            subject=config['subjects']
        ),
        # DKI_fit 
        expand( 
            "derivatives/DKI-fit/sub-{subject}/dwi/sub-{subject}_tvReg_FA.nii.gz",
            subject=config['subjects']
        ),  
        # Recombined LTE and STE data
        expand(
            "derivatives/preproc/sub-{subject}/sub-{subject}_diff_preproc.nii.gz",
            subject=config['subjects']
        ),     
        # Output for uFA     
        expand(
            "derivatives/uFA-fit/sub-{subject}/sub-{subject}_uFA.nii.gz",
            subject=config['subjects']
        ),     
        # HippUnfold output
        expand(
            "derivatives/hippunfold/sub-{subject}/surf/sub-{subject}_space-T1w_den-0p5mm_label-hipp_surfaces.spec",
            subject=config['subjects']
        ),      
        expand(
            "derivatives/hippunfold/sub-{subject}/surf/sub-{subject}_space-T1w_den-0p5mm_label-dentate_surfaces.spec",
            subject=config['subjects']
        ),    
        expand(
            "derivatives/hippunfold/sub-{subject}/anat/sub-{subject}_hemi-L_space-T1w_desc-subfields_atlas-multihist7_dseg.nii.gz",
            subject=config['subjects']
        ),      
        expand(
            "derivatives/hippunfold/sub-{subject}/anat/sub-{subject}_hemi-R_space-T1w_desc-subfields_atlas-multihist7_dseg.nii.gz",
            subject=config['subjects']
        ),    
        expand(
            "derivatives/hippunfold/sub-{subject}/surf/sub-{subject}_hemi-L_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii",
            subject=config['subjects']
        ),  
        expand(
            "derivatives/hippunfold/sub-{subject}/surf/sub-{subject}_hemi-R_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii",
            subject=config['subjects']
        ),
        expand( 
            "derivatives/hippunfold/sub-{subject}/anat/sub-{subject}_desc-preproc_T1w.nii.gz",
            subject=config['subjects']
        ),
        # Generate T1w Brainmask 
        expand(
            "derivatives/preproc/sub-{subject}/sub-{subject}_brain_mask_acq-mde_T1w.nii.gz", 
            subject=config['subjects']
        ), 
        # b0 to T1w affine 
        expand( 
            "derivatives/registration/sub-{subject}/warps/sub-{subject}_run-from-b0_to-T1w_mode-image_affine.mat", 
            subject=config['subjects']
        ), 
        # Apply affine 
        expand( 
            "derivatives/registration/sub-{subject}/anat/sub-{subject}_space-T1w_b0.nii.gz",
            subject=config['subjects']
        ), 
        # Register uFA 
        expand( 
            "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_uFA.nii.gz",
            subject=config['subjects']
        ), 
        expand(
            "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_uFA_noFWE.nii.gz",
            subject=config['subjects']
        ),
        expand(
            "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_D.nii.gz", 
            subject=config['subjects']
        ),
        expand(
            "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_D_noFWE.nii.gz", 
            subject=config['subjects']
        ),
        expand(
            "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kiso.nii.gz",
            subject=config['subjects']
        ),
        expand( 
            "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kiso_noFWE.nii.gz", 
            subject=config['subjects']
        ),
        expand( 
            "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_FA.nii.gz",
            subject=config['subjects']
        ), 
        expand( 
            "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kpar.nii.gz",
            subject=config['subjects']
        ),
        expand( 
            "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kperp.nii.gz", 
            subject=config['subjects']
        ), 
        expand( 
            "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Wmean.nii.gz",
            subject=config['subjects']
        ), 
        expand( 
            "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Klin.nii.gz",
            subject=config['subjects']
        ),
        expand( 
            "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Klin_noFWE.nii.gz",
            subject=config['subjects']
        ),
        #Get Kaniso 
        expand( 
            "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kaniso.nii.gz",
            subject=config['subjects']
        ), 

        expand( 
            "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kaniso_noFWE.nii.gz",
            subject=config['subjects']
        ),

        #Dataframe
        expand( 
            "derivatives/dataframes/sub-{subject}_hippocampus_onlyDWI.csv",
            subject=config['subjects']
        ), 

        expand( 
            "derivatives/dataframes/dataframe_cat/TLE_HC_hippocampus.csv",
            subject=config['subjects']
        ),

        expand( 
            "derivatives/dataframes/dataframe_cat/TLE_HC_hippocampus_volumes.tsv",
            subject=config['subjects']
        )


# Rule to preprocess diffusion MRI data 
rule preprocess_dMRI:
    input:  
        diff = config['input_dwi'],
        phase = config['input_dwi_phase']
    output: 
        "derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_desc-preproc_dwi.nii.gz"
    threads: workflow.cores
    shell: 
        """
        mrcalc {input.diff} {input.phase} -polar tmp1_{wildcards.subject}.nii -force
        dwidenoise tmp1_{wildcards.subject}.nii img_d_{wildcards.subject}.nii.gz -force
        mrdegibbs img_d_{wildcards.subject}.nii.gz img_de_{wildcards.subject}.nii.gz  -force
        mrcalc -force img_de_{wildcards.subject}.nii.gz -abs {output} -force
        rm tmp1_{wildcards.subject}.nii img_d_{wildcards.subject}.nii.gz img_de_{wildcards.subject}.nii.gz
        """
        
# Rule to generate brain mask 
rule generate_brain_mask: 
    input:
        bval = config['input_bval'],
        bvec = config['input_bvec'],
        preproc = "derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_desc-preproc_dwi.nii.gz"
    output:        
        b0mean="derivatives/preproc/sub-{subject}/sub-{subject}_b0mean.nii.gz",
        mask="derivatives/preproc/sub-{subject}/sub-{subject}_brain_mask_acq-mde_dwi.nii.gz",
    threads: workflow.cores
    shell:
        """
        dwiextract {input.preproc} -fslgrad {input.bvec} {input.bval} - -bzero | mrmath - mean {output.b0mean} -axis 3
        mri_synthstrip --i {output.b0mean} --o {wildcards.subject}_mask.nii.gz
        
        mv {wildcards.subject}_mask.nii.gz {output.mask}
        """

# Rule to extract data subsets
rule create_index_subsets:
    input:
        bmat = config['input_bmat'],
        bval = config['input_bval'],
        bvec = config['input_bvec'],
    output:
        LTE_index="derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_LTE_inds.txt",
        STE_index="derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_STE_inds.txt",
        STE_bvec="derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_STE.bvec",
        STE_bval="derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_STE.bval",
        LTE_bvec="derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_LTE.bvec",
        LTE_bval="derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_LTE.bval"
    threads: workflow.cores/4
    shell:
        """
        cp {input.bmat} {input.bval} {input.bvec} "derivatives/preproc/sub-{wildcards.subject}"
        matlab -batch "addpath('/cifs/baron/pipelines/2024_snsx_recons/matmri/dMRI'); niiExtractLTE([], 'derivatives/preproc/sub-{wildcards.subject}/sub-{wildcards.subject}_acq-mde_dwi.bmat', 'derivatives/preproc/sub-{wildcards.subject}/sub-{wildcards.subject}_acq-mde_dwi.bval', 'derivatives/preproc/sub-{wildcards.subject}/sub-{wildcards.subject}_acq-mde_dwi.bvec', 3);"
        """
        
# Rule to extract subsets 
rule extract_subsets:
    input:
        preproc="derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_desc-preproc_dwi.nii.gz",
        LTE_index="derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_LTE_inds.txt",
        STE_index="derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_STE_inds.txt"
    output:
        LTE="derivatives/preproc/sub-{subject}/sub-{subject}_tmp_LTE.nii",
        STE="derivatives/preproc/sub-{subject}/sub-{subject}_tmp_STE.nii"
    threads: workflow.cores
    shell:
        """
        mrconvert {input.preproc} {output.LTE} -coord 3 $(cat {input.LTE_index})
        mrconvert {input.preproc} {output.STE} -coord 3 $(cat {input.STE_index})
        """
        
# Rule to run eddy for STE
rule eddy_STE:
    input:
        STE="derivatives/preproc/sub-{subject}/sub-{subject}_tmp_STE.nii",
        mask="derivatives/preproc/sub-{subject}/sub-{subject}_brain_mask_acq-mde_dwi.nii.gz",
        bvecs="derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_STE.bvec",
        bvals="derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_STE.bval"
    threads: workflow.cores
    output:
        eddy_out="derivatives/preproc/sub-{subject}/sub-{subject}_STE_eddy.nii.gz",
        bvec_out="derivatives/preproc/sub-{subject}/sub-{subject}_diff_preproc_STE.bvec",
        bval_out="derivatives/preproc/sub-{subject}/sub-{subject}_diff_preproc_STE.bval"
    shell:
        """
        tmp_outputFolder="{resources.tmpdir}/sub-{wildcards.subject}/sub-{wildcards.subject}_eddySTE"
        mkdir -p $tmp_outputFolder
        
        printf "0 1 0 0.05" > {wildcards.subject}_acqparams_STE.txt
        SIZES=$(mrinfo -size {input.STE})
        SIZES2=($SIZES)
        indx=""
        for ((i=1; i<=${{SIZES2[3]}}; i+=1)); do 
            indx="$indx 1"; 
        done
        echo "$indx" > {wildcards.subject}_index_STE.txt
        eddy --imain={input.STE} --mask={input.mask} --acqp={wildcards.subject}_acqparams_STE.txt --index={wildcards.subject}_index_STE.txt --bvecs={input.bvecs} --bvals={input.bvals} --out=$tmp_outputFolder --flm=movement --niter=3 --data_is_shelled --verbose
        cp {resources.tmpdir}/sub-{wildcards.subject}/sub-{wildcards.subject}_eddySTE.eddy_rotated_bvecs {output.bvec_out}
        cp {input.bvals} {output.bval_out}
        mrcalc {resources.tmpdir}/sub-{wildcards.subject}/sub-{wildcards.subject}_eddySTE.nii.gz -abs {output.eddy_out} -force
        rm {wildcards.subject}_acqparams_STE.txt {wildcards.subject}_index_STE.txt
        
        rm -rf $tmp_outputFolder
        rm -r {resources.tmpdir}/sub-{wildcards.subject}
        """
       
# Rule to run eddy for LTE
rule eddy_LTE:
    input:
        LTE="derivatives/preproc/sub-{subject}/sub-{subject}_tmp_LTE.nii",
        mask="derivatives/preproc/sub-{subject}/sub-{subject}_brain_mask_acq-mde_dwi.nii.gz",
        bvecs="derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_LTE.bvec",
        bvals="derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_LTE.bval"
    threads: workflow.cores
    output:
        preproc="derivatives/preproc/sub-{subject}/sub-{subject}_diff_preproc_LTE.nii.gz",
        bvec_out="derivatives/preproc/sub-{subject}/sub-{subject}_diff_preproc_LTE.bvec",
        bval_out="derivatives/preproc/sub-{subject}/sub-{subject}_diff_preproc_LTE.bval"
    shell:
        """
        tmp_outputFolder="{resources.tmpdir}/sub-{wildcards.subject}/sub-{wildcards.subject}_eddyLTE"
        mkdir -p $tmp_outputFolder
        
        printf "0 1 0 0.05" > {wildcards.subject}_acqparams_LTE.txt
        SIZES=$(mrinfo -size {input.LTE})
        # Extracting the size from mrinfo output
        SIZES2=($SIZES)
        indx=""
        for ((i=1; i<=${{SIZES2[3]}}; i+=1)); do 
            indx="$indx 1"
        done
        echo $indx > {wildcards.subject}_index_LTE.txt
        eddy --imain={input.LTE} --mask={input.mask} --acqp={wildcards.subject}_acqparams_LTE.txt --index={wildcards.subject}_index_LTE.txt --bvecs={input.bvecs} --bvals={input.bvals} --out=$tmp_outputFolder --flm=movement --niter=3 --data_is_shelled --verbose
        cp {resources.tmpdir}/sub-{wildcards.subject}/sub-{wildcards.subject}_eddyLTE.eddy_rotated_bvecs {output.bvec_out}
        cp {input.bvals} {output.bval_out}
        mrcalc {resources.tmpdir}/sub-{wildcards.subject}/sub-{wildcards.subject}_eddyLTE.nii.gz -abs {output.preproc} -force
        rm {wildcards.subject}_acqparams_LTE.txt {wildcards.subject}_index_LTE.txt
        
        rm -rf $tmp_outputFolder
        rm -r {resources.tmpdir}/sub-{wildcards.subject}
        """   

# DKI fit parameters 
rule fit_DKI: 
    input: 
        dwi = "derivatives/preproc/sub-{subject}/sub-{subject}_diff_preproc_LTE.nii.gz",
        bvec = "derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_LTE.bvec",
        bval = "derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_LTE.bval", 
        mask = "derivatives/preproc/sub-{subject}/sub-{subject}_brain_mask_acq-mde_dwi.nii.gz"
    params: 
        matmri_path = config['matmri_path'],
        tvReg = config['tvReg']
    output: 
        Dmean = "derivatives/DKI-fit/sub-{subject}/dwi/sub-{subject}_tvReg_Dmean.nii.gz",
        Dpar = "derivatives/DKI-fit/sub-{subject}/dwi/sub-{subject}_tvReg_Dpar.nii.gz",
        Dperp = "derivatives/DKI-fit/sub-{subject}/dwi/sub-{subject}_tvReg_Dperp.nii.gz",
        Dpowder = "derivatives/DKI-fit/sub-{subject}/dwi/sub-{subject}_tvReg_Dpowder.nii.gz",
        FA = "derivatives/DKI-fit/sub-{subject}/dwi/sub-{subject}_tvReg_FA.nii.gz",
        FAvec = "derivatives/DKI-fit/sub-{subject}/dwi/sub-{subject}_tvReg_FAvec.nii.gz",
        Kpar = "derivatives/DKI-fit/sub-{subject}/dwi/sub-{subject}_tvReg_Kpar.nii.gz",
        Kperp = 'derivatives/DKI-fit/sub-{subject}/dwi/sub-{subject}_tvReg_Kperp.nii.gz',
        Wmean = 'derivatives/DKI-fit/sub-{subject}/dwi/sub-{subject}_tvReg_Wmean.nii.gz',
        Wpowder = 'derivatives/DKI-fit/sub-{subject}/dwi/sub-{subject}_tvReg_Wpowder.nii.gz'
    threads: config['threads']
    shell: 
        """
        DKI_tmp_outputFolder="{resources.tmpdir}/sub-{wildcards.subject}/sub-{wildcards.subject}_DKI"
        mkdir -p $DKI_tmp_outputFolder
        outputFile="$DKI_tmp_outputFolder/sub-{wildcards.subject}_tvReg"
        bvec_bval=$(echo "{input.bval}" | cut -f 1 -d '.')

        matlab -batch "addpath(genpath('{params.matmri_path}')); opt.tvReg = {params.tvReg};nii2kurt('{input.dwi}', '$bvec_bval', [], '{input.mask}', opt, '$outputFile');"
        
        mv $DKI_tmp_outputFolder/sub-{wildcards.subject}_tvReg_Dmean.nii.gz {output.Dmean}
        mv $DKI_tmp_outputFolder/sub-{wildcards.subject}_tvReg_Dpar.nii.gz {output.Dpar}
        mv $DKI_tmp_outputFolder/sub-{wildcards.subject}_tvReg_Dperp.nii.gz {output.Dperp}
        mv $DKI_tmp_outputFolder/sub-{wildcards.subject}_tvReg_Dpowder.nii.gz {output.Dpowder}
        mv $DKI_tmp_outputFolder/sub-{wildcards.subject}_tvReg_FA.nii.gz {output.FA}
        mv $DKI_tmp_outputFolder/sub-{wildcards.subject}_tvReg_FAvec.nii.gz {output.FAvec}
        mv $DKI_tmp_outputFolder/sub-{wildcards.subject}_tvReg_Kpar.nii.gz {output.Kpar}
        mv $DKI_tmp_outputFolder/sub-{wildcards.subject}_tvReg_Kperp.nii.gz {output.Kperp}
        mv $DKI_tmp_outputFolder/sub-{wildcards.subject}_tvReg_Wmean.nii.gz {output.Wmean}
        mv $DKI_tmp_outputFolder/sub-{wildcards.subject}_tvReg_Wpowder.nii.gz {output.Wpowder}

        rm -r $DKI_tmp_outputFolder
        """

# Rule to recombine LTE and STE
rule recombine_LTE_STE:
    input:
        LTE="derivatives/preproc/sub-{subject}/sub-{subject}_diff_preproc_LTE.nii.gz",
        STE="derivatives/preproc/sub-{subject}/sub-{subject}_STE_eddy.nii.gz",
        LTE_inds="derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_LTE_inds.txt",
        STE_inds="derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_STE_inds.txt",
        LTE_bvec="derivatives/preproc/sub-{subject}/sub-{subject}_diff_preproc_LTE.bvec",
        STE_bvec="derivatives/preproc/sub-{subject}/sub-{subject}_acq-mde_dwi_b0_STE.bvec" 
    threads: workflow.cores
    output:
        combined="derivatives/preproc/sub-{subject}/sub-{subject}_diff_preproc.nii.gz"
    shell:
        """
        matlab -batch "addpath('/cifs/baron/pipelines/2024_snsx_recons/matmri/dMRI'); niiCombine('{input.LTE}', '{input.STE}', '{input.LTE_inds}', '{input.STE_inds}', 'tmp_eddyout.nii', '{input.LTE_bvec}', '{input.STE_bvec}');"
        mrcalc -force tmp_eddyout.nii -abs {output.combined}
        rm tmp_eddyout.nii
        """

# Rule uFA       
rule nii2uFA_fwe:
    input:
        dwi = "derivatives/preproc/sub-{subject}/sub-{subject}_diff_preproc.nii.gz",
        bval = config['input_bval'],
        bmat = config['input_bmat'],
        maskfile = "derivatives/preproc/sub-{subject}/sub-{subject}_brain_mask_acq-mde_dwi.nii.gz"
    params:
        matmri_path = config['matmri_path']
    threads: workflow.cores
    output:
        uFA = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_uFA.nii.gz",
        Kiso = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_Kiso.nii.gz",
        Klin = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_Klin.nii.gz",
        D = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_D.nii.gz",
        sf = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_sf.nii.gz",
        uFA_noFWE = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_uFA_noFWE.nii.gz",
        Klin_noFWE = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_Klin_noFWE.nii.gz",
        Kiso_noFWE = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_Kiso_noFWE.nii.gz",
        D_noFWE = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_D_noFWE.nii.gz",
        sSTE = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_sSTE.nii.gz",
        sLTE = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_sLTE.nii.gz"
    shell:
        """
        uFA_tmp_outputFolder="{resources.tmpdir}/sub-{wildcards.subject}/sub-{wildcards.subject}_uFA"
        mkdir -p $uFA_tmp_outputFolder
        
        matlab -batch "addpath(genpath('{params.matmri_path}'));[uFA, Kiso, Klin, D, sf, uFA_noFWE, Klin_noFWE, Kiso_noFWE, D_noFWE, sSTE, sLTE] = nii2uFA_fwe('{input.dwi}', '{input.bval}', '{input.bmat}', '{input.maskfile}', [], [], ['$uFA_tmp_outputFolder' filesep 'sub-{wildcards.subject}']);"

        cp $uFA_tmp_outputFolder/sub-{wildcards.subject}_uFA.nii.gz {output.uFA}
        cp $uFA_tmp_outputFolder/sub-{wildcards.subject}_Kiso.nii.gz {output.Kiso}
        cp $uFA_tmp_outputFolder/sub-{wildcards.subject}_Klin.nii.gz {output.Klin}
        cp $uFA_tmp_outputFolder/sub-{wildcards.subject}_D.nii.gz {output.D}
        cp $uFA_tmp_outputFolder/sub-{wildcards.subject}_sigFrac.nii.gz {output.sf}
        cp $uFA_tmp_outputFolder/sub-{wildcards.subject}_uFA_noFWE.nii.gz {output.uFA_noFWE}
        cp $uFA_tmp_outputFolder/sub-{wildcards.subject}_Klin_noFWE.nii.gz {output.Klin_noFWE}
        cp $uFA_tmp_outputFolder/sub-{wildcards.subject}_Kiso_noFWE.nii.gz {output.Kiso_noFWE}
        cp $uFA_tmp_outputFolder/sub-{wildcards.subject}_D_noFWE.nii.gz {output.D_noFWE}
        cp $uFA_tmp_outputFolder/sub-{wildcards.subject}_sSTE.nii.gz {output.sSTE}
        cp $uFA_tmp_outputFolder/sub-{wildcards.subject}_sLTE.nii.gz {output.sLTE}
        
        rm -r $uFA_tmp_outputFolder
        """

# HippUnfold Rule
rule hippunfold:
    input: # T1w in the bids folder
        'Data/sub-{subject}'
    params:
        hippunfold_container= config['hippunfold_container'],
        modality = 'T1w'
    output: 
        anat = 'derivatives/hippunfold/sub-{subject}/anat/sub-{subject}_desc-preproc_T1w.nii.gz',
        hipp_surf_spec = 'derivatives/hippunfold/sub-{subject}/surf/sub-{subject}_space-T1w_den-0p5mm_label-hipp_surfaces.spec',
        dentate_surf_spec = 'derivatives/hippunfold/sub-{subject}/surf/sub-{subject}_space-T1w_den-0p5mm_label-dentate_surfaces.spec',
        L_dseg= 'derivatives/hippunfold/sub-{subject}/anat/sub-{subject}_hemi-L_space-T1w_desc-subfields_atlas-multihist7_dseg.nii.gz',
        L_hipp_midthickness = 'derivatives/hippunfold/sub-{subject}/surf/sub-{subject}_hemi-L_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii',
        L_hippUnf_midthickness = 'derivatives/hippunfold/sub-{subject}/surf/sub-{subject}_hemi-L_space-unfold_den-0p5mm_label-hipp_midthickness.surf.gii',
        L_hipp_outer_surf = 'derivatives/hippunfold/sub-{subject}/surf/sub-{subject}_hemi-L_space-T1w_den-0p5mm_label-hipp_outer.surf.gii',
        L_hipp_inner_surf = 'derivatives/hippunfold/sub-{subject}/surf/sub-{subject}_hemi-L_space-T1w_den-0p5mm_label-hipp_inner.surf.gii',
        R_dseg= 'derivatives/hippunfold/sub-{subject}/anat/sub-{subject}_hemi-R_space-T1w_desc-subfields_atlas-multihist7_dseg.nii.gz',
        R_hipp_midthickness = 'derivatives/hippunfold/sub-{subject}/surf/sub-{subject}_hemi-R_space-T1w_den-0p5mm_label-hipp_midthickness.surf.gii',
        R_hippUnf_midthickness = 'derivatives/hippunfold/sub-{subject}/surf/sub-{subject}_hemi-R_space-unfold_den-0p5mm_label-hipp_midthickness.surf.gii',
        R_hipp_outer_surf = 'derivatives/hippunfold/sub-{subject}/surf/sub-{subject}_hemi-R_space-T1w_den-0p5mm_label-hipp_outer.surf.gii',
        R_hipp_inner_surf = 'derivatives/hippunfold/sub-{subject}/surf/sub-{subject}_hemi-R_space-T1w_den-0p5mm_label-hipp_inner.surf.gii',

    threads: config['threads']
    shell:
        """
        input_tmpdir="/tmp/hippunfold_input_sub-{wildcards.subject}/"
        tmpdir="/tmp/hippunfold_sub-{wildcards.subject}/"
        mkdir -p $input_tmpdir 
        mkdir -p $tmpdir 

        cp  -r {input} $input_tmpdir

        singularity run -e {params.hippunfold_container} $input_tmpdir $tmpdir participant --participant_label {wildcards.subject} --cores {threads} --modality {params.modality}

        cp -r $tmpdir"hippunfold/sub-{wildcards.subject}/anat" derivatives/hippunfold/sub-{wildcards.subject}/
        cp -r $tmpdir"hippunfold/sub-{wildcards.subject}/coords" derivatives/hippunfold/sub-{wildcards.subject}/
        cp -r $tmpdir"hippunfold/sub-{wildcards.subject}/qc" derivatives/hippunfold/sub-{wildcards.subject}/
        cp -r $tmpdir"hippunfold/sub-{wildcards.subject}/surf" derivatives/hippunfold/sub-{wildcards.subject}/

        rm -r $tmpdir
        rm -r $input_tmpdir
        """ 

# Rule to generate T1w brain mask 
rule generate_T1w_brain_mask: 
    input:
        preproc = "derivatives/hippunfold/sub-{subject}/anat/sub-{subject}_desc-preproc_T1w.nii.gz"
    output:        
        mask = "derivatives/preproc/sub-{subject}/sub-{subject}_brain_mask_acq-mde_T1w.nii.gz"
    threads: workflow.cores
    shell:
        """
        mri_synthstrip --i {input.preproc} --o {wildcards.subject}_mask.nii.gz
        
        mv {wildcards.subject}_mask.nii.gz {output.mask}
        """    

# Rule to register b0 mask to T1w mask 
rule reg_b0_mask_to_T1w_mask: 
    input: 
        b0="derivatives/preproc/sub-{subject}/sub-{subject}_brain_mask_acq-mde_dwi.nii.gz",
        T1w="derivatives/preproc/sub-{subject}/sub-{subject}_brain_mask_acq-mde_T1w.nii.gz"
    params: 
        metric = 'NMI',
        iterations = '100x50x1',
        greedy_container= config['greedy_container']
    output:
        affine = "derivatives/registration/sub-{subject}/warps/sub-{subject}_run-from-b0_to-T1w_mode-image_affine.mat"
    threads: workflow.cores
    shell:
        """
        tmpdir=$(mktemp -d /tmp/affine_sub-{wildcards.subject}_XXXXXX)
        cp {input.b0}  $tmpdir/sub-{wildcards.subject}_b0.nii.gz
        cp {input.T1w} $tmpdir/sub-{wildcards.subject}_T1w.nii.gz

        mkdir -p $(dirname {output.affine})

        {params.greedy_container} -d 3 -a \
        -m {params.metric} \
        -i $tmpdir/sub-{wildcards.subject}_T1w.nii.gz $tmpdir/sub-{wildcards.subject}_b0.nii.gz \
        -o $tmpdir/sub-{wildcards.subject}_run-from-b0_to-T1w_mode-image_affine.mat \
        -n {params.iterations}

        cp $tmpdir/sub-{wildcards.subject}_run-from-b0_to-T1w_mode-image_affine.mat {output.affine}

        rm -r $tmpdir
        """

# Rule to apply affine 
rule apply_affine: 
    input: 
        b0="derivatives/preproc/sub-{subject}/sub-{subject}_brain_mask_acq-mde_dwi.nii.gz",
        T1w="derivatives/preproc/sub-{subject}/sub-{subject}_brain_mask_acq-mde_T1w.nii.gz",
        affine="derivatives/registration/sub-{subject}/warps/sub-{subject}_run-from-b0_to-T1w_mode-image_affine.mat"
    params:
        interpolation = 'LABEL 0.2vox',
        greedy_container= config['greedy_container']
    output: 
        application = "derivatives/registration/sub-{subject}/anat/sub-{subject}_space-T1w_b0.nii.gz"
    threads: config['threads']
    shell:
        """
        tmpdir=$(mktemp -d /tmp/affine_sub-{wildcards.subject}_XXXXXX)
        cp {input.b0} $tmpdir/sub-{wildcards.subject}_b0.nii.gz 
        cp {input.T1w} $tmpdir/sub-{wildcards.subject}_T1w.nii.gz 
        cp {input.affine} $tmpdir/sub-{wildcards.subject}_affine.mat

        mkdir -p $(dirname {output.application})

        {params.greedy_container} -d 3 -rf $tmpdir/sub-{wildcards.subject}_T1w.nii.gz \
        -rm $tmpdir/sub-{wildcards.subject}_b0.nii.gz $tmpdir/sub-{wildcards.subject}_space-T1w_b0.nii.gz -ri {params.interpolation} \
        -r $tmpdir/sub-{wildcards.subject}_affine.mat

        cp $tmpdir/sub-{wildcards.subject}_space-T1w_b0.nii.gz {output.application}

        rm -r $tmpdir
        """

# Rule to register uFA
rule reg_uFA:
    input: 
        T1w_ref = "derivatives/registration/sub-{subject}/anat/sub-{subject}_space-T1w_b0.nii.gz",
        uFA = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_uFA.nii.gz",
        uFA_noFWE = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_uFA_noFWE.nii.gz",
        D = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_D.nii.gz", 
        D_noFWE = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_D_noFWE.nii.gz",
        Kiso = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_Kiso.nii.gz",
        Kiso_noFWE = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_Kiso_noFWE.nii.gz",
        Klin = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_Klin.nii.gz",
        Klin_noFWE = "derivatives/uFA-fit/sub-{subject}/sub-{subject}_Klin_noFWE.nii.gz",
        FA = "derivatives/DKI-fit/sub-{subject}/dwi/sub-{subject}_tvReg_FA.nii.gz",
        warp = "derivatives/registration/sub-{subject}/warps/sub-{subject}_run-from-b0_to-T1w_mode-image_affine.mat",
        Kpar = "derivatives/DKI-fit/sub-{subject}/dwi/sub-{subject}_tvReg_Kpar.nii.gz",
        Kperp = "derivatives/DKI-fit/sub-{subject}/dwi/sub-{subject}_tvReg_Kperp.nii.gz",
        Wmean = "derivatives/DKI-fit/sub-{subject}/dwi/sub-{subject}_tvReg_Wmean.nii.gz"
    params:
        interpolation = 'LABEL 0.2vox',
        greedy_container= config['greedy_container']
    output: 
        uFA = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_uFA.nii.gz",
        uFA_noFWE = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_uFA_noFWE.nii.gz",
        D = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_D.nii.gz",
        D_noFWE = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_D_noFWE.nii.gz",
        Kiso = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kiso.nii.gz",
        Kiso_noFWE = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kiso_noFWE.nii.gz",
        Klin = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Klin.nii.gz",
        Klin_noFWE = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Klin_noFWE.nii.gz",
        FA = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_FA.nii.gz",
        Kpar = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kpar.nii.gz",
        Kperp = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kperp.nii.gz", 
        Wmean = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Wmean.nii.gz"

    threads: config['threads']
    shell:
        """
        mkdir -p $(dirname {output.uFA})
        mkdir -p $(dirname {output.uFA_noFWE})
        mkdir -p $(dirname {output.D})
        mkdir -p $(dirname {output.D_noFWE})
        mkdir -p $(dirname {output.Kiso})
        mkdir -p $(dirname {output.Kiso_noFWE})
        mkdir -p $(dirname {output.Klin})
        mkdir -p $(dirname {output.Klin_noFWE})
        mkdir -p $(dirname {output.FA})
        mkdir -p $(dirname {output.Kpar})
        mkdir -p $(dirname {output.Kperp})
        mkdir -p $(dirname {output.Wmean})

        tmpdir=$(mktemp -d /tmp/affine_sub-{wildcards.subject}_XXXXXX)
        cp {input.T1w_ref} $tmpdir/sub-{wildcards.subject}_T1w_ref.nii.gz
        cp {input.uFA} $tmpdir/sub-{wildcards.subject}_uFA.nii.gz
        cp {input.uFA_noFWE} $tmpdir/sub-{wildcards.subject}_uFA_noFWE.nii.gz
        cp {input.D} $tmpdir/sub-{wildcards.subject}_D.nii.gz
        cp {input.D_noFWE} $tmpdir/sub-{wildcards.subject}_D_noFWE.nii.gz
        cp {input.Kiso} $tmpdir/sub-{wildcards.subject}_Kiso.nii.gz
        cp {input.Kiso_noFWE} $tmpdir/sub-{wildcards.subject}_Kiso_noFWE.nii.gz
        cp {input.warp} $tmpdir/sub-{wildcards.subject}_warp.mat
        cp {input.Klin} $tmpdir/sub-{wildcards.subject}_Klin.nii.gz
        cp {input.Klin_noFWE} $tmpdir/sub-{wildcards.subject}_Klin_noFWE.nii.gz
        cp {input.FA} $tmpdir/sub-{wildcards.subject}_FA.nii.gz
        cp {input.Kpar} $tmpdir/sub-{wildcards.subject}_Kpar.nii.gz
        cp {input.Kperp} $tmpdir/sub-{wildcards.subject}_Kperp.nii.gz 
        cp {input.Wmean} $tmpdir/sub-{wildcards.subject}_Wmean.nii.gz

        {params.greedy_container} -d 3 -rf $tmpdir/sub-{wildcards.subject}_T1w_ref.nii.gz -rm $tmpdir/sub-{wildcards.subject}_uFA.nii.gz $tmpdir/sub-{wildcards.subject}_uFA_out.nii.gz -ri {params.interpolation} -r $tmpdir/sub-{wildcards.subject}_warp.mat
        {params.greedy_container} -d 3 -rf $tmpdir/sub-{wildcards.subject}_T1w_ref.nii.gz -rm $tmpdir/sub-{wildcards.subject}_uFA_noFWE.nii.gz $tmpdir/sub-{wildcards.subject}_uFA_noFWE_out.nii.gz -ri {params.interpolation} -r $tmpdir/sub-{wildcards.subject}_warp.mat
        {params.greedy_container} -d 3 -rf $tmpdir/sub-{wildcards.subject}_T1w_ref.nii.gz -rm $tmpdir/sub-{wildcards.subject}_D.nii.gz $tmpdir/sub-{wildcards.subject}_D_out.nii.gz -ri {params.interpolation} -r $tmpdir/sub-{wildcards.subject}_warp.mat
        {params.greedy_container} -d 3 -rf $tmpdir/sub-{wildcards.subject}_T1w_ref.nii.gz -rm $tmpdir/sub-{wildcards.subject}_D_noFWE.nii.gz $tmpdir/sub-{wildcards.subject}_D_noFWE_out.nii.gz -ri {params.interpolation} -r $tmpdir/sub-{wildcards.subject}_warp.mat
        {params.greedy_container} -d 3 -rf $tmpdir/sub-{wildcards.subject}_T1w_ref.nii.gz -rm $tmpdir/sub-{wildcards.subject}_Kiso.nii.gz $tmpdir/sub-{wildcards.subject}_Kiso_out.nii.gz -ri {params.interpolation} -r $tmpdir/sub-{wildcards.subject}_warp.mat
        {params.greedy_container} -d 3 -rf $tmpdir/sub-{wildcards.subject}_T1w_ref.nii.gz -rm $tmpdir/sub-{wildcards.subject}_Kiso_noFWE.nii.gz $tmpdir/sub-{wildcards.subject}_Kiso_noFWE_out.nii.gz -ri {params.interpolation} -r $tmpdir/sub-{wildcards.subject}_warp.mat
        {params.greedy_container} -d 3 -rf $tmpdir/sub-{wildcards.subject}_T1w_ref.nii.gz -rm $tmpdir/sub-{wildcards.subject}_Klin.nii.gz $tmpdir/sub-{wildcards.subject}_Klin_out.nii.gz -ri {params.interpolation} -r $tmpdir/sub-{wildcards.subject}_warp.mat
        {params.greedy_container} -d 3 -rf $tmpdir/sub-{wildcards.subject}_T1w_ref.nii.gz -rm $tmpdir/sub-{wildcards.subject}_Klin_noFWE.nii.gz $tmpdir/sub-{wildcards.subject}_Klin_noFWE_out.nii.gz -ri {params.interpolation} -r $tmpdir/sub-{wildcards.subject}_warp.mat
        {params.greedy_container} -d 3 -rf $tmpdir/sub-{wildcards.subject}_T1w_ref.nii.gz -rm $tmpdir/sub-{wildcards.subject}_FA.nii.gz $tmpdir/sub-{wildcards.subject}_FA_out.nii.gz -ri {params.interpolation} -r $tmpdir/sub-{wildcards.subject}_warp.mat
        {params.greedy_container} -d 3 -rf $tmpdir/sub-{wildcards.subject}_T1w_ref.nii.gz -rm $tmpdir/sub-{wildcards.subject}_Kpar.nii.gz $tmpdir/sub-{wildcards.subject}_Kpar_out.nii.gz -ri {params.interpolation} -r $tmpdir/sub-{wildcards.subject}_warp.mat
        {params.greedy_container} -d 3 -rf $tmpdir/sub-{wildcards.subject}_T1w_ref.nii.gz -rm $tmpdir/sub-{wildcards.subject}_Kperp.nii.gz $tmpdir/sub-{wildcards.subject}_Kperp_out.nii.gz -ri {params.interpolation} -r $tmpdir/sub-{wildcards.subject}_warp.mat
        {params.greedy_container} -d 3 -rf $tmpdir/sub-{wildcards.subject}_T1w_ref.nii.gz -rm $tmpdir/sub-{wildcards.subject}_Wmean.nii.gz $tmpdir/sub-{wildcards.subject}_Wmean_out.nii.gz -ri {params.interpolation} -r $tmpdir/sub-{wildcards.subject}_warp.mat
        
        cp $tmpdir/sub-{wildcards.subject}_uFA_out.nii.gz {output.uFA}
        cp $tmpdir/sub-{wildcards.subject}_uFA_noFWE_out.nii.gz {output.uFA_noFWE}
        cp $tmpdir/sub-{wildcards.subject}_D_out.nii.gz {output.D}
        cp $tmpdir/sub-{wildcards.subject}_D_noFWE_out.nii.gz {output.D_noFWE}
        cp $tmpdir/sub-{wildcards.subject}_Kiso_out.nii.gz {output.Kiso}
        cp $tmpdir/sub-{wildcards.subject}_Kiso_noFWE_out.nii.gz {output.Kiso_noFWE}
        cp $tmpdir/sub-{wildcards.subject}_FA_out.nii.gz {output.FA}
        cp $tmpdir/sub-{wildcards.subject}_Klin_out.nii.gz {output.Klin}
        cp $tmpdir/sub-{wildcards.subject}_Klin_noFWE_out.nii.gz {output.Klin_noFWE}
        cp $tmpdir/sub-{wildcards.subject}_Kpar_out.nii.gz {output.Kpar}
        cp $tmpdir/sub-{wildcards.subject}_Kperp_out.nii.gz {output.Kperp}
        cp $tmpdir/sub-{wildcards.subject}_Wmean_out.nii.gz {output.Wmean}

        rm -r $tmpdir

        """
#Kaniso rule        
rule get_Kaniso:
    input:
        Kiso = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kiso.nii.gz",
        Kiso_noFWE = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kiso_noFWE.nii.gz",
        Klte = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Klin.nii.gz",
        Klte_noFWE = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Klin_noFWE.nii.gz",
    output:
        Kaniso = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kaniso.nii.gz",
        Kaniso_noFWE = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kaniso_noFWE.nii.gz"
    threads: config['threads']
    shell: 
        """
        mrcalc {input.Klte} {input.Kiso} -subtract {output.Kaniso}
        mrcalc {input.Klte_noFWE} {input.Kiso_noFWE} -subtract {output.Kaniso_noFWE}
        """

#DATA Frame
rule subject_onlyDWI_dataframe:
    input:
        # uFA inputs
        uFA = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_uFA.nii.gz",
        uFA_noFWE = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_uFA_noFWE.nii.gz",
        D = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_D.nii.gz",
        D_noFWE = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_D_noFWE.nii.gz",
        Kaniso = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kaniso.nii.gz",
        Kaniso_noFWE = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kaniso_noFWE.nii.gz",
        Kiso = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kiso.nii.gz",
        Kiso_noFWE = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kiso_noFWE.nii.gz",
        # DKI inputs
        Klin = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Klin.nii.gz",
        Klin_noFWE = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Klin_noFWE.nii.gz",
        FA = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_FA.nii.gz",
        Kpar = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kpar.nii.gz",
        Kperp = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Kperp.nii.gz", 
        Wmean = "derivatives/registration/sub-{subject}/dwi/sub-{subject}_space-T1w_Wmean.nii.gz",
        # Hippunfold segmentation inputs
        L_dseg = 'derivatives/hippunfold/sub-{subject}/anat/sub-{subject}_hemi-L_space-T1w_desc-subfields_atlas-multihist7_dseg.nii.gz',
        R_dseg = 'derivatives/hippunfold/sub-{subject}/anat/sub-{subject}_hemi-R_space-T1w_desc-subfields_atlas-multihist7_dseg.nii.gz',
        # participants list
        participants = 'data/participants.tsv'
    output:
        "derivatives/dataframes/sub-{subject}_hippocampus_onlyDWI.csv"
    conda:
        "dataframes"
    threads: config['threads']/2
    script:
        "scripts/subject_onlyDWI_dataframe.py"

rule concatenate_dataframes:
    input:
        expand(
            "derivatives/dataframes/sub-{subject}_hippocampus_onlyDWI.csv",          
            subject=config['subjects'],
            run=config['runs'],
            tvReg=config['tvReg']
        )
    output:
        "derivatives/dataframes/dataframe_cat/TLE_HC_hippocampus.csv"
    shell:
        "awk 'NR==1 || FNR!=1' {input}   > {output}"

rule concatenate_volumes_dataframes:
    input:
        expand(
            "derivatives/hippunfold/sub-{subject}/anat/sub-{subject}_space-cropT1w_desc-subfields_atlas-multihist7_volumes.tsv",          
            subject=config['subjects'] 
        )
    output:
        "derivatives/dataframes/dataframe_cat/TLE_HC_hippocampus_volumes.tsv"
    shell:
        "awk 'NR==1 || FNR!=1' {input} > {output}"
