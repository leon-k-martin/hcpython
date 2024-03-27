################################
# MRTrix3                      #
################################


rule create_mif:
    """Create .mif for further processing with MRTrix."""
    input:
        data=rules.threshold_data.output.data,
        bvals=rules.rotate_bvecs_to_struct.output.bvals_struct,
        bvecs=rules.rotate_bvecs_to_struct.output.bvec_struct,
    output:
        dwi=join(MRTrix_out, "DWI.mif"),
        log=join(logdir, "create_mif.done"),
    log:
        stdout=MRTrix_log,
        sterr=MRTrix_log,
    group:
        "preprocessing"
    threads: 1
    resources:
        mem_mb=8 * 1024,
        runtime=10,  # in minutes
    shell:
        """
        mrconvert {input.data} {output.dwi} --fslgrad {input.bvecs} {input.bvals} -datatype float32 -stride 0,0,0,1 -nthreads {threads} -force -info
        echo $(date) > {output.log}
        """


rule tissue_segmentation:
    """
    Freesurfer-based tissue segmentation algorithm
    ==============================================
    Freesurfer's aparc+aseg.mgz is used as input for 5ttgen freesurfer. Works way better than FSL's 5ttgen fsl.
    """
    input:
        t1_brain=join(dir_T1w, "T1w_acpc_dc_restore_brain.nii.gz"),
        t2_brain=join(dir_T1w, "T2w_acpc_dc_restore_brain.nii.gz"),
        aseg=join(dir_T1w, "aparc+aseg.nii.gz"),
    output:
        seg=join(MRTrix_out, "5TT.mif"),
        log=join(logdir, "tissue_segmentation.done"),
    group:
        "preprocessing"
    threads: 16
    resources:
        mem_mb=8 * 1024,
        runtime=10,  # in minutes
    shell:
        """
        # 5ttgen fsl {input.t1_brain} {output.seg} -t2 {input.t2_brain} -premasked -nthreads {threads}

        5ttgen freesurfer {input.aseg} {output.seg} -lut $(python -c 'from hcpy import data; print(data.fs_lut)')

        echo $(date) > {output.log}
        """


rule bias_correction:
    """ANTS-based bias correction."""
    input:
        dwi=rules.create_mif.output.dwi,
    output:
        dwi_bias=join(MRTrix_out, "DWI_bias_ants.mif"),
        dwi_bias_field=join(MRTrix_out, "DWI_bias_ants_field.mif"),
        log=join(logdir, "bias_correction.done"),
    group:
        "preprocessing"
    threads: 1
    resources:
        mem_mb=8 * 1024,
        runtime=10,  # in minutes
    container:
        config["containers"]["mrtrix"]
    shell:
        r"""
        dwibiascorrect ants {input.dwi} {output.dwi_bias} \
            -bias {output.dwi_bias_field} \
            -nthreads {threads} \
            -force -info

        echo $(date) > {output.log}
        """


rule create_mean_b0:
    """Mean image"""
    input:
        dwi_bias=rules.bias_correction.output.dwi_bias,
    output:
        mean_b0=join(MRTrix_out, "meanb0.mif"),
        log=join(logdir, "create_mean_b0.done"),
    group:
        "preprocessing"
    threads: 1
    resources:
        mem_mb=8 * 1024,
        runtime=5,  # in minutes
    shell:
        """
        dwiextract {input.dwi_bias} - -bzero | mrmath - mean {output.mean_b0} -axis 3 -force -info
        echo $(date) > {output.log}
        """


rule extract_response_functions:
    """Extract response function. Uses -stride 0,0,0,1"""
    input:
        dwi_bias=rules.bias_correction.output.dwi_bias,
    output:
        response_wm=join(MRTrix_out, "response_wm.txt"),
        response_gm=join(MRTrix_out, "response_gm.txt"),
        response_csf=join(MRTrix_out, "response_csf.txt"),
        voxels=join(MRTrix_out, "RF_voxels.mif"),
        log=join(logdir, "extract_response_functions.done"),
    group:
        "preprocessing"
    threads: 1
    # container:
    #     config["containers"]["mrtrix"]
    resources:
        mem_mb=8 * 1024,
        runtime=10,  # in minutes
    shell:
        """
        dwi2response dhollander {input.dwi_bias} {output.response_wm} {output.response_gm} {output.response_csf} -voxels {output.voxels} -nthreads {threads} -force -info

        echo $(date) > {output.log}
        """


rule generate_mask:
    """Create DWI brain mask,"""
    input:
        dwi_bias=rules.bias_correction.output.dwi_bias,
    output:
        dwi_mask=join(MRTrix_out, "DWI_mask.mif"),
        log=join(logdir, "generate_mask.done"),
    group:
        "preprocessing"
    threads: 1
    # container:
    #     config["containers"]["mrtrix"]
    resources:
        mem_mb=8 * 1024,
        runtime=10,  # in minutes
    shell:
        """
        dwi2mask {input.dwi_bias} {output.dwi_mask} -nthreads {threads} -force -info
        echo $(date) > {output.log}
        """


rule generate_fod:
    """Estimate fibre orientation distributions from diffusion data using spherical deconvolution"""
    input:
        dwi_bias=rules.bias_correction.output.dwi_bias,
        dwi_mask=rules.generate_mask.output.dwi_mask,
        response_wm=rules.extract_response_functions.output.response_wm,
        response_gm=rules.extract_response_functions.output.response_gm,
        response_csf=rules.extract_response_functions.output.response_csf,
    output:
        fod_wm=join(MRTrix_out, "wmfod.mif"),
        fod_gm=join(MRTrix_out, "gmfod.mif"),
        fod_csf=join(MRTrix_out, "csffod.mif"),
        vf=join(MRTrix_out, "volume_fraction.mif"),
        log=join(logdir, "generate_fod.done"),
    group:
        "preprocessing"
    threads: 1
    # container:
    #     config["containers"]["mrtrix"]
    resources:
        mem_mb=8 * 1024,
        runtime=60,  # in minutes
    shell:
        r"""
        dwi2fod msmt_csd {input.dwi_bias} \
            {input.response_wm} {output.fod_wm} \
            {input.response_gm} {output.fod_gm} \
            {input.response_csf} {output.fod_csf} \
            -mask {input.dwi_mask} \
            -nthreads {threads} -force -info

        mrconvert -coord 3 0 {output.fod_wm} - | mrcat {output.fod_csf} {output.fod_gm} - {output.vf}
        echo $(date) > {output.log}
        """


rule normalize_fod:
    """
    Intensity Normalization
    =======================
    Correct  for  global  intensity. Important when performing group studies!
    """
    input:
        dwi_mask=rules.generate_mask.output.dwi_mask,
        fod_wm=rules.generate_fod.output.fod_wm,
        fod_gm=rules.generate_fod.output.fod_gm,
        fod_csf=rules.generate_fod.output.fod_csf,
    output:
        fod_wm_norm=join(MRTrix_out, "wmfod_norm.mif"),
        fod_gm_norm=join(MRTrix_out, "gmfod_norm.mif"),
        fod_csf_norm=join(MRTrix_out, "csffod_norm.mif"),
        check_norm=join(MRTrix_out, "mtnormalise_norm.mif"),
        check_mask=join(MRTrix_out, "mtnormalise_mask.mif"),
        log=join(logdir, "normalize_fod.done"),
    group:
        "preprocessing"
    threads: 1
    # container:
    #     config["containers"]["mrtrix"]
    resources:
        mem_mb=8 * 1024,
        runtime=10,  # in minutes
    shell:
        r"""
        mtnormalise \
            {input.fod_wm} {output.fod_wm_norm} \
            {input.fod_gm} {output.fod_gm_norm} \
            {input.fod_csf} {output.fod_csf_norm} \
            -mask {input.dwi_mask} \
            -check_norm {output.check_norm} \
            -check_mask {output.check_mask} \
            -nthreads {threads} -force -info

        echo $(date) > {output.log}
        """


rule run_tractography:
    """Create Tractogram with MRTrix tckgen"""
    input:
        fod_wm_norm=rules.normalize_fod.output.fod_wm_norm,
        seg=rules.tissue_segmentation.output.seg,
    output:
        tracks=join(MRTrix_out, "100M_tracks.tck"),
        donwsampled_tracks=join(MRTrix_out, "200k_tracks.tck"),
        seeds=join(MRTrix_out, "seeds.txt"),
        done=join(logdir, "run_tractography.done"),
    group:
        "tractography"
    threads: 32
    # container:
    #     config["containers"]["mrtrix"]
    resources:
        mem_mb=32 * 1024,
        runtime=20 * 60,  # in minutes
    shell:
        r"""
        tckgen \
        -force \
        -nthreads {threads} \
        -algorithm iFOD2 \
        -select 1M \
        -cutoff 0.03 \
        -minlength 2.5 \
        -maxlength 250.0 \
        -act {input.seg} \
        -backtrack \
        -crop_at_gmwmi \
        -max_attempts_per_seed 1000 \
        -seed_dynamic {input.fod_wm_norm} \
        -output_seeds {output.seeds} \
        {input.fod_wm_norm} {output.tracks}

        tckedit {output.tracks} -number 200k {output.donwsampled_tracks}
        echo $(date) > {output.done}
        """


rule qc_fod:
    input:
        tracks=rules.run_tractography.output.tracks,
        t1=join(dir_T1w, "T1w_acpc_dc_restore_brain.nii.gz"),
        fod_wm=rules.generate_fod.output.fod_wm,
        fod_wm_norm=rules.normalize_fod.output.fod_wm_norm,
    output:
        fig=join(qcdir, "FOD.png"),
        log=join(logdir, "qc_fod.done"),
    shell:
        """
        mrconvert -force {input.fod_wm} {input.fod_wm}.nii.gz
        mrconvert -force {input.fod_wm_norm} {input.fod_wm_norm}.nii.gz
        python -c '
from hcpy import plotting
plotting.plot_fod_histogram(fod="{input.fod_wm}.nii.gz", fod_norm="{input.fod_wm_norm}.nii.gz", output_fig="{output.fig}")'
        echo $(date) > {output.log}
        """


rule run_sift:
    """Select tracts using spherical deconvolution"""
    input:
        tracks=rules.run_tractography.output.tracks,
        fod_wm_norm=rules.normalize_fod.output.fod_wm_norm,
        seg=rules.tissue_segmentation.output.seg,
    output:
        sift_weights=join(MRTrix_out, "sift_weights.txt"),
        done=join(logdir, "run_sift.done"),
    group:
        "preprocessing"
    threads: 32
    # container:
    #     config["containers"]["mrtrix"]
    resources:
        mem_mb=8 * 1024,
        runtime=10,  # in minutes
    shell:
        r"""
        tcksift2 \
            {input.tracks} {input.fod_wm_norm} {output.sift_weights} \
            -act {input.seg} \
            -nthreads {threads} -force -info
        echo $(date) > {output.done}
        """


rule qc_tractography:
    input:
        tracks=rules.run_tractography.output.tracks,
        t1=join(dir_T1w, "T1w_acpc_dc_restore_brain.nii.gz"),
        dwi_bias=rules.bias_correction.output.dwi_bias,
        sift_weights=rules.run_sift.output.sift_weights,
    output:
        gif=join(qcdir, "tractography.gif"),
        log=join(logdir, "qc_tractography.done"),
    params:
        fig=join(qcdir, "tractogram"),
    shell:
        """
        mrconvert -force {input.dwi_bias} {input.dwi_bias}.nii.gz

        THRESHOLD=$(python -c 'from hcpy import utils; print(utils.get_sift_threshold("{input.sift_weights}"))')

        tckedit -force {input.tracks} {input.tracks}.sift.tck -tck_weights_in {input.sift_weights} -minweight $THRESHOLD

        python -c '
from hcpy import plotting
plotting.rotating_tractogram("{input.tracks}.sift.tck", "{input.dwi_bias}.nii.gz", "{input.t1}", output_gif="{output.gif}")

plotting.plot_tractogram("{input.tracks}.sift.tck", "{input.dwi_bias}.nii.gz", "{input.t1}", out_path="{params.fig}")'
        """


################################
# Parcellation                 #
################################


rule mmp2surf:
    input:
        mmp_lh=data.mmp_annot_lh,
        mmp_rh=data.mmp_annot_rh,
        fsdir=config["freesurferdir"],
        fshome=config["freesurferhome"],
    output:
        mmp_lh=join(config["freesurferdir"], "{subid}", "label", "lh.hcpmmp1.annot"),
        mmp_rh=join(config["freesurferdir"], "{subid}", "label", "rh.hcpmmp1.annot"),
        done=join(logdir, "mmp2surf.done"),
    # container:
    #     config["containers"]["qunex"]
    shell:
        """
        export SUBJECTS_DIR={input.fsdir}
        ln -sf {input.fshome}/subjects/fsaverage {input.fsdir}/fsaverage

        mri_surf2surf --srcsubject fsaverage --trgsubject {wildcards.subid} --hemi lh --sval-annot {input.mmp_lh} --tval {output.mmp_lh}

        mri_surf2surf --srcsubject fsaverage --trgsubject {wildcards.subid} --hemi rh --sval-annot {input.mmp_rh} --tval {output.mmp_rh}

        echo $(date) > {output.done}
        """


rule mmp2vol:
    input:
        annot_lh=rules.mmp2surf.output.mmp_lh,
        annot_rh=rules.mmp2surf.output.mmp_rh,
        rawavg=join(config["freesurferdir"], "{subid}", "mri", "rawavg.mgz"),
        fsdir=config["freesurferdir"],
    output:
        orig=join(parcdir, "orig_hcpmmp1.mgz"),
        rawavg=join(parcdir, "rawavg_hcpmmp1.nii.gz"),
        done=join(logdir, "mmp2vol.done"),
    # container:
    #     config["containers"]["qunex"]
    shell:
        """
        export SUBJECTS_DIR={input.fsdir}
        mri_aparc2aseg --old-ribbon --s {wildcards.subid} --annot hcpmmp1 --o {output.orig}

        mri_label2vol --seg {output.orig} --temp {input.rawavg} --o {output.rawavg} --regheader {output.orig}

        echo $(date) > {output.done}
        """


rule mmp_vol_to_mrtrix:
    """
    Convert Freesurfer labels to MRTrix format
    """
    input:
        rawavg=rules.mmp2vol.output.rawavg,
        original_order=data.mmp_mrtrix_original,
        mrtrix_order=data.mmp_mrtrix_ordered,
    output:
        labels=join(parcdir, "hcpmmp1.mif"),
        ordered_labels=join(parcdir, "hcpmmp1_ordered.mif"),
        parc=join(parcdir, "sub-{subid}_parc-hcpmmp1_space-T1acpc.nii.gz"),
        done=join(logdir, "mmp_vol_to_mrtrix.done"),
    # container:
    #     config["containers"]["mrtrix"]
    shell:
        """
        mrconvert -datatype uint32 {input.rawavg} {output.labels}

        labelconvert {output.labels} {input.original_order} {input.mrtrix_order} {output.ordered_labels}

        mrconvert -datatype uint32 {output.ordered_labels} {output.parc}

        echo $(date) > {output.done}
        """


rule desikan_vol_to_mrtrix:
    """
    Convert Freesurfer labels to MRTrix format
    """
    input:
        aparc=join(dir_T1w, "aparc+aseg.nii.gz"),
        LUT=join(config["freesurferhome"], "FreeSurferColorLUT.txt"),
        dk_ordered=data.dk_mrtrix_ordered,
    output:
        parc=join(parcdir, "sub-{subid}_parc-desikan_space-T1acpc.nii.gz"),
        done=join(logdir, "desikan_vol_to_mrtrix.done"),
    # container:
    #     config["containers"]["mrtrix"]
    shell:
        """
        labelconvert {input.aparc} {input.LUT} {input.dk_ordered} {output.parc}
        echo $(date) > {output.done}
        """


rule qc_parc:
    input:
        t1=join(dir_T1w, "T1w_acpc_dc_restore_brain.nii.gz"),
        parc=rules.mmp_vol_to_mrtrix.output.parc,
        parc_dk=rules.desikan_vol_to_mrtrix.output.parc,
    output:
        fig=join(qcdir, "qc_mmp.png"),
        fig_dk=join(qcdir, "qc_dk.png"),
        done=join(logdir, "qc_parc.done"),
    shell:
        """
        python -c '
from hcpy import plotting
plotting.plot_mmp2subject("{wildcards.subid}", "{input.t1}", "{input.parc}", "{output.fig}")
plotting.plot_dk2subject("{wildcards.subid}", "{input.t1}", "{input.parc_dk}", "{output.fig_dk}")
'
        echo $(date) > {output.done}
        """


rule create_sc:
    input:
        parc=rules.mmp_vol_to_mrtrix.output.parc,
        tracks=rules.run_tractography.output.tracks,
        sift_weights=rules.run_sift.output.sift_weights,
    output:
        weights=join(MRTrix_out, "sub-{subid}_parc-mmp1_sc_weights.csv"),
        lengths=join(MRTrix_out, "sub-{subid}_parc-mmp1_sc_lengths.csv"),
        assignments=join(MRTrix_out, "assignments.txt"),
        done=join(logdir, "create_sc.done"),
    threads: config["nthreads"]
    # container:
    #     config["containers"]["mrtrix"]
    shell:
        """
        tck2connectome {input.tracks} {input.parc} {output.weights} -tck_weights_in {input.sift_weights} -out_assignments {output.assignments} -symmetric -zero_diagonal -scale_invnodevol -nthreads {threads}

        tck2connectome {input.tracks} {input.parc} {output.lengths} -tck_weights_in {input.sift_weights} -symmetric -zero_diagonal -scale_length -stat_edge mean -nthreads {threads}
        echo $(date) > {output.done}
        """


rule create_sc_dk:
    input:
        parc=rules.desikan_vol_to_mrtrix.output.parc,
        tracks=rules.run_tractography.output.tracks,
        sift_weights=rules.run_sift.output.sift_weights,
    output:
        weights=join(MRTrix_out, "sub-{subid}_parc-DeskianKilliany_sc_weights.csv"),
        lengths=join(MRTrix_out, "sub-{subid}_parc-DeskianKilliany_sc_lengths.csv"),
        done=join(logdir, "create_sc_dk.done"),
    threads: config["nthreads"]
    # container:
    #     config["containers"]["mrtrix"]
    shell:
        """
        tck2connectome {input.tracks} {input.parc} {output.weights} -tck_weights_in {input.sift_weights} -symmetric -zero_diagonal -nthreads {threads}
        tck2connectome {input.tracks} {input.parc} {output.lengths} -tck_weights_in {input.sift_weights} -symmetric -zero_diagonal -scale_length -stat_edge mean -nthreads {threads}
        echo $(date) > {output.done}
        """


# rule roi_tracks:
#     input:
#         assignments=rules.create_sc.output.assignments,
#         tracks=rules.run_tractography.output.tracks,
#         sift_weights=rules.run_sift.output.sift_weights,
#     output:
#         roi_tracks=join(MRTrix_out, "roi_tracks"),
#         done=join(logdir, "roi_tracks.done"),
#     shell:
#         """
#         connectome2tck 100M_tracks.tck assignments.txt node -files per_node -tck_weights_in sift_weights.txt
#         """


rule qc_sc:
    input:
        parc=rules.mmp_vol_to_mrtrix.output.parc,
        weights=rules.create_sc.output.weights,
        lengths=rules.create_sc.output.lengths,
        weights_dk=rules.create_sc_dk.output.weights,
        lengths_dk=rules.create_sc_dk.output.lengths,
        t1=join(dir_T1w, "T1w_acpc_dc_restore_brain.nii.gz"),
        dwi=rules.bias_correction.output.dwi_bias,
    output:
        fig=join(qcdir, "qc_parc-mmp_sc.png"),
        fig_dk=join(qcdir, "qc_parc-dk_sc.png"),
        done=join(logdir, "qc_sc.done"),
    shell:
        """
        python -c '
from hcpy import plotting
plotting.plot_sc("{wildcards.subid}", "{input.weights}", "{input.lengths}", "{output.fig}")
plotting.plot_sc("{wildcards.subid}", "{input.weights_dk}", "{input.lengths_dk}", "{output.fig_dk}")
'
        echo $(date) > {output.done}
        """
