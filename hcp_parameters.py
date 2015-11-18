# -*- coding: utf-8 -*-

"""
Definition of the values for each parameter of the hcp_pipeline.

"""

## Location of the folder for the data of a single subject
dir_hcp_subj = '/Users/paolo/Datasets/HCP/RAW_DATA/100307/'
dir_hcp_subj = '/Users/paolo/Datasets/HCP/TEST/124422/'
# 100307 124422 161731 199655 201111 239944 245333 366446 528446 856766

## Location of the subfolder for structural MRI data
dir_hcp_t1 = 'T1w'

## Location of subfolder for diffusion MRI data
dir_hcp_dwi = 'T1w/Diffusion'

## Location of folder for results of hcp_pipeline
dir_hcp_out = 'Tractography'

## b value for shell extraction from multishell diffusion data
par_b_shell = 1000

## file name tag for denotation of single shell extraction
par_b_tag = 'b1k'

## b value threshold for single shell ectraction
par_b0_threshold = 10

## Resampling resolution of target diffusion data
par_dim_vox = 2.0

## file name tag for resampling resolution
par_dim_tag = '2mm'

##
par_fa_tag = 'FA'

##
par_md_tag = 'MD'

##
par_wmfa_tag = 'fa'

##
par_wmp_tag = 'parc'

## Auto response radius
par_ar_radius = 10

## Auto response FA threshold
par_ar_fa_th = 0.7

## Sh order for constrained spherical deconvolution
par_csd_sh_order = 8

## Relative peak threshold for constrained spherical deconvolution
par_csd_peak_th = 0.5

## Minimum angle for constrained spherical deconvolution
par_csd_min_angle = 25.0

## Number of peaks for constrained spherical deconvolution
par_csd_num_peak = 5

## Parallel computation of constrained spherical deconvolution
par_csd_parallel = False

##
par_eudx_seeds = 100000

##
par_eudx_threshold = 0.2

##
par_eudx_tag = 'eudx100k'

##
par_rec_tag = 'dti'

##
par_trk_density = [2, 2, 2]

##
par_trk_fa_th = 0.1

##
par_trk_max_angle = 30.0

##
par_trk_step_size = 0.5

##
par_trk_odf_tag = 'odf'

##
par_trk_prob_tag = 'prob'

##
par_trk_fa_tag = 'fa'

##
par_trk_tag = 'eudx'

##
par_wmql_dir = 'wmql_tracts'

##
par_wmql_opt = ''
