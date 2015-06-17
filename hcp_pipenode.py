import os
import sys
import inspect
import numpy as np
import nibabel as nib
from os.path import join as pjoin
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from dipy.align.aniso2iso import resample
from dipy.align.aniso2iso import reslice
from dipy.reconst.dti import TensorModel
from dipy.data import get_sphere
from dipy.tracking.utils import seeds_from_mask
from dipy.reconst.dti import quantize_evecs
from dipy.reconst.csdeconv import ConstrainedSphericalDeconvModel, auto_response
from dipy.direction.probabilistic_direction_getter import ProbabilisticDirectionGetter
from dipy.direction import DeterministicMaximumDirectionGetter
from dipy.tracking.eudx import EuDX
from dipy.tracking.local import ThresholdTissueClassifier, BinaryTissueClassifier
from dipy.tracking import utils
from dipy.tracking.local import LocalTracking
from dipy.tracking.metrics import length
from dipy.io.trackvis import save_trk
from hcp_utils import load_nifti, save_nifti, tract_querier
from hcp_parameters import *
from hcp_wm_mask import *

def single_shell_extraction(dir_src, dir_out, verbose=False):

    fbval = pjoin(dir_src, 'bvals')
    fbvec = pjoin(dir_src, 'bvecs')
    fmask = pjoin(dir_src, 'nodif_brain_mask.nii.gz')
    fdwi = pjoin(dir_src, 'data.nii.gz')

    bvals, bvecs = read_bvals_bvecs(fbval, fbvec)
    data, affine = load_nifti(fdwi, verbose)

    if par_b_shell == 1000:
        sind = (bvals < 10) | ((bvals < 1100) & (bvals > 900))
    elif par_b_shell == 2000:
        sind = (bvals < 10) | ((bvals < 2100) & (bvals > 1900))
    elif par_b_shell == 3000:
        sind = (bvals < 10) | ((bvals < 3100) & (bvals > 2900))

    shell_data = data[..., sind]
    shell_gtab = gradient_table(bvals[sind], bvecs[sind, :],
                                  b0_threshold=par_b0_threshold)

    fname = 'data_' + par_b_tag + '.nii.gz'
    save_nifti(pjoin(dir_out, fname), shell_data, affine)
    np.savetxt(pjoin(dir_out, 'bvals_' + par_b_tag), shell_gtab.bvals)
    np.savetxt(pjoin(dir_out, 'bvecs_' + par_b_tag), shell_gtab.bvecs.T)


def resample_data_resolution(dir_src, dir_out, verbose=False):
    
    fmask = pjoin(dir_src, 'nodif_brain_mask.nii.gz')
    fdwi = pjoin(dir_out, 'data_' + par_b_tag + '.nii.gz')
    
    data, affine = load_nifti(fdwi, verbose)
    mask, _ = load_nifti(fmask, verbose)

    data2, affine2 = reslice(data, affine, (1.25,) * 3, (par_dim_vox,) * 3)
    mask2, _ = reslice(mask, affine, (1.25,) * 3, (par_dim_vox,) * 3, order=0)

    fname = 'data_' + par_b_tag + '_' + par_dim_tag + '.nii.gz'
    save_nifti(pjoin(dir_out, fname), data2, affine2)
    fname = 'nodif_brain_mask_' + par_dim_tag + '.nii.gz'
    save_nifti(pjoin(dir_out, fname), mask2, affine2)
    
    fwmparc = pjoin(dir_src, '../wmparc.nii.gz')
    data, affine = load_nifti(fwmparc, verbose)
    data2, affine2 = reslice(data, affine, (0.7,) * 3, (par_dim_vox,) * 3, order=0)
    fname = 'wmparc_' + par_dim_tag + '.nii.gz'
    save_nifti(pjoin(dir_out, fname), data2, affine2)

    ft1w = pjoin(dir_src, '../T1w_acpc_dc_restore_brain.nii.gz')
    data, affine = load_nifti(ft1w, verbose)
    data2, affine2 = reslice(data, affine, (0.7,) * 3, (par_dim_vox,) * 3, order=0, mode='constant')
    fname = 't1w_acpc_dc_restore_' + par_dim_tag + '.nii.gz'
    save_nifti(pjoin(dir_out, fname), data2, affine2)
    

def compute_tensor_model(dir_src, dir_out, verbose=False):

    fbval = pjoin(dir_src, 'bvals_' + par_b_tag)
    fbvec = pjoin(dir_src, 'bvecs_' + par_b_tag)
    fdwi =  pjoin(dir_src, 'data_' + par_b_tag + '_' + par_dim_tag + '.nii.gz')
    fmask = pjoin(dir_src, 'nodif_brain_mask_' + par_dim_tag + '.nii.gz')

    bvals, bvecs = read_bvals_bvecs(fbval, fbvec)
    gtab = gradient_table(bvals, bvecs, b0_threshold=par_b0_threshold)
    data, affine = load_nifti(fdwi, verbose)
    mask, _ = load_nifti(fmask, verbose)

    ten_model = TensorModel(gtab)
    ten_fit = ten_model.fit(data, mask)

    FA = ten_fit.fa
    MD = ten_fit.md
    EV = ten_fit.evecs.astype(np.float32)

    fa_name = 'data_' + par_b_tag + '_' + par_dim_tag + '_FA.nii.gz'
    save_nifti(pjoin(dir_out, fa_name), FA, affine)
    md_name = 'data_' + par_b_tag + '_' + par_dim_tag + '_MD.nii.gz'
    save_nifti(pjoin(dir_out, md_name), MD, affine)
    ev_name = 'data_' + par_b_tag + '_' + par_dim_tag + '_EV.nii.gz'
    save_nifti(pjoin(dir_out, ev_name), EV, affine)

    
def white_matter_mask_FA(dir_src, dir_out, verbose=False):

    src_fa =  'data_' + par_b_tag + '_' + par_dim_tag + '_FA.nii.gz'
    src_md =  'data_' + par_b_tag + '_' + par_dim_tag + '_MD.nii.gz'
    FA, affine = load_nifti(pjoin(dir_src, src_fa), verbose)
    MD, _ = load_nifti(pjoin(dir_src, src_md), verbose)

    wm_mask = (np.logical_or(FA >= 0.4, 
                             (np.logical_and(FA >= 0.15, MD >= 0.0011/2.))))

    out_wm = 'wm_mask_' + par_b_tag + '_' + par_dim_tag + '.nii.gz'
    save_nifti(pjoin(dir_out, out_wm), wm_mask.astype('f4'), affine)


def white_matter_mask_wmparc(dir_src, dir_out, verbose=False):

    fwmparc_src = pjoin(dir_src, 'wmparc.nii.gz')
    fribbon_src = pjoin(dir_src, 'ribbon.nii.gz')
    wmparc, affine = load_nifti(fwmparc_src)
    ribbon, affine = load_nifti(fribbon_src)
    mask = np.zeros_like(wmparc)

    for label in ribbon_structures:
        mask[ribbon == label] = 1

    for label in wmparc_structures + wmparc_cc_structures:
        mask[wmparc == label] = 1

    for label in wmparc_del_structures + wmparc_del_structures2:
        mask[wmparc == label] = 0

    mask = mask.astype('f8')
    mask2, affine2 = resample(mask, affine,
                              (0.7,) * 3, (par_dim_vox,) * 3, order=0)

    wm_out = "wm_mask_%s_%s_%s.nii.gz" % (par_b_tag, par_dim_tag, par_wmp_tag)
    save_nifti(pjoin(dir_out, wm_out), mask2, affine2)



def constrained_spherical_deconvolution(dir_src, dir_out, verbose=False):

    # Load data
    fbval = pjoin(dir_src, 'bvals_' + par_b_tag)
    fbvec = pjoin(dir_src, 'bvecs_' + par_b_tag)
    fdwi =  pjoin(dir_src, 'data_' + par_b_tag + '_' + par_dim_tag + '.nii.gz')
    #fmask = pjoin(dir_src, 'nodif_brain_mask_' + par_dim_tag + '.nii.gz')
    fmask = pjoin(dir_src, 'wm_mask_' + par_b_tag + '_' + par_dim_tag + '.nii.gz')

    bvals, bvecs = read_bvals_bvecs(fbval, fbvec)
    gtab = gradient_table(bvals, bvecs, b0_threshold=par_b0_threshold)
    data, affine = load_nifti(fdwi, verbose)
    mask, _ = load_nifti(fmask, verbose)

    sphere = get_sphere('symmetric724')

    response, ratio = auto_response(gtab, data, roi_radius=par_ar_radius, 
                                    fa_thr=par_ar_fa_th)
    # print('Response function', response)

    # Model fitting
    csd_model = ConstrainedSphericalDeconvModel(gtab, response)
    csd_fit = csd_model.fit(data, mask=mask)

    # Saving Spherical Harmonic Coefficient
    out_peaks = 'sh_' + par_b_tag + '_' + par_dim_tag + '.nii.gz'
    save_nifti(pjoin(dir_out, out_peaks), csd_fit.shm_coeff, affine)
    #out_B = 'B_' + par_b_tag + '_' + par_dim_tag + '.txt'
    #np.savetxt(pjoin(dir_out, out_B), peaks.B)



def tracking_eudx(dir_src, dir_out, verbose=False):

    # Loading FA and evecs data
    fa_name = 'data_' + par_b_tag + '_' + par_dim_tag + '_FA.nii.gz'
    FA, affine = load_nifti(pjoin(dir_src, fa_name), verbose)

    evecs_name = 'data_' + par_b_tag + '_' + par_dim_tag + '_EV.nii.gz'
    evecs, _ = load_nifti(pjoin(dir_src, evecs_name), verbose)

    # Computation of streamlines
    sphere = get_sphere('symmetric724') 
    peak_indices = quantize_evecs(evecs, sphere.vertices)
    streamlines = EuDX(FA.astype('f8'),
                       ind=peak_indices, 
                       seeds=par_eudx_seeds,
                       odf_vertices= sphere.vertices,
                       a_low=par_eudx_threshold)

    # Saving tractography
    voxel_size =  (par_dim_vox,) * 3
    dims = FA.shape[:3]
    hdr = nib.trackvis.empty_header()
    hdr['voxel_size'] = voxel_size
    hdr['voxel_order'] = 'LAS'
    hdr['dim'] = dims
    hdr['vox_to_ras'] = affine
    strm = ((sl, None, None) for sl in streamlines)
    trk_name = 'tractogram_' + par_b_tag + '_' + par_dim_tag + '_' + par_rec_tag + '_' + par_eudx_tag + '.trk'
    trk_out = os.path.join(dir_out, trk_name)
    nib.trackvis.write(trk_out, strm, hdr, points_space='voxel')    


def tracking_maxodf(dir_src, dir_out, verbose=False):

    wm_name = 'wm_mask_' + par_b_tag + '_' + par_dim_tag + '.nii.gz'
    wm_mask, affine = load_nifti(pjoin(dir_src, wm_name), verbose)

    sh_name = 'sh_' + par_b_tag + '_' + par_dim_tag + '.nii.gz'
    sh, _ = load_nifti(pjoin(dir_src, sh_name), verbose)

    sphere = get_sphere('symmetric724') 

    classifier = ThresholdTissueClassifier(wm_mask.astype('f8'), .5)
    classifier = BinaryTissueClassifier(wm_mask)
    max_dg = DeterministicMaximumDirectionGetter.from_shcoeff(sh, max_angle=par_trk_max_angle, sphere=sphere)
    seeds = utils.seeds_from_mask(wm_mask, density=2, affine=affine)
    streamlines = LocalTracking(max_dg, classifier, seeds, affine, step_size=par_trk_step_size)
    streamlines = list(streamlines)

    trk_name = 'tractogram_' + par_b_tag + '_' + par_dim_tag + '_' + par_trk_odf_tag + '.trk'
    save_trk(pjoin(dir_out, trk_name), streamlines, affine, wm_mask.shape)


def compute_tract_query(dir_src, dir_out, subj, verbose=False):

    trk_name = "tractogram_%s_%s_%s_%s.trk" % (par_b_tag, par_dim_tag, 
                                                par_rec_tag, par_trk_tag)
    trk_file = os.path.join(dir_out, trk_name)
    wmparc_name = 'wmparc_' + par_dim_tag + '.nii.gz'
    wmparc_file = os.path.join(dir_out, wmparc_name)
    # Query file in the directory of module 
    wmql_root = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe()))[0])) 
    wmql_file = os.path.join(wmql_root, 'wmql_paper.qry')
    wmql_out = os.path.join(dir_out, par_wmql_dir)
    if not os.path.exists(wmql_out):
        os.makedirs(wmql_out)
    wmql_prefix = pjoin(wmql_out, subj)
    print wmql_prefix
    tract_querier(trk_file, wmparc_file, wmql_file, wmql_prefix, par_wmql_opt)
