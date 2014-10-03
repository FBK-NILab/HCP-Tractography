# -*- coding: utf-8 -*-
"""
Copyright (c) 2014, Diana Porro and Paolo Avesani
Distributed under the BSD 3-clause license. See COPYING.txt.
"""


import numpy as np
import nibabel as nib
import dipy.reconst.dti as dti
from dipy.core.gradients import gradient_table
from dipy.io.dpy import Dpy
from dipy.data import get_sphere
from dipy.tracking.eudx import EuDX


def tractography_rec(imag, bvals, bvecs, seed, threshold):
    ''' Script to generate tractography. Uses the EuDX function from dipy. Returns tractography and FA.
    
    Parameters
    ----------
    imag: NiftiImage object 
    bvals: bvals array
    bvecs: bvecs array
    seed: int or ndarray (Parameter for the EuDX function)
    threshold : float (Parameter for the EuDX function)
    '''

    print "Retrieving data and affine"
    data = img.get_data()
    affine = img.get_affine()

    #new version of dipy
    print "Computing tensor model"
    gradients = gradient_table(bvals,bvecs)
    tensor_model = dti.TensorModel(gradients)
    tensors = tensor_model.fit(data)

    print "Computing FA"
    FA = dti.fractional_anisotropy(tensors.evals)
    FA[np.isnan(FA)] = 0

    print "Computing evecs"
    evecs_img = nib.Nifti1Image(tensors.evecs.astype(np.float32), affine)
    evecs = evecs_img.get_data()

    sphere = get_sphere('symmetric724')
    peak_indices = dti.quantize_evecs(evecs, sphere.vertices)
    
    print "Computing EuDX reconstruction."
    streamlines = EuDX(FA.astype('f8'),
                        ind=peak_indices,
			seeds=seed,
                        odf_vertices= sphere.vertices,
                        a_low=threshold)

    return streamlines, FA


def save_trk(streamlines, voxel_size, dimensions, filename):
    '''Save tractography to a .trk file'''
    
    print "Save tracks as .trk"
    hdr = nib.trackvis.empty_header()
    hdr['voxel_size'] = voxel_size
    hdr['voxel_order'] = 'LAS'
    hdr['dim'] = dimensions
    strm = ((sl, None, None) for sl in streamlines)

    nib.trackvis.write(filename, strm, hdr, points_space='voxel')

def save_dpy(streamlines, filename):
    ''' Save tractography to a .dpy file'''

    print "Save tracks as .dpy"	
    tracks = [track for track in streamlines]
    dpw = Dpy(filename, 'w')
    dpw.write_tracks(tracks)
    dpw.close()



if __name__ == '__main__':

    dirname = 'HCP/124422/T1w/Diffusion/'
    nii_filename = dirname + 'data.nii.gz'
    bval_filename = dirname + 'bvals'
    bvec_filename = dirname + 'bvecs'

    print "Loading data"
    img = nib.load(nii_filename)
    bvals = np.loadtxt(bval_filename)
    bvecs = np.loadtxt(bvec_filename).T

    print "Setting parameters"
    seed = 100000
    threshold = .1

    print "Calling method to reconstruct tractography"
    streamlines, FA = tractography_rec(img, bvals, bvecs, seed, threshold)

    print "Saving data"
    voxel_size = img.get_header().get_zooms()[:3]
    dims = FA.shape[:3]

    save_filename_trk = dirname + 'tracks_dti_%s.trk' % seed
    save_trk(streamlines, voxel_size, dims, save_filename_trk) 

    save_filename_dpy = dirname + 'tracks_dti_%s.dpy' % seed
    save_dpy(streamlines, save_filename_dpy)

    print "Save FA"
    mapfile = dirname+'DTI/FA_map.nii.gz'
    nib.save(nib.Nifti1Image(FA, img.get_affine()), mapfile)
