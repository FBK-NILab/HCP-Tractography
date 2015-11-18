import nibabel as nib
import numpy as np
from subprocess import Popen, PIPE
from dipy.tracking.distances import bundles_distances_mam


def load_nifti(fname, verbose=False):
    img = nib.load(fname)
    data = img.get_data()
    affine = img.get_affine()
    if verbose:
        print('Loading...')
        print(fname)
        print(data.shape)
        print(affine)
        print(img.get_header().get_zooms()[:3])
        print(nib.aff2axcodes(affine))
        print
    return data, affine


def save_nifti(fname, data, affine, verbose=False):
    if verbose:
        print('Saving...')
        print(fname)
    nib.save(nib.Nifti1Image(data, affine), fname)


def pipe(cmd, print_sto=True, print_ste=True):
    """Open a pipe to a subprocess where execute an external command.
    """
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    sto = p.stdout.readlines()
    ste = p.stderr.readlines()
    if print_sto :
        print(sto)
    if print_ste :
        print(ste)

def tract_querier(ftrk, fatlas, fquery, fout, options='-f'):
    cmd = "tract_querier -t %s -a %s -q %s -o %s %s" % (ftrk, fatlas, fquery, fout, options)
    print cmd
    pipe(cmd, print_sto=False, print_ste=False)

def tract_converter(fsrc, fout):
    cmd = "TractConverter.py -i %s -o %s" % (fsrc, fout)
    pipe(cmd, print_sto=False, print_ste=False)


def save_peaks(peaks_dir, peaks, tag=''):

    save_nifti(pjoin(peaks_dir, 'peak_dirs' + tag + '.nii.gz'),
               peaks.peak_dirs)

    save_nifti(pjoin(peaks_dir, 'peak_values' + tag + '.nii.gz'),
               peaks.peak_values)

    save_nifti(pjoin(peaks_dir, 'peak_indices' + tag + '.nii.gz'),
               peaks.peak_indices)


def load_peaks(dname, peaks):

    pass


def shannon_entropy(tract):
    '''
    compute the Shannon Entropy of a set of tracks as defined by Lauren
    H(A) = (-1/|A|) * sum{log(1/|A|)* sum[p(f_i|f_j)]} 
    where p(f_i|f_j) = exp(-d(f_i,f_j)*d(f_i,f_j))
    '''

    dm = bundles_distances_mam(tract, tract)
    dm = np.array(dm, dtype =float)
    A = len(tract)
    theta = 10.
    pr_all = np.exp((-dm**2)/theta)    
    pr_i = (1./A) * np.array([sum(pr_all[i]) for i in np.arange(A)])
    entropy = (-1./A) * sum([np.log(pr_i[i]) for i in np.arange(A)])

    return entropy
