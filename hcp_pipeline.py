#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module implementing the main sequence of pipeline analysis to process
diffusion MRI from data distributed by Human Connectome Project (HCP).

"""

"""
The available steps of pipeline analysis:
1. Single shell extraction
2. Resampling data resolution
3. Computation of tensor model
"""

import os
import sys
import getopt
from hcp_parameters import *
from hcp_pipenode import single_shell_extraction, resample_data_resolution, compute_tensor_model, white_matter_mask_FA, white_matter_mask_wmparc, constrained_spherical_deconvolution, tracking_eudx, tracking_maxodf, compute_tract_query


do_step = [1] * 10
verbose = False


def run_pipeline():

    print "*** BEGINNING OF HCP PIPELINE COMPUTATION ***"

    if not os.path.isdir(dir_hcp_subj):
        print "FAIL: setting parameters - FILE: data directory not found!"
        sys.exit()

    subj = os.path.basename(os.path.realpath(dir_hcp_subj))
    step = 0

    print "Subject: ", subj
    print "Step %i: Setting parameters....." % step

    ## Setting directories

    dir_t1_src = os.path.join(dir_hcp_subj, dir_hcp_t1)
    dir_dwi_src = os.path.join(dir_hcp_subj, dir_hcp_dwi)
    dir_res_out = os.path.join(dir_hcp_subj, dir_hcp_out)
    if not os.path.exists(dir_res_out):
        os.makedirs(dir_res_out)

    print "DONE!"
    step += 1
    
    ## Preprocessing of HCP data

    print "Step %i: Shell extraction..." % step
    if do_step[step]:
        single_shell_extraction(dir_dwi_src, dir_res_out, verbose)
        print "DONE!"
    else:
        print "Skipped."
    step += 1

    print "Step %i: Resampling data resolution..." % step
    if do_step[step]:
        resample_data_resolution(dir_dwi_src, dir_res_out, verbose)
        print "DONE!"
    else:
        print "Skipped."
    step += 1

    print "Step %i: Computation of tensor model..." % step
    if do_step[step]:
        compute_tensor_model(dir_res_out, dir_res_out, verbose)
        print "DONE!"
    else:
        print "Skipped."
    step += 1

    print "Step %i: Computation of constrained spherical deconvolution..." % step
    if do_step[step]:
        constrained_spherical_deconvolution(dir_res_out, dir_res_out, verbose)
        print "DONE!"
    else:
        print "Skipped."
    step += 1

    print "Step %i: Computation white matter mask from FA..." % step
    if do_step[step]:
        white_matter_mask_FA(dir_res_out, dir_res_out, verbose)
        print "DONE!"
    else:
        print "Skipped."
    step += 1

    print "Step %i: Computation white matter mask from wmparc..." % step
    if do_step[step]:
        white_matter_mask_wmparc(dir_t1_src, dir_res_out, verbose)
        print "DONE!"
    else:
        print "Skipped."
    step += 1

    print "Step %i: Computation of EuDX Tracking..." % step
    if do_step[step]:
        tracking_eudx(dir_res_out, dir_res_out, verbose)
        print "DONE!"
    else:
        print "Skipped."
    step += 1

    print "Step %i: Computation of Max ODF Tracking ..." % step
    if do_step[step]:
        tracking_maxodf(dir_res_out, dir_res_out, verbose)
        print "DONE!"
    else:
        print "Skipped."
    step += 1

    print "Step %i: Tract Dissection by White Matter Query Language..." % step
    if do_step[step]:
        compute_tract_query(dir_res_out, dir_res_out, subj, verbose)
        print "DONE!"
    else:
        print "Skipped."
    step += 1

    print "*** END OF PIPELINE ***"


if __name__ == '__main__':
    
    for opt, arg in getopt.getopt(sys.argv[1:], "hi:s:v")[0]:
        
        if opt == '-s':
            do_step =   [0] * 12
            arg_step = map(int, arg.split())
            for s in arg_step: do_step[s]=1

        if opt == '-i' and os.path.isdir(arg):
            dir_hcp_subj = os.path.abspath(arg)

        if opt == '-p' and os.path.isfile(arg):
           print "load parameters not yet implemented"
           sys.exit()

        if opt == '-v':
            verbose = True
           
        if opt == '-h':
            print "Usage:"
            print "   pipeline.py <args>* "
            print "Arguments:"
            print "   path: -i <pathname>"
            print "         Directory where are located the data"
            print "   path: -p <pathname>"
            print "         File where are defined the parameters"
            print "   step: -s 'number number ...'"
            print "         1. Single shell extraction"
            print "         2. Resampling data resolution"
            print "         3. Computation of tensor model"
            print "         4. Computation of constrained spherical deconvolution"
            print "         5. Computation of white matter mask from FA"
            print "         6. Computation of white matter mask from wmparc"
            print "         7. Computation of EuDX Tracking"
            print "         8. Computation of Max ODF Tracking"
            print "         9. Tract Dissection by White Matter Query Language"
            print "   help: -h"
            print "         this help"
            print "Examples:"
            print "   pipeline.py -i /path/to/my/data"
            print "   pipeline.py -i /path/to/my/data -s '1 2 7'"
            print "   pipeline.py -s '1 2 7'"
            print "   pipeline.py -p /path/to/my/parameters"
            sys.exit()
            
    run_pipeline()
