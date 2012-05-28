# Copyright (C) 2009-2011, Ecole Polytechnique Federale de Lausanne (EPFL) and
# Hospital Center and University of Lausanne (UNIL-CHUV), Switzerland
# All rights reserved.
#
#  This software is distributed under the open-source license Modified BSD.

import gzip
import os, os.path as op
from time import time
from ...logme import *
from glob import glob
import cmp.util as util
import nibabel as nib
import numpy as np
from multiprocessing import Pool

def runCmdDefaultLog(cmd):
    runCmd(cmd, log)

def convert_wm_mask():
    
    log.info("Convert WM MASK to 8 bit/pixel")
    log.info("==============================")
    
    infile = op.join(gconf.get_cmp_tracto_mask_tob0(), 'fsmask_1mm.nii.gz')
    outfile = op.join(gconf.get_cmp_tracto_mask_tob0(), 'fsmask_1mm__8bit.nii.gz')
    resampout = op.join(gconf.get_cmp_tracto_mask_tob0(), 'fsmask_1mm_resamp2x2x2.nii.gz')
    
    fsl_cmd = 'fslmaths %s %s -odt char' % (infile, outfile) 
    runCmd( fsl_cmd, log )
    
    # XXX: resample the white matter mask to b0 native space
    # make it usable dti_tracker
    mri_cmd = 'mri_convert -vs 2 2 2 %s %s ' % (infile, resampout)
    runCmd( mri_cmd, log )
        
    log.info("[ DONE ]")
    

def decompress(f):
    log.info("Decompress %s" % f)
    fin=gzip.open(f, 'rb')
    fc = fin.read()
    fout=open(f.rstrip('.gz'), 'wb')
    fout.write(fc)
    fout.close()
    fin.close()

def decompress_fsmask_nifti():
    log.info("Decompress Nifti files")
    log.info("======================")
    
    # need not decompress dtk output for now, store it
    # output as .nii
#    odf_out_path = gconf.get_cmp_rawdiff_reconout()
#    fi = glob(op.join(odf_out_path, '*.nii.gz'))
#    for f in fi:
#        decompress(f)
#        # remove .nii.gz
#        os.remove(f)

    # do not remove mask nii.gz
    mask = op.join(gconf.get_cmp_tracto_mask_tob0(), 'fsmask_1mm__8bit.nii.gz')
    decompress(mask)


def fiber_tracking_dsi():
    
    log.info("Run STREAMLINE tractography (DSI)")
    log.info("=================================")
    
    fibers_path = gconf.get_cmp_fibers()
    odf_out_path = gconf.get_cmp_rawdiff_reconout()
    
    # streamline tractography
    if not gconf.streamline_param == '':
        param = gconf.streamline_param
    else:
        param = '--angle 60 --seeds 32'

    cmd = op.join(gconf.get_cmp_binary_path(), 'DTB_streamline')
    dtb_cmd = '%s --dir %s --wm %s  --out %s %s' % (cmd, op.join(odf_out_path, 'dsi_dir.nii'),
                            # use the white matter mask after registration!
                            op.join(gconf.get_cmp_tracto_mask_tob0(), 'fsmask_1mm__8bit.nii'),
                            op.join(fibers_path, 'streamline.trk'), param )
    
    runCmd( dtb_cmd, log )
        
    if not op.exists(op.join(fibers_path, 'streamline.trk')):
        log.error('No streamline.trk created')    
    
    log.info("[ DONE ]")


    
def fiber_tracking_dsi_old_streamline():
    
    log.info("Run STREAMLINE tractography")
    log.info("===========================")
    
    fibers_path = gconf.get_cmp_fibers()
    odf_out_path = gconf.get_cmp_rawdiff_reconout()
    
    # streamline tractography
    if not gconf.streamline_param == '':
        param = gconf.streamline_param
    else:
        param = '--angle 60 --seeds 32'

    cmd = op.join(gconf.get_cmp_binary_path(), 'DTB_streamline')
    dtb_cmd = '%s --odf %s --wm %s --odfdir %s --out %s %s' % (cmd, op.join(odf_out_path, 'dsi_'),
                            # use the white matter mask after registration!
                            op.join(gconf.get_cmp_tracto_mask_tob0(), 'fsmask_1mm__8bit.nii'),
                            gconf.get_dtb_streamline_vecs_file(),
                            op.join(fibers_path, 'streamline'), param )
    
    runCmd( dtb_cmd, log )
        
    if not op.exists(op.join(fibers_path, 'streamline.trk')):
        log.error('No streamline.trk created')    
    
    log.info("[ DONE ]")

def fiber_tracking_dti():

    log.info("Run STREAMLINE tractography (DTI)")
    log.info("=================================")
    
    fibers_path = gconf.get_cmp_fibers()
    odf_out_path = gconf.get_cmp_rawdiff_reconout()
    
    # streamline tractography
    # streamline tractography
    if not gconf.streamline_param == '':
        param = gconf.streamline_param
    else:
        param = '--angle 60 --seeds 32'
        
    cmd = op.join(gconf.get_cmp_binary_path(), 'DTB_streamline')
    dtb_cmd = '%s --dir %s --wm %s --out %s %s' % (cmd, op.join(odf_out_path, 'dti_dir.nii'),
                            # use the white matter mask after registration!
                            op.join(gconf.get_cmp_tracto_mask_tob0(), 'fsmask_1mm__8bit.nii'),
                            op.join(fibers_path, 'streamline.trk'), param )
    # don't use runCmd, because otherwise slow like hell!
    #runCmd( dtb_cmd, log)
    os.system(dtb_cmd)
        
    if not op.exists(op.join(fibers_path, 'streamline.trk')):
        log.error('No streamline.trk created')
    
    log.info("[ DONE ]")

def fiber_tracking_probtrackx():

    log.info("Run probabilistic tractography")
    log.info("==============================")
    
    log.info(op.dirname(gconf.get_cmp_rawdiff()))
    
    bedpostxdir = gconf.get_cmp_rawdiff() + '.bedpostX'
    roidir = op.join(gconf.get_cmp(), 'fs_output', 'HR__registered-TO-b0',
                      'freesurferaparc')
    roifile = op.join(roidir, 'ROIv_HR_th.nii.gz')
    roifile_resamp = op.join(roidir, 'ROIv_HR_th_resampled.nii.gz')
    
    probtrack_dir = op.join(gconf.get_cmp(),'probtrack')
    
    mkdir_cmd = 'mkdir -p ' + probtrack_dir
    runCmd(mkdir_cmd,log)

    resampleROI_cmd = 'mri_convert -rl %s -rt nearest %s %s' % (op.join(bedpostxdir,'nodif_brain_mask.nii.gz'),
                                                                roifile,
                                                                roifile_resamp)
    runCmd(resampleROI_cmd, log)

    for roi in range(1,84):
        makeROIs_cmd = 'fslmaths %s -thr %s -uthr %s %s' % (roifile_resamp,str(roi),str(roi),op.join(roidir,'target' + str(roi) +'.nii.gz'))
        runCmd(makeROIs_cmd, log)

    rm_cmd = 'rm -f %s' % (op.join(probtrack_dir,'targetmasks'))
    runCmd(rm_cmd,log)

    createTargetMasksFile_cmd = 'touch %s' % (op.join(probtrack_dir,'targetmasks'))
    runCmd(createTargetMasksFile_cmd,log)
    for roi in range(1,84):
        fillTargetMaskFile_cmd = "echo %s/target%s.nii.gz >> %s" %(roidir,str(roi),op.join(probtrack_dir,'targetmasks'))
        runCmd(fillTargetMaskFile_cmd,log)

    pt_cmds = []
    for roi in range(1,84):
        # concatenate command
        pt_cmds.append('probtrackx -s %s -x %s  -l -c 0.2 -S 2000 --steplength=0.5 -P 5000 --forcedir --opd --mask=%s --os2t --mode=seedmask --targetmasks=%s --dir=%s --s2tastext' % (op.join(bedpostxdir,'merged'), op.join(roidir,'target'+str(roi)+'.nii.gz'),op.join(bedpostxdir,'nodif_brain_mask.nii.gz'),op.join(probtrack_dir,'targetmasks'),op.join(probtrack_dir,'seed'+str(roi))))

    # initialize pool
    pool = Pool(processes=4)
    results = pool.map(runCmdDefaultLog, pt_cmds)

    log.info("[ DONE ]")

def fiber_tracking_qball():

    log.info("Run STREAMLINE tractography (qball)")
    log.info("===================================")

    fibers_path = gconf.get_cmp_fibers()
    odf_out_path = gconf.get_cmp_rawdiff_reconout()

    # streamline tractography
    if not gconf.streamline_param == '':
        param = gconf.streamline_param
    else:
        param = '--angle 60  --seeds 32'

    cmd = op.join(gconf.get_cmp_binary_path(), 'DTB_streamline')
    dtb_cmd = '%s --dir %s --wm %s  --out %s %s' % (cmd, op.join(odf_out_path, 'hardi_dir.nii'),
                            # use the white matter mask after registration!
                            op.join(gconf.get_cmp_tracto_mask_tob0(), 'fsmask_1mm__8bit.nii'),
                            op.join(fibers_path, 'streamline.trk'), param )
    runCmd( dtb_cmd, log )

    if not op.exists(op.join(fibers_path, 'streamline.trk')):
        log.error('No streamline.trk created')

    log.info("[ DONE ]")
    
def inspect(gconf):
    """ Inspect the results of this stage """
    log = gconf.get_logger()
    trkcmd = 'trackvis %s' % op.join(gconf.get_cmp_fibers(), 'streamline.trk')
    runCmd( trkcmd, log )

def run(conf):
    """ Run the tractography step
    
    Parameters
    ----------
    conf : PipelineConfiguration object
        
    """
    # setting the global configuration variable
    globals()['gconf'] = conf
    globals()['log'] = gconf.get_logger() 
    start = time()
    
    convert_wm_mask()
    
    if gconf.diffusion_imaging_model == 'DSI':
        decompress_fsmask_nifti()
        fiber_tracking_dsi()
    elif gconf.diffusion_imaging_model == 'DTI':
        decompress_fsmask_nifti()
        if gconf.tracktography_mode == 'streamline':
            fiber_tracking_dti()
        if gconf.tracktography_mode == 'probabilistic':
            fiber_tracking_probtrackx()
    elif gconf.diffusion_imaging_model == 'QBALL':
        decompress_fsmask_nifti()
        fiber_tracking_qball()
    
    log.info("Module took %s seconds to process." % (time()-start))

    if not len(gconf.emailnotify) == 0:
        msg = ["Tractography", int(time()-start)]
        send_email_notification(msg, gconf, log)  

def declare_inputs(conf):
    """Declare the inputs to the stage to the PipelineStatus object"""
    
    stage = conf.pipeline_status.GetStage(__name__)
    if conf.tracktography_mode == 'streamline':
        diffusion_out_path = conf.get_cmp_rawdiff_reconout()  
    if conf.tracktography_mode == 'probabilistic':
        diffusion_out_path = conf.get_cmp_rawdiff() + '.bedpostX'
    conf.pipeline_status.AddStageInput(stage, conf.get_cmp_tracto_mask_tob0(), 'fsmask_1mm.nii.gz', 'fsmask_1mm-nii-gz')

    if conf.diffusion_imaging_model == 'DSI':
        conf.pipeline_status.AddStageInput(stage, diffusion_out_path, 'dsi_odf.nii', 'dsi_odf-nii')
    elif conf.diffusion_imaging_model == 'DTI':
        if conf.tracktography_mode == 'streamline':
            conf.pipeline_status.AddStageInput(stage, diffusion_out_path, 'dti_tensor.nii', 'dti_tensor-nii')
        if conf.tracktography_mode == 'probabilistic':
            conf.pipeline_status.AddStageInput(stage, diffusion_out_path, 'merged_th1samples.nii.gz', 'merged_th1samples-nii-gz')  
    elif conf.diffusion_imaging_model == 'QBALL':
        conf.pipeline_status.AddStageInput(stage, diffusion_out_path, 'hardi_odf.nii', 'hardi_odf-nii')
    
def declare_outputs(conf):
    """Declare the outputs to the stage to the PipelineStatus object"""
    
    stage = conf.pipeline_status.GetStage(__name__)
    fibers_path = conf.get_cmp_fibers()
        
    conf.pipeline_status.AddStageOutput(stage, conf.get_cmp_tracto_mask_tob0(), 'fsmask_1mm__8bit.nii.gz', 'fsmask_1mm__8bit-nii-gz')
    
    if conf.diffusion_imaging_model == 'DSI':
        conf.pipeline_status.AddStageOutput(stage, fibers_path, 'streamline.trk', 'streamline-trk')
    elif conf.diffusion_imaging_model == 'DTI':
        if conf.tracktography_mode == 'streamline':
            conf.pipeline_status.AddStageOutput(stage, fibers_path, 'streamline.trk', 'streamline-trk')
    elif conf.diffusion_imaging_model == 'QBALL':
        conf.pipeline_status.AddStageOutput(stage, fibers_path, 'streamline.trk', 'streamline-trk')
          
