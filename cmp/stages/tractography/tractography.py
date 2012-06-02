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
    dtb_cmd = '%s --dir %s --wm %s  --out %s %s' % (cmd, op.join(odf_out_path, 'dsi_dir.nii'),
                            # use the white matter mask after registration!
                            op.join(gconf.get_cmp_tracto_mask_tob0(), 'fsmask_1mm__8bit.nii'),
                            op.join(fibers_path, 'streamline.trk'), param )
    
    # set third argument to zero, in order to allow fast processing
    runCmd( dtb_cmd, log, 0.0 )
        
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
    
    # set third argument to zero, in order to allow fast processing
    runCmd( dtb_cmd, log, 0.0 )
        
    if not op.exists(op.join(fibers_path, 'streamline.trk')):
        log.error('No streamline.trk created')    
    
    log.info("[ DONE ]")

def fiber_tracking_dti():

    log.info("Run STREAMLINE tractography")
    log.info("===========================")
    
    fibers_path = gconf.get_cmp_fibers()
    odf_out_path = gconf.get_cmp_rawdiff_reconout()
    
    # streamline tractography
    # streamline tractography
    if not gconf.streamline_param == '':
        param = gconf.streamline_param
    else:
        param = '--angle 60  --seeds 32'
        
    cmd = op.join(gconf.get_cmp_binary_path(), 'DTB_streamline')
    dtb_cmd = '%s --dir %s --wm %s  --out %s %s' % (cmd, op.join(odf_out_path, 'dti_dir.nii'),
                            # use the white matter mask after registration!
                            op.join(gconf.get_cmp_tracto_mask_tob0(), 'fsmask_1mm__8bit.nii'),
                            op.join(fibers_path, 'streamline.trk'), param )
    # set third argument to zero, in order to allow fast processing
    runCmd( dtb_cmd, log, 0.0 )
        
    if not op.exists(op.join(fibers_path, 'streamline.trk')):
        log.error('No streamline.trk created')    
    
    log.info("[ DONE ]")

def probtrackx_tracking_dti():

    log.info("Run PROBTRACKX tractography")
    log.info("===========================")
    
    fibers_path = gconf.get_cmp_fibers()
    odf_out_path = gconf.get_cmp_rawdiff_reconout()

    # convert some imeges
    #convert_cmd = 'mri_convert %s/mri/orig.mgz %s/mri/nifti/orig.nii.gz' % (gconf.get_fs(),gconf.get_fs())
    #runCmd(convert_cmd, log)
    #convert_cmd = 'mri_convert %s/mri/rawavg.mgz %s/mri/nifti/rawavg.nii.gz' % (gconf.get_fs(),gconf.get_fs())
    #runCmd(convert_cmd, log)
    #convert_cmd = 'mri_convert %s/mri/brain.mgz %s/mri/nifti/brain.nii.gz' % (gconf.get_fs(),gconf.get_fs())
    #runCmd(convert_cmd, log)

    # convert surfaces
    convert_cmd = 'mris_convert %s/surf/lh.white %s/surf/lh.white.asc' % (gconf.get_fs(),gconf.get_fs())
    runCmd(convert_cmd, log)
    convert_cmd = 'mris_convert %s/surf/rh.white %s/surf/rh.white.asc' % (gconf.get_fs(),gconf.get_fs())
    runCmd(convert_cmd, log)
 
    # compute linear transform from Freesurfer to T1
    FS_to_T1_transform = op.join(gconf.get_nifti_trafo(), 'FS-TO-T1.mat')
    FS_to_T1_cmd = 'tkregister2 --mov %s/mri/orig.mgz --targ %s/mri/rawavg.mgz --regheader --reg /tmp/junk --fslregout %s --noedit' % (gconf.get_fs(),gconf.get_fs(), FS_to_T1_transform)
    runCmd(FS_to_T1_cmd, log)

    # concatenate transformations to obtain Freesurfer to DTI
    trafo_dir = gconf.get_nifti_trafo()
    nifti_dir = gconf.get_nifti()
    T1_to_FS_transform = op.join(gconf.get_nifti_trafo(), 'T1-TO-FS.mat')
    convert_cmd = 'convert_xfm -omat %s -inverse %s' % (T1_to_FS_transform,FS_to_T1_transform)
    runCmd(convert_cmd, log)

    FS_to_T2_transform = op.join(gconf.get_nifti_trafo(), 'FS-TO-T2.mat')
    convert_cmd = 'convert_xfm -omat %s -concat  %s %s' % (FS_to_T2_transform,op.join(gconf.get_nifti_trafo(),'T1-TO-T2.mat'), FS_to_T1_transform)
    runCmd(convert_cmd, log)

    T2_to_FS_transform = op.join(gconf.get_nifti_trafo(), 'T2-TO-FS.mat')
    convert_cmd = 'convert_xfm -omat %s -inverse %s' % (T2_to_FS_transform, FS_to_T2_transform)
    runCmd(convert_cmd, log)

    FS_to_b0_warp = op.join(gconf.get_nifti(), 'FS-TO-b0_warp.nii.gz')
    convert_cmd = 'convertwarp -o %s -r %s -m %s -w %s' % (FS_to_b0_warp, op.join(gconf.get_nifti(), 'DTI_first.nii.gz'), FS_to_T2_transform, op.join(gconf.get_nifti(),'T2-TO-b0_warp.nii.gz'))
    #runCmd(convert_cmd, log)
    b0_to_FS_warp = op.join(gconf.get_nifti(), 'b0-TO-FS_warp.nii.gz')
    convert_cmd = 'convertwarp -o %s -r %s/mri/fsmask_1mm.nii.gz -w %s --postmat=%s' % (b0_to_FS_warp, gconf.get_fs(), op.join(gconf.get_nifti(),'b0-TO-T2_warp.nii.gz'), T2_to_FS_transform)
    #runCmd(convert_cmd, log)

    # create surface labels
    if gconf.parcellation_scheme == 'Destrieux':
        labels_path = gconf.get_lausanne_parcellation_path('destrieuxaparc') + '/label'
        runCmd('mkdir -p %s' %(labels_path), log)
        convert_cmd = 'mri_annotation2label --subject FREESURFER --hemi lh --annotation %s/label/lh.aparc.a2009s.annot --outdir %s/label' % (gconf.get_fs(), gconf.get_lausanne_parcellation_path('destrieuxaparc'))
        runCmd(convert_cmd, log)
        convert_cmd = 'mri_annotation2label --subject FREESURFER --hemi rh --annotation %s/label/lh.aparc.a2009s.annot --outdir %s/label' % (gconf.get_fs(), gconf.get_lausanne_parcellation_path('destrieuxaparc'))
        runCmd(convert_cmd, log)
        rm_cmd = 'rm %s/?h.Unknown.label' % (gconf.get_lausanne_parcellation_path('destrieuxaparc'))
        runCmd(rm_cmd, log)
        targetvols_path = op.join(gconf.get_lausanne_parcellation_path('destrieuxaparc'),'targetvols')
        rm_mkdir_cmd = 'rm -rf %s; mkdir -p %s' % (targetvols_path, targetvols_path)
        runCmd(rm_mkdir_cmd, log)
    elif gconf.parcellation_scheme == 'NativeFreesurfer':
        labels_path = gconf.get_lausanne_parcellation_path('freesurferaparc') + '/label'
        runCmd('mkdir -p %s' %(labels_path), log)
        convert_cmd = 'mri_annotation2label --subject FREESURFER --hemi lh --annotation %s/label/lh.aparc.annot --outdir %s/label' % (gconf.get_fs(), gconf.get_lausanne_parcellation_path('freesurferaparc'))
        runCmd(convert_cmd, log)
        convert_cmd = 'mri_annotation2label --subject FREESURFER --hemi rh --annotation %s/label/lh.aparc.annot --outdir %s/label' % (gconf.get_fs(), gconf.get_lausanne_parcellation_path('freesurferaparc'))
        runCmd(convert_cmd, log)
        rm_cmd = 'rm %s/?h.Unknown.label' % (gconf.get_lausanne_parcellation_path('freesurferaparc'))
        runCmd(rm_cmd, log)
    else:
        log.error("Incompatible parcellation scheme: " + gconf.parcellation_scheme)
    
    identity_reg_cmd = 'printf "FREESURFER\n1\n1\n1\n1 0 0 0\n0 1 0 0\n0 0 1 0\n0 0 0 1\n" > %s/identitiy.dat' % (gconf.get_nifti_trafo())

    # convert labels of cortical labels to volumes
    pool = Pool(processes=gconf.nb_parallel_processes)
    roi_cmds = []
    for hemi in ['lh','rh']:
        for label in glob( op.join(labels_path, hemi + '*.label')):
            labelname = op.basename(label)
            roi_cmds.append('mri_label2vol --label %s --temp %s/mri/fsmask_1mm.nii.gz \
                                --subject FREESURFER --hemi %s --o %s/%s.nii.gz \
                                --proj frac 0 .5 0.1 --fillthresh 0.1 --reg ./register.dat; \
                             fslmaths %s/%s.nii.gz -sub %s/mri/fsmask_1mm.nii.gz \
                               -bin %s/%s.nii.gz' % (label, gconf.get_fs(), hemi, targetvols_path, labelname, targetvols_path, labelname, gconf.get_fs(), targetvols_path, labelname))
            
    for i in range(35,42) + range(76,84):
        if gconf.parcellation_scheme == 'NativeFreesurfer':
            parc = 'freesurferaparc'
        elif gconf.parcellation_scheme == 'Destrieux':
            parc = 'destrieuxaparc'
        roi_cmds.append('fslmaths %s/mri/fsmask_1mm.nii.gz  -kernel 3D -dilF \
               -mul %s/mri/ROIv_%s.nii.gz -thr %i \
               -uthr %i -bin ./targetvols/target%i.nii.gz' % (gconf.get_fs(),gconf.get_fs(),parc,i,i,i))

    result = pool.map(runCmdDefaultLog, roi_cmds)

    log.info("[ DONE ]")

def fiber_tracking_qball():

    log.info("Run STREAMLINE tractography")
    log.info("===========================")

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
        if conf.tractography_mode == 'streamline':
            fiber_tracking_dti()
        elif conf.tractography_mode == 'probabilistic':
            probtrackx_tracking_dti()
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
    diffusion_out_path = conf.get_cmp_rawdiff_reconout()
    
    conf.pipeline_status.AddStageInput(stage, conf.get_cmp_tracto_mask_tob0(), 'fsmask_1mm.nii.gz', 'fsmask_1mm-nii-gz')

    if conf.diffusion_imaging_model == 'DSI':
        conf.pipeline_status.AddStageInput(stage, diffusion_out_path, 'dsi_odf.nii', 'dsi_odf-nii')
    elif conf.diffusion_imaging_model == 'DTI':
        if conf.tractography_mode == 'streamline':
            conf.pipeline_status.AddStageInput(stage, diffusion_out_path, 'dti_tensor.nii', 'dti_tensor-nii')      
        elif conf.tractography_mode == 'probabilistic':
            conf.pipeline_status.AddStageInput(stage, diffusion_out_path, 'dyads1.nii.gz', 'dyads1-nii')
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
        conf.pipeline_status.AddStageOutput(stage, fibers_path, 'streamline.trk', 'streamline-trk')
    elif conf.diffusion_imaging_model == 'QBALL':
        conf.pipeline_status.AddStageOutput(stage, fibers_path, 'streamline.trk', 'streamline-trk')
          
