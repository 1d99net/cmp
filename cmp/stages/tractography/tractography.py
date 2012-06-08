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
    runCmd(cmd, log, 0.05)

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

    pool = Pool(processes=gconf.nb_parallel_processes)
    convertwarp_cmds = []
    if gconf.registration_mode == 'Nonlinear':
        FS_to_b0_warp = op.join(gconf.get_nifti(), 'FS-TO-b0_warp.nii.gz')
        convertwarp_cmds.append('convertwarp -o %s -r %s -m %s -w %s' % (FS_to_b0_warp, op.join(gconf.get_nifti(), 'DTI_first.nii.gz'), FS_to_T2_transform, op.join(gconf.get_nifti(),'T2-TO-b0_warp.nii.gz')))
        b0_to_FS_warp = op.join(gconf.get_nifti(), 'b0-TO-FS_warp.nii.gz')
        convertwarp_cmds.append('convertwarp -o %s -r %s/mri/fsmask_1mm.nii.gz -w %s --postmat=%s' % (b0_to_FS_warp, gconf.get_fs(), op.join(gconf.get_nifti(),'b0-TO-T2_warp.nii.gz'), T2_to_FS_transform))

    result = pool.map(runCmdDefaultLog, convertwarp_cmds)
    
    # create surface labels
    if gconf.parcellation_scheme == 'Destrieux':
        roi_path = op.join(gconf.get_cmp_tracto_mask(),'destrieuxaparc')
        labels_path = op.join(roi_path, 'label')
        runCmd('mkdir -p %s' %(labels_path), log)
        convert_cmd = 'mri_annotation2label --sd %s --subject FREESURFER --hemi lh --annotation %s/label/lh.aparc.a2009s.annot --outdir %s' % (gconf.get_subj_dir(),gconf.get_fs(), labels_path)
        runCmd(convert_cmd, log)
        convert_cmd = 'mri_annotation2label --sd %s --subject FREESURFER --hemi rh --annotation %s/label/rh.aparc.a2009s.annot --outdir %s' % (gconf.get_subj_dir(),gconf.get_fs(), labels_path)
        runCmd(convert_cmd, log)
        rm_cmd = 'rm %s/?h.Unknown.label' % (labels_path)
        runCmd(rm_cmd, log)
        targetvols_path = op.join(gconf.get_cmp_tracto_mask(),'destrieuxaparc','targetvols')
        rm_mkdir_cmd = 'mkdir -p %s' % (targetvols_path)
        runCmd(rm_mkdir_cmd, log)

    elif gconf.parcellation_scheme == 'NativeFreesurfer':
        roi_path = op.join(gconf.get_cmp_tracto_mask(),'freesurferaparc')
        labels_path = op.join(roi_path, 'label')
        runCmd('mkdir -p %s' %(labels_path), log)
        convert_cmd = 'mri_annotation2label --subject FREESURFER --hemi lh --annotation %s/label/lh.aparc.annot --outdir %s' % (gconf.get_fs(), labels_path)
        runCmd(convert_cmd, log)
        convert_cmd = 'mri_annotation2label --subject FREESURFER --hemi rh --annotation %s/label/lh.aparc.annot --outdir %s' % (gconf.get_fs(), labels_path)
        runCmd(convert_cmd, log)
        rm_cmd = 'rm %s/?h.Unknown.label' % (labels_path)
        runCmd(rm_cmd, log)
        targetvols_path = op.join(gconf.get_cmp_tracto_mask(),'freesurferaparc','targetvols') 
        rm_mkdir_cmd = 'mkdir -p %s' % (targetvols_path)
        runCmd(rm_mkdir_cmd, log)
    else:
        log.error("Incompatible parcellation scheme: " + gconf.parcellation_scheme)
    
    identity_reg_cmd = 'printf "FREESURFER\n1\n1\n1\n1 0 0 0\n0 1 0 0\n0 0 1 0\n0 0 0 1\n" > %s/identity.dat' % (gconf.get_nifti_trafo())
    runCmd(identity_reg_cmd, log)

    # convert labels of cortical labels to volumes
    pool = Pool(processes=gconf.nb_parallel_processes)
    roi_cmds = []
    for hemi in ['lh','rh']:
        for label in glob( op.join(labels_path, hemi + '*.label')):
            labelname = op.basename(label)
            roi_cmds.append('SUBJECTS_DIR=%s; mri_label2vol --label %s --temp %s/mri/fsmask_1mm.nii.gz \
                                --subject FREESURFER --hemi %s --o %s/%s.nii.gz \
                                --proj frac 0 .5 0.1 --fillthresh 0.1 --reg %s/identity.dat; \
                             fslmaths %s/%s.nii.gz -sub %s/mri/fsmask_1mm.nii.gz \
                               -bin %s/%s.nii.gz' % (gconf.get_subj_dir(), label, gconf.get_fs(), hemi, targetvols_path, labelname, gconf.get_nifti_trafo(), targetvols_path, labelname, gconf.get_fs(), targetvols_path, labelname))
    
    if gconf.parcellation_scheme == 'NativeFreesurfer':
        parc = 'freesurferaparc'
        rng = range(35,42) + range(76,84)
    elif gconf.parcellation_scheme == 'Destrieux':
        parc = 'destrieuxaparc'
        rng = range(75,82) + range(156,164)

        
    for i in rng:
        roi_cmds.append('fslmaths %s/mri/fsmask_1mm.nii.gz  -kernel 3D -dilF \
               -mul %s/mri/ROIv_%s.nii.gz -thr %i \
               -uthr %i -bin %s/target%i.nii.gz' % (gconf.get_fs(),gconf.get_fs(),parc,i,i,targetvols_path,i))

    result = pool.map(runCmdDefaultLog, roi_cmds)

    # union of ROIs
    rois = glob(op.join(targetvols_path, '*.nii.gz'))
    roi_union_cmd = 'fslmaths %s -bin %s' % (' -add '.join(rois), op.join(roi_path,'ROI_union.nii.gz'))
    runCmd(roi_union_cmd, log)
    avoid_mask_cmd = 'fslmaths %s/mri/fsmask_1mm.nii.gz -add %s -bin -mul -1 -add 1 -bin %s' % (gconf.get_fs(), op.join(roi_path, 'ROI_union.nii.gz'), op.join(roi_path, 'fsmask_1mm_avoid.nii.gz'))
    runCmd(avoid_mask_cmd, log)
    waypoint_mask_cmd = 'fslmaths %s/mri/fsmask_1mm.nii.gz -bin -kernel 3D -dilM -bin %s' % (gconf.get_fs(), op.join(roi_path, 'fsmask_1mm_waypoint.nii.gz'))
    runCmd(waypoint_mask_cmd, log)

    if gconf.parcellation_scheme == 'Destrieux':
        fin = open(op.join(gconf.get_lausanne_parcellation_path('destrieuxaparc'),'targets.txt'))

    elif gconf.parcellation_scheme == 'NativeFreesurfer':
        fin = open(op.join(gconf.get_lausanne_parcellation_path('freesurferaparc'),'targets.txt'))

    tracto_targets = fin.read().split()
    fin.close()

    fout = open(targetvols_path + '.txt', 'w')
    for target in tracto_targets:
        fout.write(op.join(targetvols_path,target + '.nii.gz\n'))
        
    fout.close()

    # Tractography
    if gconf.registration_mode == 'BBregister':
        xfm = op.join(gconf.get_nifti_bbregister(),'orig-TO-b0.mat')
        invxfm = op.join(gconf.get_nifti_bbregister(), 'b0-TO-orig.mat')
    elif gconf.registration_mode == 'Nonlinear':
        xfm = op.join(gconf.get_nifti(),'FS-TO-b0_warp.nii.gz')
        invxfm = op.join(gconf.get_nifti(), 'b0-TO-FS_warp.nii.gz')        
    elif gconf.registration_mode == 'Linear':
        xfm = op.join(gconf.get_nifti_trafo(),'FS-TO-b0.mat')
        invxfm = op.join(gconf.get_nifti_trafo(), 'b0-TO-FS.mat')
    else:
        log.error('incompatible registration method: %s' % (gconf.registration_mode))

    pool = Pool(processes=gconf.nb_parallel_processes)        
    probtrackx_cmds = []
    if not op.exists(op.join(gconf.get_fs(),'tmp')):
        os.mkdir(op.join(gconf.get_fs(),'tmp'))

    probtrackx_cmds_tmp = []
    for hemi in ['lh','rh']:
        for label in glob(op.join(labels_path,hemi + '.*.label')):
            stopmask = op.join(gconf.get_fs(),'tmp',label + '_stop.nii.gz')
            probtrackx_cmds.append('fslmaths %s -sub %s %s; \
                                    probtrackx --samples=%s \
                                               --mask=%s \
                                               --seed=%s \
                                               --verbose=1 \
                                               --mode=seedmask \
                                               --targetmasks=%s \
                                               --mesh=%s \
                                               --seedref=%s \
                                               --dir=%s \
                                               --forcedir --opd --os2t --loopcheck \
                                               --out=fdt_paths.nii.gz \
                                               --avoid=%s \
                                               --waypoints=%s \
                                               --xfm=%s --invxfm=%s \
                                               --nsamples=%i --nsteps=%i --distthresh=%i --cthr=%f --steplength=%f \
                                               --s2tastext \
                                               --stop=%s %s' % (op.join(roi_path,'ROI_union.nii.gz'),
                                                                op.join(targetvols_path,op.basename(label) + '.nii.gz'),
                                                                stopmask,
                                                                op.join(gconf.get_cmp_rawdiff_reconout(),'merged'),
                                                                op.join(gconf.get_cmp_rawdiff_reconout(),'nodif_brain_mask.nii.gz'),
                                                                label,
                                                                targetvols_path + '.txt',
                                                                op.join(gconf.get_fs(),'surf',hemi + '.white.asc'),
                                                                op.join(gconf.get_fs(),'mri','fsmask_1mm.nii.gz'),
                                                                op.join(gconf.get_cmp(),'probtractography',op.basename(label)),
                                                                op.join(roi_path, 'fsmask_1mm_avoid.nii.gz'),
                                                                op.join(roi_path, 'fsmask_1mm_waypoint.nii.gz'),
                                                                xfm,
                                                                invxfm,
                                                                int(gconf.probtrackx_options_nsamples),
                                                                int(gconf.probtrackx_options_nsteps),
                                                                int(gconf.probtrackx_options_distthresh),
                                                                float(gconf.probtrackx_options_cthr),
                                                                float(gconf.probtrackx_options_steplength),
                                                                stopmask,
                                                                gconf.probtrackx_options_other))
            
    for seedroi in rng:
        stopmask = op.join(gconf.get_fs(),'tmp','target' +  str(seedroi) + '_stop.nii.gz')
        probtrackx_cmds.append('fslmaths %s -sub %s %s; \
                                probtrackx --samples=%s \
                                           --mask=%s \
                                           --seed=%s \
                                           --verbose=1 \
                                           --mode=seedmask \
                                           --targetmasks=%s \
                                           --seedref=%s \
                                           --dir=%s \
                                           --forcedir --opd --os2t --loopcheck \
                                           --out=fdt_paths.nii.gz \
                                           --avoid=%s \
                                           --waypoints=%s \
                                           --xfm=%s --invxfm=%s \
                                           --nsamples=%i --nsteps=%i --distthresh=%i --cthr=%f --steplength=%f \
                                           --s2tastext \
                                           --stop=%s %s' % (op.join(roi_path,'ROI_union.nii.gz'),
                                                            op.join(targetvols_path,'target' + str(seedroi) + '.nii.gz'),
                                                            stopmask,
                                                            op.join(gconf.get_cmp_rawdiff_reconout(),'merged'),
                                                            op.join(gconf.get_cmp_rawdiff_reconout(),'nodif_brain_mask.nii.gz'),
                                                            op.join(targetvols_path,'target' + str(seedroi) + '.nii.gz'),
                                                            targetvols_path + '.txt',
                                                            op.join(gconf.get_fs(),'mri','fsmask_1mm.nii.gz'),
                                                            op.join(gconf.get_cmp(),'probtractography','target' + str(seedroi)),
                                                            op.join(roi_path, 'fsmask_1mm_avoid.nii.gz'),
                                                            op.join(roi_path, 'fsmask_1mm_waypoint.nii.gz'),
                                                            xfm, invxfm, 
                                                            int(gconf.probtrackx_options_nsamples),
                                                            int(gconf.probtrackx_options_nsteps),
                                                            int(gconf.probtrackx_options_distthresh),
                                                            float(gconf.probtrackx_options_cthr),
                                                            float(gconf.probtrackx_options_steplength),
                                                            stopmask,
                                                            gconf.probtrackx_options_other))

    result = pool.map(runCmdDefaultLog, probtrackx_cmds)

    # construct connectivity matrix
    if gconf.parcellation_scheme == 'Destrieux':
        conmatrix = op.join(gconf.get_cmp(),'probconmatrix_destrieuxaparc.txt')
        fin = open(targetvols_path + '.txt')
        numregions = 163
    elif gconf.parcellation_scheme == 'NativeFreesurfer':
        connmatrix = op.join(gconf.get_cmp(),'probconmatrix_freesurferaparc.txt')
        fin = open(targetvols_path + '.txt')
        numregiens = 83
    elif gconf.parcellation_scheme == 'Lausanne2008':
        log.error("Parcellation Lausanne2008 non supported for probabilistic tractography.")
    else:
        log.error("Parcellation scheme not recognized:" + gconf.parcellation_scheme)

    tracto_targets = fin.read().split()
    fin.close()

    matrix = []
    for seed in tracto_targets:
        s2t_file = op.basename(seed)[:-7]
        s2t_matrix = np.loadtxt(op.join(gconf.get_cmp(),'probtractography',s2t_file,'matrix_seeds_to_all_targets'))
        s2t_sum = np.sum(s2t_matrix,axis=0)
        if matrix == []:
            matrix = s2t_sum
        else:
            matrix = np.vstack((matrix,s2t_sum))

        matrix.shape

    np.savetxt(conmatrix, matrix)
        
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
        if conf.tractography_mode == 'streamline':
            conf.pipeline_status.AddStageOutput(stage, fibers_path, 'streamline.trk', 'streamline-trk')
        elif conf.tractography_mode == 'probabilistic':
            if conf.parcellation_scheme == 'Destrieux':
                conf.pipeline_status.AddStageOutput(stage, conf.get_cmp(), 'probconmatrix_destrieuxaparc.txt', 'probconmatrix')
            elif cconf.parcellation_scheme == 'NativeFreesurfer':
                conf.pipeline_status.AddStageOutput(stage, conf.get_cmp(), 'probconmatrix_freesurferaparc.txt', 'probconmatrix')
    elif conf.diffusion_imaging_model == 'QBALL':
        conf.pipeline_status.AddStageOutput(stage, fibers_path, 'streamline.trk', 'streamline-trk')
          
