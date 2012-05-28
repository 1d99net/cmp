# Copyright (C) 2009-2011, Ecole Polytechnique Federale de Lausanne (EPFL) and
# Hospital Center and University of Lausanne (UNIL-CHUV), Switzerland
# All rights reserved.
#
#  This software is distributed under the open-source license Modified BSD.

import os, os.path as op
from glob import glob
from time import time
import shutil
import subprocess
import sys
from ...logme import *

import nibabel as ni
import networkx as nx
import numpy as np
import math
from cmp.util import mymove

def create_annot_label():

    log.info("Create the cortical labels necessary for our ROIs")
    log.info("=================================================")

    subjects_dir = gconf.get_subj_dir()

    fs_label_dir = op.join(gconf.get_fs(), 'label')
    fs_dir = gconf.get_fs()

    paths = []

    for p in gconf.parcellation.keys():
        for hemi in ['lh', 'rh']:
            spath = gconf.parcellation[p]['fs_label_subdir_name'] % hemi
            paths.append(spath)

    # XXX: this makedirs should go in the preprocessing module
    for p in paths:
        try:
            os.makedirs(op.join(fs_label_dir, p))
        except:
            pass

    comp = [
    ('rh', 'myatlas_36_rh.gcs', 'rh.myaparc_36.annot', 'regenerated_rh_36','myaparc_36'),
    ('rh', 'myatlasP1_16_rh.gcs','rh.myaparcP1_16.annot','regenerated_rh_500','myaparcP1_16'),
    ('rh', 'myatlasP17_28_rh.gcs','rh.myaparcP17_28.annot','regenerated_rh_500','myaparcP17_28'),
    ('rh', 'myatlasP29_36_rh.gcs','rh.myaparcP29_36.annot','regenerated_rh_500','myaparcP29_36'),
    ('rh','myatlas_60_rh.gcs','rh.myaparc_60.annot','regenerated_rh_60','myaparc_60'),
    ('rh','myatlas_125_rh.gcs','rh.myaparc_125.annot','regenerated_rh_125','myaparc_125'),
    ('rh','myatlas_250_rh.gcs','rh.myaparc_250.annot','regenerated_rh_250','myaparc_250'),
    ('lh', 'myatlas_36_lh.gcs', 'lh.myaparc_36.annot', 'regenerated_lh_36','myaparc_36'),
    ('lh', 'myatlasP1_16_lh.gcs','lh.myaparcP1_16.annot','regenerated_lh_500','myaparcP1_16'),
    ('lh', 'myatlasP17_28_lh.gcs','lh.myaparcP17_28.annot','regenerated_lh_500','myaparcP17_28'),
    ('lh', 'myatlasP29_36_lh.gcs','lh.myaparcP29_36.annot','regenerated_lh_500','myaparcP29_36'),
    ('lh','myatlas_60_lh.gcs','lh.myaparc_60.annot','regenerated_lh_60', 'myaparc_60'),
    ('lh','myatlas_125_lh.gcs','lh.myaparc_125.annot','regenerated_lh_125','myaparc_125'),
    ('lh','myatlas_250_lh.gcs','lh.myaparc_250.annot','regenerated_lh_250','myaparc_250'),
    ]

    for out in comp:      
        mris_cmd = 'mris_ca_label "FREESURFER" %s "%s/surf/%s.sphere.reg" "%s" "%s" ' % (out[0], 
            fs_dir, out[0], gconf.get_lausanne_atlas(out[1]), op.join(fs_label_dir, out[2]))
        runCmd( mris_cmd, log )
        log.info('-----------')

        annot = '--annotation "%s"' % out[4]

        mri_an_cmd = 'mri_annotation2label --subject "FREESURFER" --hemi %s --outdir "%s" %s' % (out[0], op.join(fs_label_dir, out[3]), annot)
        runCmd( mri_an_cmd, log )
        log.info('-----------')

    # extract cc and unknown to add to tractography mask, we do not want this as a region of interest
    # in FS 5.0, unknown and corpuscallosum are not available for the 35 scale (why?),
    # but for the other scales only, take the ones from _60
    rhun = op.join(fs_label_dir, 'rh.unknown.label')
    lhun = op.join(fs_label_dir, 'lh.unknown.label')
    rhco = op.join(fs_label_dir, 'rh.corpuscallosum.label')
    lhco = op.join(fs_label_dir, 'lh.corpuscallosum.label')
    shutil.copy(op.join(fs_label_dir, 'regenerated_rh_60', 'rh.unknown.label'), rhun)
    shutil.copy(op.join(fs_label_dir, 'regenerated_lh_60', 'lh.unknown.label'), lhun)
    shutil.copy(op.join(fs_label_dir, 'regenerated_rh_60', 'rh.corpuscallosum.label'), rhco)
    shutil.copy(op.join(fs_label_dir, 'regenerated_lh_60', 'lh.corpuscallosum.label'), lhco)

    mri_cmd = """mri_label2vol --label "%s" --label "%s" --label "%s" --label "%s" --temp "%s" --o  "%s" --identity """ % (rhun, lhun, rhco, lhco, op.join(fs_dir, 'mri', 'orig.mgz'), op.join(fs_dir, 'label', 'cc_unknown.nii.gz') )
    runCmd( mri_cmd, log )

    runCmd( 'mris_volmask "FREESURFER"', log)

    mri_cmd = 'mri_convert -i "%s/mri/ribbon.mgz" -o "%s/mri/ribbon.nii.gz"' % (fs_dir, fs_dir)
    runCmd( mri_cmd, log )

    mri_cmd = 'mri_convert -i "%s/mri/aseg.mgz" -o "%s/mri/aseg.nii.gz"' % (fs_dir, fs_dir)
    runCmd( mri_cmd, log )

    log.info("[ DONE ]")  


def create_roi():
    """ Creates the ROI_%s.nii.gz files using the given parcellation information
    from networks. Iteratively create volume. """

    log.info("Create the ROIs:")
    fs_dir = gconf.get_fs()

    # load aseg volume
    aseg = ni.load(op.join(fs_dir, 'mri', 'aseg.nii.gz'))
    asegd = aseg.get_data()	# numpy.ndarray

    # identify cortical voxels, right (3) and left (42) hemispheres
    idxr = np.where(asegd == 3)
    idxl = np.where(asegd == 42)
    xx = np.concatenate((idxr[0],idxl[0]))
    yy = np.concatenate((idxr[1],idxl[1]))
    zz = np.concatenate((idxr[2],idxl[2]))

    # initialize variables necessary for cortical ROIs dilation
    # dimensions of the neighbourhood for rois labels assignment (choose odd dimensions!)
    shape = (25,25,25)
    center = np.array(shape) // 2
    # dist: distances from the center of the neighbourhood
    dist = np.zeros(shape, dtype='float32')
    for x in range(shape[0]):
        for y in range(shape[1]):
            for z in range(shape[2]):
                distxyz = center - [x,y,z]
                dist[x,y,z] = math.sqrt(np.sum(np.multiply(distxyz,distxyz)))

    # LOOP throughout all the SCALES
    # (from the one with the highest number of region to the one with the lowest number of regions)
    parkeys = gconf.parcellation.keys()
    values = list()
    for i in range(len(parkeys)):
        values.append(gconf.parcellation[parkeys[i]]['number_of_regions'])
    temp = zip(values, parkeys)
    temp.sort(reverse=True)
    values, parkeys = zip(*temp)
    roisMax = np.zeros( (256, 256, 256), dtype=np.int16 ) # numpy.ndarray
    for i,parkey in enumerate(parkeys):
        parval = gconf.parcellation[parkey]

        log.info("Working on parcellation: " + parkey)
        log.info("========================")
        pg = nx.read_graphml(parval['node_information_graphml'])

        # each node represents a brain region
        # create a big 256^3 volume for storage of all ROIs
        rois = np.zeros( (256, 256, 256), dtype=np.int16 ) # numpy.ndarray

        for brk, brv in pg.nodes_iter(data=True):   # slow loop

            if brv['dn_hemisphere'] == 'left':
                hemi = 'lh'
            elif brv['dn_hemisphere'] == 'right':
                hemi = 'rh'

            if brv['dn_region'] == 'subcortical':

                log.info("---------------------")
                log.info("Work on brain region: %s" % (brv['dn_region']) )
                log.info("Freesurfer Name: %s" %  brv['dn_fsname'] )
                log.info("---------------------")

                # if it is subcortical, retrieve roi from aseg
                idx = np.where(asegd == int(brv['dn_fs_aseg_val']))
                rois[idx] = int(brv['dn_correspondence_id'])

            elif brv['dn_region'] == 'cortical':
                log.info("---------------------")
                log.info("Work on brain region: %s" % (brv['dn_region']) )
                log.info("Freesurfer Name: %s" %  brv['dn_fsname'] )
                log.info("---------------------")

                labelpath = op.join(fs_dir, 'label', parval['fs_label_subdir_name'] % hemi)

                # construct .label file name
                fname = '%s.%s.label' % (hemi, brv['dn_fsname'])

                # execute fs mri_label2vol to generate volume roi from the label file
                # store it in temporary file to be overwritten for each region (slow!)
                mri_cmd = 'mri_label2vol --label "%s" --temp "%s" --o "%s" --identity' % (op.join(labelpath, fname),
                        op.join(fs_dir, 'mri', 'orig.mgz'), op.join(labelpath, 'tmp.nii.gz'))
                runCmd( mri_cmd, log )

                tmp = ni.load(op.join(labelpath, 'tmp.nii.gz'))
                tmpd = tmp.get_data()

                # find voxel and set them to intensity value in rois
                idx = np.where(tmpd == 1)
                rois[idx] = int(brv['dn_correspondence_id'])

        newrois = rois.copy()
        # store scale500 volume for correction on multi-resolution consistency
        if i == 0:
            log.info("Storing ROIs volume maximal resolution...")
            roisMax = rois.copy()
            idxMax = np.where(roisMax > 0)
            xxMax = idxMax[0]
            yyMax = idxMax[1]
            zzMax = idxMax[2]
        # correct cortical surfaces using as reference the roisMax volume (for consistency between resolutions)
        else:
            log.info("Adapt cortical surfaces...")
            adaptstart = time()
            idxRois = np.where(rois > 0)
            xxRois = idxRois[0]
            yyRois = idxRois[1]
            zzRois = idxRois[2]
            # correct voxels labeled in current resolution, but not labeled in highest resolution
            for j in range(xxRois.size):
                if ( roisMax[xxRois[j],yyRois[j],zzRois[j]]==0 ):
                    newrois[xxRois[j],yyRois[j],zzRois[j]] = 0;
            # correct voxels not labeled in current resolution, but labeled in highest resolution
            for j in range(xxMax.size):
                if ( newrois[xxMax[j],yyMax[j],zzMax[j]]==0 ):
                    local = extract(rois, shape, position=(xxMax[j],yyMax[j],zzMax[j]), fill=0)
                    mask = local.copy()
                    mask[np.nonzero(local>0)] = 1
                    thisdist = np.multiply(dist,mask)
                    thisdist[np.nonzero(thisdist==0)] = np.amax(thisdist)
                    value = np.int_(local[np.nonzero(thisdist==np.amin(thisdist))])
                    if value.size > 1:
                        counts = np.bincount(value)
                        value = np.argmax(counts)
                    newrois[xxMax[j],yyMax[j],zzMax[j]] = value
            log.info("Cortical ROIs adaptation took %s seconds to process." % (time()-adaptstart))

        # store volume eg in ROI_scale33.nii.gz
        out_roi = op.join(fs_dir, 'label', 'ROI_%s.nii.gz' % parkey)
        # update the header
        hdr = aseg.get_header()
        hdr2 = hdr.copy()
        hdr2.set_data_dtype(np.uint16)
        log.info("Save output image to %s" % out_roi)
        img = ni.Nifti1Image(newrois, aseg.get_affine(), hdr2)
        ni.save(img, out_roi)

        # dilate cortical regions
        log.info("Dilating cortical regions...")
        dilatestart = time()
        # loop throughout all the voxels belonging to the aseg GM volume
        for j in range(xx.size):
            if newrois[xx[j],yy[j],zz[j]] == 0:
                local = extract(rois, shape, position=(xx[j],yy[j],zz[j]), fill=0)
                mask = local.copy()
                mask[np.nonzero(local>0)] = 1
                thisdist = np.multiply(dist,mask)
                thisdist[np.nonzero(thisdist==0)] = np.amax(thisdist)
                value = np.int_(local[np.nonzero(thisdist==np.amin(thisdist))])
                if value.size > 1:
                    counts = np.bincount(value)
                    value = np.argmax(counts)
                newrois[xx[j],yy[j],zz[j]] = value
        log.info("Cortical ROIs dilation took %s seconds to process." % (time()-dilatestart))

        # store volume eg in ROIv_scale33.nii.gz
        out_roi = op.join(fs_dir, 'label', 'ROIv_%s.nii.gz' % parkey)
        log.info("Save output image to %s" % out_roi)
        img = ni.Nifti1Image(newrois, aseg.get_affine(), hdr2)
        ni.save(img, out_roi)

    log.info("[ DONE ]")  


def create_wm_mask():

    log.info("Create white matter mask")

    fs_dir = gconf.get_fs()
    fs_cmd_dir = gconf.get_cmp_fsout()
    reg_path = gconf.get_cmp_tracto_mask()

    # load ribbon as basis for white matter mask
    fsmask = ni.load(op.join(fs_dir, 'mri', 'ribbon.nii.gz'))
    fsmaskd = fsmask.get_data()

    wmmask = np.zeros( fsmask.get_data().shape )

    # these data is stored and could be extracted from fs_dir/stats/aseg.txt

    # extract right and left white matter 
    idx_lh = np.where(fsmaskd == 120)
    idx_rh = np.where(fsmaskd == 20)

    wmmask[idx_lh] = 1
    wmmask[idx_rh] = 1

    # remove subcortical nuclei from white matter mask
    aseg = ni.load(op.join(fs_dir, 'mri', 'aseg.nii.gz'))
    asegd = aseg.get_data()

    try:
        import scipy.ndimage.morphology as nd
    except ImportError:
        raise Exception('Need scipy for binary erosion of white matter mask')

    # need binary erosion function
    imerode = nd.binary_erosion

    # ventricle erosion    
    csfA = np.zeros( asegd.shape )
    csfB = np.zeros( asegd.shape )

    # structuring elements for erosion
    se1 = np.zeros( (3,3,5) )
    se1[1,:,2] = 1; se1[:,1,2] = 1; se1[1,1,:] = 1
    se = np.zeros( (3,3,3) )
    se[1,:,1] = 1; se[:,1,1] = 1; se[1,1,:] = 1

    # lateral ventricles, thalamus proper and caudate
    # the latter two removed for better erosion, but put back afterwards
    idx = np.where( (asegd == 4) |
                    (asegd == 43) |
                    (asegd == 11) |
                    (asegd == 50) |
                    (asegd == 31) |
                    (asegd == 63) |
                    (asegd == 10) |
                    (asegd == 49) )
    csfA[idx] = 1
    csfA = imerode(imerode(csfA, se1),se)

    # thalmus proper and cuadate are put back because they are not lateral ventricles
    idx = np.where( (asegd == 11) |
                    (asegd == 50) |
                    (asegd == 10) |
                    (asegd == 49) )
    csfA[idx] = 0

    # REST CSF, IE 3RD AND 4TH VENTRICULE AND EXTRACEREBRAL CSF    
    idx = np.where( (asegd == 5) |
                    (asegd == 14) |
                    (asegd == 15) |
                    (asegd == 24) |
                    (asegd == 44) |
                    (asegd == 72) |
                    (asegd == 75) |
                    (asegd == 76) |
                    (asegd == 213) |
                    (asegd == 221))    
    # 43 ??, 4??  213?, 221?
    # more to discuss.
    for i in [5,14,15,24,44,72,75,76,213,221]:
        idx = np.where(asegd == i)
        csfB[idx] = 1

    # do not remove the subthalamic nucleus for now from the wm mask
    # 23, 60
    # would stop the fiber going to the segmented "brainstem"

    # grey nuclei, either with or without erosion
    gr_ncl = np.zeros( asegd.shape )

    # with erosion
    for i in [10,11,12,49,50,51]:
        idx = np.where(asegd == i)
        # temporary volume
        tmp = np.zeros( asegd.shape )
        tmp[idx] = 1
        tmp = imerode(tmp,se)
        idx = np.where(tmp == 1)
        gr_ncl[idx] = 1

    # without erosion
    for i in [13,17,18,26,52,53,54,58]:
        idx = np.where(asegd == i)
        gr_ncl[idx] = 1

    # remove remaining structure, e.g. brainstem
    remaining = np.zeros( asegd.shape )
    idx = np.where( asegd == 16 )
    remaining[idx] = 1

    # now remove all the structures from the white matter
    idx = np.where( (csfA != 0) | (csfB != 0) | (gr_ncl != 0) | (remaining != 0) )
    wmmask[idx] = 0
    log.info("Removing lateral ventricles and eroded grey nuclei and brainstem from white matter mask")

    # ADD voxels from 'cc_unknown.nii.gz' dataset
    ccun = ni.load(op.join(fs_dir, 'label', 'cc_unknown.nii.gz'))
    ccund = ccun.get_data()
    idx = np.where(ccund != 0)
    log.info("Add corpus callosum and unknown to wm mask")
    wmmask[idx] = 1
    # XXX add unknown dilation for connecting corpus callosum?
#    se2R = zeros(15,3,3); se2R(8:end,2,2)=1;
#    se2L = zeros(15,3,3); se2L(1:8,2,2)=1;
#    temp = (cc_unknown.img==1 | cc_unknown.img==2);
#    fsmask.img(imdilate(temp,se2R))    =  1;
#    fsmask.img(imdilate(temp,se2L))    =  1;
#    fsmask.img(cc_unknown.img==3)    =  1;
#    fsmask.img(cc_unknown.img==4)    =  1;

    # XXX: subtracting wmmask from ROI. necessary?
    for parkey, parval in gconf.parcellation.items():

        # check if we should subtract the cortical rois from this parcellation
        if parval.has_key('subtract_from_wm_mask'):
            if not bool(int(parval['subtract_from_wm_mask'])):
                continue
        else:
            continue

        log.info("Loading %s to subtract cortical ROIs from white matter mask" % ('ROI_%s.nii.gz' % parkey) )
        roi = ni.load(op.join(gconf.get_fs(), 'label', 'ROI_%s.nii.gz' % parkey))
        roid = roi.get_data()

        assert roid.shape[0] == wmmask.shape[0]

        pg = nx.read_graphml(parval['node_information_graphml'])

        for brk, brv in pg.nodes_iter(data=True):

            if brv['dn_region'] == 'cortical':

                log.info("Subtracting region %s with intensity value %s" % (brv['dn_region'], brv['dn_correspondence_id']))

                idx = np.where(roid == int(brv['dn_correspondence_id']))
                wmmask[idx] = 0

    # output white matter mask. crop and move it afterwards
    wm_out = op.join(fs_dir, 'mri', 'fsmask_1mm.nii.gz')
    img = ni.Nifti1Image(wmmask, fsmask.get_affine(), fsmask.get_header() )
    log.info("Save white matter mask: %s" % wm_out)
    ni.save(img, wm_out)

    
def crop_and_move_datasets():

    fs_dir = gconf.get_fs()
    fs_cmd_dir = gconf.get_cmp_fsout()
    reg_path = gconf.get_cmp_tracto_mask()

    log.info("Cropping and moving datasets to %s" % reg_path)

    # datasets to crop and move: (from, to)
    ds = [
          (op.join(fs_dir, 'mri', 'aseg.nii.gz'), op.join(reg_path, 'aseg.nii.gz') ),
          (op.join(fs_dir, 'mri', 'ribbon.nii.gz'), op.join(reg_path, 'ribbon.nii.gz') ),
          (op.join(fs_dir, 'mri', 'fsmask_1mm.nii.gz'), op.join(reg_path, 'fsmask_1mm.nii.gz') ),
          (op.join(fs_dir, 'label', 'cc_unknown.nii.gz'), op.join(reg_path, 'cc_unknown.nii.gz') )
          ]

    for p in gconf.parcellation.keys():
        ds.append( (op.join(fs_dir, 'label', 'ROI_%s.nii.gz' % p), op.join(reg_path, p, 'ROI_HR_th.nii.gz')) )
        ds.append( (op.join(fs_dir, 'label', 'ROIv_%s.nii.gz' % p), op.join(reg_path, p, 'ROIv_HR_th.nii.gz')) )

    orig = op.join(fs_dir, 'mri', 'orig', '001.mgz')

    for d in ds:
        log.info("Processing %s:" % d[0])

        # does it exist at all?
        if not op.exists(d[0]):
            raise Exception('File %s does not exist.' % d[0])
        # reslice to original volume because the roi creation with freesurfer
        # changed to 256x256x256 resolution
        mri_cmd = 'mri_convert -rl "%s" -rt nearest "%s" -nc "%s"' % (orig, d[0], d[1])
        runCmd( mri_cmd,log )


def generate_WM_and_GM_mask():

    log.info("Create the WM and GM mask")
    fs_dir = gconf.get_fs()
    reg_path = gconf.get_cmp_tracto_mask()

    # need to convert
    mri_cmd = 'mri_convert -i "%s/mri/aparc+aseg.mgz" -o "%s/mri/aparc+aseg.nii.gz"' % (fs_dir, fs_dir)
    runCmd( mri_cmd, log )
    mri_cmd = 'mri_convert -i "%s/mri/aparc.a2009s+aseg.mgz" -o "%s/mri/aparc.a2009s+aseg.nii.gz"' % (fs_dir, fs_dir)
    runCmd( mri_cmd, log )

    fout = op.join(fs_dir, 'mri', 'aparc+aseg.nii.gz')    
    niiAPARCimg = ni.load(fout)
    niiAPARCdata = niiAPARCimg.get_data()
    
    # override NativeFreesurfer
    fout = op.join(fs_dir, 'mri', 'aparc.a2009s+aseg.nii.gz')    
    niiAPARCimg = ni.load(fout)
    niiAPARCdata = niiAPARCimg.get_data()

    # mri_convert aparc+aseg.mgz aparc+aseg.nii.gz
    WMout = op.join(fs_dir, 'mri', 'fsmask_1mm.nii.gz')

    #%% label mapping
    # Using FreesurferColorLUT.txt
    # mappings are stored in mappings.ods

#    CORTICAL = {1 : [ 1, 2, 3, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34],
#                2 : [31,13, 9,21,27,25,19,29,15,23, 1,24, 4,30,26,11, 6, 2, 5,22,16,14,10,20,12, 7, 8,18,30,17, 3,28,33]}
#
#
#    SUBCORTICAL = {1:[48,49,50,51,52,53,54,58,59,60, 9,10,11,12,13,17,18,26,27,28],
#                   2:[34,34,35,36,37,40,41,38,39,39,75,75,76,77,78,81,82,79,80,80]}
#
#    OTHER = {1:[16],
#             2:[83]}
    
    # Mapping of NativeFreesurfer
    MAPPING = [[1,2012],[2,2019],[3,2032],[4,2014],[5,2020],[6,2018],[7,2027],[8,2028],[9,2003],[10,2024],[11,2017],[12,2026],
               [13,2002],[14,2023],[15,2010],[16,2022],[17,2031],[18,2029],[19,2008],[20,2025],[21,2005],[22,2021],[23,2011],
               [24,2013],[25,2007],[26,2016],[27,2006],[28,2033],[29,2009],[30,2015],[31,2001],[32,2030],[33,2034],[34,2035],
               [35,49],[36,50],[37,51],[38,52],[39,58],[40,53],[41,54],[42,1012],[43,1019],[44,1032],[45,1014],[46,1020],[47,1018],
               [48,1027],[49,1028],[50,1003],[51,1024],[52,1017],[53,1026],[54,1002],[55,1023],[56,1010],[57,1022],[58,1031],
               [59,1029],[60,1008],[61,1025],[62,1005],[63,1021],[64,1011],[65,1013],[66,1007],[67,1016],[68,1006],[69,1033],
               [70,1009],[71,1015],[72,1001],[73,1030],[74,1034],[75,1035],[76,10],[77,11],[78,12],[79,13],[80,26],[81,17],
               [82,18],[83,16]]

    # Mapping of Destrieux
    MAPPING = [[1,11106],[2,11107],[3,11108],[4,11101],[5,11102],[6,11103],[7,11104],[8,11105],[9,11109],[10,11110],[11,11111],[12,11112],[13,11113],[14,11114],[15,11115],[16,11116],[17,11117],[18,11118],[19,11119],[20,11120],[21,11121],[22,11122],[23,11123],[24,11124],[25,11127],[26,11125],[27,11126],[28,11128],[29,11129],[30,11130],[31,11131],[32,11132],[33,11137],[34,11138],[35,11133],[36,11134],[37,11135],[38,11136],[39,11139],[40,11140],[41,11141],[42,11143],[43,11144],[44,11145],[45,11146],[46,11147],[47,11148],[48,11149],[49,11150],[50,11151],[51,11152],[52,11153],[53,11154],[54,11155],[55,11156],[56,11157],[57,11160],[58,11158],[59,11159],[60,11161],[61,11162],[62,11165],[63,11163],[64,11164],[65,11166],[66,11167],[67,11168],[68,11169],[69,11170],[70,11171],[71,11172],[72,11173],[73,11174],[74,11175],[75,],[76,],[77,],[78,],[79,],[80,],[81,],[82,12106],[83,12107],[84,12108],[85,12101],[86,12102],[87,12103],[88,12104],[89,12105],[90,12109],[91,12110],[92,12111],[93,12112],[94,12113],[95,12114],[96,12115],[97,12116],[98,12117],[99,12118],[100,12119],[101,12120],[102,12121],[103,12122],[104,12123],[105,12124],[106,12127],[107,12125],[108,12126],[109,12128],[110,12129],[111,12130],[112,12131],[113,12132],[114,12137],[115,12138],[116,12133],[117,12134],[118,12135],[119,12136],[120,12139],[121,12140],[122,12141],[123,12143],[124,12144],[125,12145],[126,12146],[127,12147],[128,12148],[129,12149],[130,12150],[131,12151],[132,12152],[133,12153],[134,12154],[135,12155],[136,12156],[137,12157],[138,12160],[139,12158],[140,12159],[141,12161],[142,12162],[143,12165],[144,12163],[145,12164],[146,12166],[147,12167],[148,12168],[149,12169],[150,12170],[151,12171],[152,12172],[153,12173],[154,12174],[155,12175],[156,],[157,],[158,],[159,],[160,],[161,],[162,],[163,16]]

    WM = [2, 29, 32, 41, 61, 64, 59, 60, 27, 28] +  range(77,86+1) + range(100, 117+1) + range(155,158+1) + range(195,196+1) + range(199,200+1) + range(203,204+1) + [212, 219, 223] + range(250,255+1)
    # add
    # 59  Right-Substancia-Nigra
    # 60  Right-VentralDC
    # 27  Left-Substancia-Nigra
    # 28  Left-VentralDC

    log.info("WM mask....")
    #%% create WM mask    
    niiWM = np.zeros( niiAPARCdata.shape, dtype = np.uint8 )

    for i in WM:
         niiWM[niiAPARCdata == i] = 1

    # we do not add subcortical regions
#    for i in SUBCORTICAL[1]:
#         niiWM[niiAPARCdata == i] = 1

    img = ni.Nifti1Image(niiWM, niiAPARCimg.get_affine(), niiAPARCimg.get_header())
    log.info("Save to: " + WMout)
    ni.save(img, WMout)

    log.info("GM mask....")    
    #%% create GM mask (CORTICAL+SUBCORTICAL)
    #%  -------------------------------------
    for park in gconf.parcellation.keys():
        log.info("Parcellation: " + park)
        GMout = op.join(fs_dir, 'mri', 'ROIv_%s.nii.gz' % park)

        niiGM = np.zeros( niiAPARCdata.shape, dtype = np.uint8 )

        for ma in MAPPING:
            niiGM[ niiAPARCdata == ma[1]] = ma[0]

#        # % 33 cortical regions (stored in the order of "parcel33")
#        for idx,i in enumerate(CORTICAL[1]):
#            niiGM[ niiAPARCdata == (2000+i)] = CORTICAL[2][idx] # RIGHT
#            niiGM[ niiAPARCdata == (1000+i)] = CORTICAL[2][idx] + 41 # LEFT
#
#        #% subcortical nuclei
#        for idx,i in enumerate(SUBCORTICAL[1]):
#            niiGM[ niiAPARCdata == i ] = SUBCORTICAL[2][idx]
#
#        # % other region to account for in the GM
#        for idx, i in enumerate(OTHER[1]):
#            niiGM[ niiAPARCdata == i ] = OTHER[2][idx]

        log.info("Save to: " + GMout)        
        img = ni.Nifti1Image(niiGM, niiAPARCimg.get_affine(), niiAPARCimg.get_header())
        ni.save(img, GMout)

    log.info("[DONE]")


def crop_and_move_WM_and_GM():

    fs_dir = gconf.get_fs()
    reg_path = gconf.get_cmp_tracto_mask()

    log.info("Cropping and moving datasets to %s" % reg_path)

    # datasets to crop and move: (from, to)
    ds = [
          (op.join(fs_dir, 'mri', 'fsmask_1mm.nii.gz'), op.join(reg_path, 'fsmask_1mm.nii.gz') ),
          ]

    for p in gconf.parcellation.keys():
        ds.append( (op.join(fs_dir, 'mri', 'ROIv_%s.nii.gz' % p), op.join(reg_path, p, 'ROIv_HR_th.nii.gz')) )

    orig = op.join(fs_dir, 'mri', 'orig', '001.mgz')

    for d in ds:
        log.info("Processing %s:" % d[0])

        # does it exist at all?
        if not op.exists(d[0]):
            raise Exception('File %s does not exist.' % d[0])
        # reslice to original volume because the roi creation with freesurfer
        # changed to 256x256x256 resolution
        mri_cmd = 'mri_convert -rl "%s" -rt nearest "%s" -nc "%s"' % (orig, d[0], d[1])
        runCmd( mri_cmd,log )


def inspect(gconf):
    """ Inspect the results of this stage """
    log = gconf.get_logger()
    log.info("Check parcellation")

    if gconf.parcellation_scheme == "NativeFreesurfer":
        reg_path = gconf.get_cmp_tracto_mask()

        cmdstr = op.join(reg_path, 'fsmask_1mm.nii.gz ')

        for p in gconf.parcellation.keys():
            cmdstr += op.join(reg_path, p, 'ROIv_HR_th.nii.gz:colormap=lut:lut=%s ' % gconf.get_freeview_lut('NativeFreesurfer')['freesurferaparc'])

        freeview_view_cmd = 'freeview %s' % cmdstr
        runCmd( freeview_view_cmd, log )

    elif gconf.parcellation_scheme == "Lausanne2008":
        reg_path = gconf.get_cmp_tracto_mask()

        cmdstr = op.join(reg_path, 'fsmask_1mm.nii.gz ')

        for p in gconf.parcellation.keys():
            cmdstr += " " + op.join(reg_path, p, 'ROIv_HR_th.nii.gz')

        freeview_view_cmd = 'freeview %s' % cmdstr
        runCmd( freeview_view_cmd, log )


def extract(Z, shape, position, fill):
    """ Extract voxel neighbourhood
    Parameters
    ----------
    Z: the original data
    shape: tuple containing neighbourhood dimensions
    position: tuple containing central point indexes
    fill: value for the padding of Z
    Returns
    -------
    R: the neighbourhood of the specified point in Z
    """	
    R = np.ones(shape, dtype=Z.dtype) * fill # initialize output block to the fill value
    P = np.array(list(position)).astype(int) # position coordinates(numpy array)
    Rs = np.array(list(R.shape)).astype(int) # output block dimensions (numpy array)
    Zs = np.array(list(Z.shape)).astype(int) # original volume dimensions (numpy array)

    R_start = np.zeros(len(shape)).astype(int)
    R_stop = np.array(list(shape)).astype(int)
    Z_start = (P - Rs // 2)
    Z_start_cor = (np.maximum(Z_start,0)).tolist() # handle borders
    R_start = R_start + (Z_start_cor - Z_start)
    Z_stop = (P + Rs // 2) + Rs % 2
    Z_stop_cor = (np.minimum(Z_stop,Zs)).tolist() # handle borders
    R_stop = R_stop - (Z_stop - Z_stop_cor)

    R[R_start[0]:R_stop[0], R_start[1]:R_stop[1], R_start[2]:R_stop[2]] = Z[Z_start_cor[0]:Z_stop_cor[0], Z_start_cor[1]:Z_stop_cor[1], Z_start_cor[2]:Z_stop_cor[2]]

    return R


def run(conf):
    """ Run the first mask creation step

    Parameters
    ----------
    conf : PipelineConfiguration object
    subject_tuple : tuple, (subject_id, timepoint)
        Process the given subject

    """
    # setting the global configuration variable
    globals()['gconf'] = conf
    globals()['log'] = gconf.get_logger() 
    start = time()

    log.info("ROIv_HR_th.nii.gz / fsmask_1mm.nii.gz CREATION")
    log.info("=============================================")

    from os import environ
    env = environ
    env['SUBJECTS_DIR'] = gconf.get_subj_dir()

    if gconf.parcellation_scheme == "Lausanne2008":
        create_annot_label()
        create_roi()
        create_wm_mask()
        crop_and_move_datasets()
    elif gconf.parcellation_scheme == "NativeFreesurfer":
        generate_WM_and_GM_mask()
        crop_and_move_WM_and_GM()
    elif gconf.parcellation_scheme == "Destrieux":
        generate_WM_and_GM_mask()
        crop_and_move_WM_and_GM()

    log.info("Module took %s seconds to process." % (time()-start))

    if not len(gconf.emailnotify) == 0:
        msg = ["Mask creation", int(time()-start)]
        send_email_notification(msg, gconf, log)


def declare_inputs(conf):
    """Declare the inputs to the stage to the PipelineStatus object"""

    stage = conf.pipeline_status.GetStage(__name__)    
    fs_dir = conf.get_fs()    

    hemispheres = [ 'rh', 'lh' ]
    for hemi in hemispheres:
        conf.pipeline_status.AddStageInput(stage, op.join(fs_dir, 'surf'), '%s.sphere.reg' % (hemi), '%s.sphere.reg' % (hemi))

    conf.pipeline_status.AddStageInput(stage, op.join(fs_dir, 'mri'), 'aseg.mgz', 'aseg-mgz')


def declare_outputs(conf):
    """Declare the outputs to the stage to the PipelineStatus object"""

    stage = conf.pipeline_status.GetStage(__name__)
    reg_path = conf.get_cmp_tracto_mask()

    if conf.parcellation_scheme == "Lausanne2008":
        conf.pipeline_status.AddStageOutput(stage, reg_path, 'aseg.nii.gz', 'aseg-nii-gz')
        conf.pipeline_status.AddStageOutput(stage, reg_path, 'ribbon.nii.gz', 'ribbon-nii-gz')    
        conf.pipeline_status.AddStageOutput(stage, reg_path, 'fsmask_1mm.nii.gz', 'fsmask_1mm-nii-gz')
        conf.pipeline_status.AddStageOutput(stage, reg_path, 'cc_unknown.nii.gz', 'cc_unknown-nii-gz')

        for p in conf.parcellation.keys():
            conf.pipeline_status.AddStageOutput(stage, op.join(reg_path, p), 'ROI_HR_th.nii.gz', 'ROI_HR_th_%s-nii-gz' % (p))
            conf.pipeline_status.AddStageOutput(stage, op.join(reg_path, p), 'ROIv_HR_th.nii.gz', 'ROIv_HR_th_%s-nii-gz' % (p))

    elif conf.parcellation_scheme == "NativeFreesurfer" or conf.parcellation_scheme == "Destrieux":
        conf.pipeline_status.AddStageOutput(stage, reg_path, 'fsmask_1mm.nii.gz', 'fsmask_1mm-nii-gz')
        for p in conf.parcellation.keys():
            conf.pipeline_status.AddStageOutput(stage, op.join(reg_path, p), 'ROIv_HR_th.nii.gz', 'ROIv_HR_th_%s-nii-gz' % (p))


