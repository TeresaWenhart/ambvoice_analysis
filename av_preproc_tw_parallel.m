clear all
close all
clc


% initialize parallel pool
%--------------------------------------------------------------------------
nwork = 2;
try delete(gcp('nocreate')); end
parpool(nwork)


% folders
%--------------------------------------------------------------------------
STUDYdir  = ['/home/d.uzh.ch/twenha/data/MRI/1_nifti_prepared'];
FSdir     = ['/home/d.uzh.ch/twenha/data/MRI/2_freesurfer'];
ROIdir1   = [STUDYdir,'/rois_dti_vox'];
ROIdir2   = [STUDYdir,'/rois_dti_masks'];
FACTdir   = [STUDYdir,'/group_analysis/flexible_factorial_1'];
T1dir     = [STUDYdir,'/mean_t1'];
VBMdir    = [STUDYdir,'/vbm'];
NORMdir   = [VBMdir,'/mri'];

SPMdir    = ['/home/d.uzh.ch/twenha/Documents/MATLAB/spm12']; % path to your SPM12 folder
TPMdir    = [SPMdir,'/tpm'];
TEMPdir   = [SPMdir,'/toolbox/FieldMap/'];
CATdir    = [SPMdir,'/toolbox/cat12']; % you need to add CAT12 to your SPM toolboxes

mkdir(STUDYdir)
mkdir(FACTdir)
mkdir(T1dir)
mkdir(VBMdir)
mkdir(FSdir)

addpath(SPMdir)
addpath(['C/home/d.uzh.ch/twenha/data/MRI/scripts']) % path to your script folder
spm fmri % start SPM12 /or/ load defaults


% general settings
%--------------------------------------------------------------------------
bound = [-78 -112 -70; 78 76 85];


% runs
%--------------------------------------------------------------------------
run = {...
    'baseline';...
    'ins';...
    'gen';...
    'tim';...
    'voc';...
    'tva';...
    };


% subjects
%--------------------------------------------------------------------------
subj = {...
%     'av06a';...
%     'av09b';...
%     'av11a';...
%     'av13b';...
%     'av14b';...
%     'av15a';...
%     'av16a';...
%     'av17b';...
%     'av18b';...
%     'av19a';...
%     'av20a';...
%     'av21b';...
%     'av22a';...
%     'av25a';...
%     'av26b';...
%     'av27a';...
%     'av28a';...
%     'av30b';...
%     'av31a';...'
%     'av32a';...
%     'av33b';...%error in normalization; worked in a singel subject run.
    %'av34b';... normalization failed (could not read header of
    %file),worked in signe subject run
    %'av37b';...
    %'av38b';...
    %'av39b';...
    %'av42a';...
    %'av45b';...no wsanlm...ti.nii file -->normalisation failed?
    'av46a';...no wsanlm...ti.nii file -->normalisation failed?
    };





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % convert PAR/REC to NII
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% folders = {...
%     'tvaloc';...
%     'fieldmap';...
%     't1';...
%     'run01';...
%     'run02';...
%     'run03';...
%     'run04';...
%     'run05';...
%     'run06';...
%     };
% 
% for s = 1:size(subj,1)
%     for r = 1:size(folders,1)
%         mkdir([STUDYdir,'/',subj{s,1}]);
%         subjdir = [STUDYdir,'/',subj{s,1}];
%         mkdir([subjdir,'/',folders{1}]);
%         rundir  = [subjdir,'/',folders{r}];
%         cd(rundir); disp(subj{s,1}); disp(folders{r})
%         av_convert_rec(pwd)
%     end
% end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AC-PC position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!!
% !!! ... go to t1-folder ... manually re-orient T1 image to AC-PC position
% !!!





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-process T1 images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !!! ALWAYS RUN THIS FIRST LOOP --> defines the folders
% define: T1 reorient + filter
%-------------------------------------------------------------------------
clear SUBJdir RUNdir ANATdir ALLdir DARTELdir filt1 filt2 fname* path*
for s = 1:size(subj,1)
    for r = 1:size(run,1)
        mkdir([STUDYdir,'/',subj{s,1}]);
        SUBJdir{s,r} = [STUDYdir,'/',subj{s,1}];
        mkdir([SUBJdir{s,r},'/',run{1}]);
        mkdir([STUDYdir,'/dartel/',run{r}]);
        mkdir([FSdir,'/',subj{s,1},'/t1'])
        
        DARTELdir{s,r} = [STUDYdir,'/dartel/',run{r,1}];
        RUNdir{s,r}    = [SUBJdir{s,r},'/',run{r,1}];
        ANATdir{s,r}   = [SUBJdir{s,r},'/t1'];
        ALLdir{s,r}    = [SUBJdir{s,r},'/',run{r,1}];
        fname1{s,r}    = ['sanlm_s_',subj{s,1},'_t1.nii'];
        fname2{s,r}    = ['sanlm_rs_',subj{s,1},'_t1.nii'];
        path1{s,r}     = [RUNdir{s,r},'/sanlm_s_',subj{s,1},'_t1.nii'];
        path2{s,r}     = [FSdir,'/',subj{s,1},'/t1/sanlm_s_',subj{s,1},'_t1.nii'];
        path3{s,r}     = [RUNdir{s,r},'/sanlm_rs_',subj{s,1},'_t1.nii'];
        path4{s,r}     = [FSdir,'/',subj{s,1},'/t1/',subj{s,1},'.t1.nii'];
    end
    
    for r = 1
        cd(RUNdir{s,r})
        copyfile([ANATdir{s,r},'/s_',subj{s,1},'_t1.nii'],RUNdir{s,r})
        ima{s,r} = ['s_',subj{s,1},'_t1.nii'];
        
        matlabbatch = [];
        matlabbatch{1}.spm.tools.cat.tools.sanlm.data = cellstr(['s_',subj{s,1},'_t1.nii']);
        save filter1 matlabbatch
        filt1{s,r} = [pwd,'/filter1.mat'];
        
        matlabbatch = [];
        matlabbatch{1}.spm.tools.cat.tools.sanlm.data = cellstr(['rs_',subj{s,1},'_t1.nii']);
        save filter2 matlabbatch
        filt2{s,r} = [pwd,'/filter2.mat'];
    end
end

% run: T1 reorient + filter
parfor s = 1:size(subj,1)%change to "parfor" for parallel processing
    for r = 1
        cd(RUNdir{s,r})
        reorient_mri_image(ima{s,r})
        
        temp = [];
        temp = load(filt1{s,r});
        spm_jobman('run',temp.matlabbatch);
        
        temp = [];
        temp = load(filt2{s,r});
        spm_jobman('run',temp.matlabbatch);
        
        for i = 2:6
            copyfile(fname1{s,i},ALLdir{s,i})
            copyfile(fname2{s,i},ALLdir{s,i})
        end
        
        copyfile(path1{s,r},path2{s,r})
        copyfile(path3{s,r},path4{s,r})
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. step pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define: realign
%--------------------------------------------------------------------------
clear realign
for s = 1:size(subj,1)
    for r = 1:size(run,1)
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        matlabbatch = [];
        ima  = [];
        ima  = cellstr(spm_select('FPlistRec',pwd,'^f.*.nii$'));
        temp = cellstr(repmat(',1',length(ima),1));
        ima  = strcat(ima,temp);
        matlabbatch{1}.spm.spatial.realign.estwrite.data{1,1}        = ima;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep     = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm    = 5;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm     = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp  = 7;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight  = '';
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which   = [2 1];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp  = 7;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask    = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix  = 'r';
        save realign matlabbatch
        realign{s,r} = [pwd,'/realign.mat'];
    end
end

% run: realign
parfor s = 1:size(subj,1) %change to "parfor" for parallel processing
    for r = 1:size(run,1)
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        temp = [];
        temp = load(realign{s,r});
        spm_jobman('run',temp.matlabbatch);
    end
end


% define: slicetime
%--------------------------------------------------------------------------
clear slicetime
for s = 1:size(subj,1)
    for r = 1:size(run,1)
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        matlabbatch = [];
        matlabbatch{1}.spm.temporal.st.scans{1} = cellstr(spm_select('FPlistRec',pwd,'^rf.*.nii$'));
        matlabbatch{1}.spm.temporal.st.nslices  = 28;
        matlabbatch{1}.spm.temporal.st.tr       = 1.60;
        matlabbatch{1}.spm.temporal.st.ta       = 1.60-(1.60/28);
        matlabbatch{1}.spm.temporal.st.so       = [1:1:28];
        matlabbatch{1}.spm.temporal.st.refslice = 14;
        matlabbatch{1}.spm.temporal.st.prefix   = 'a';
        save slicetime matlabbatch
        slicetime{s,r} = [pwd,'/slicetime.mat'];
    end
end

% run: slicetime
parfor s = 1:size(subj,1)%change to "parfor" for parallel processing
    for r = 1:size(run,1)
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        temp = [];
        temp = load(slicetime{s,r});
        spm_jobman('run',temp.matlabbatch);
    end
end


% define: coregister
%--------------------------------------------------------------------------
clear coregister
for s = 1:size(subj,1)
    for r = 1:size(run,1)
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        matlabbatch = []; ima = []; ima1 = []; ima2 = [];
        ima1 = cellstr(spm_select('FPlistRec',pwd,'^arf.*.nii$'));
        ima2 = cellstr(spm_select('FPlistRec',pwd,'^rf.*.nii$'));
        ima = [ima1;ima2];
        matlabbatch{1}.spm.spatial.coreg.estimate.ref               = cellstr(['sanlm_rs_',subj{s,1},'_t1.nii']);
        matlabbatch{1}.spm.spatial.coreg.estimate.source            = cellstr(spm_select('FPlistRec',pwd,'^meanf.*.nii$'));
        matlabbatch{1}.spm.spatial.coreg.estimate.other             = ima;
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = char('nmi');
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol      = [...
            0.02  0.02  0.02 ...
            0.001 0.001 0.001 ...
            0.01  0.01  0.01 ...
            0.001 0.001 0.001];
        save coregister matlabbatch
        coregister{s,r} = [pwd,'/coregister.mat'];
    end
end

% run: coregister
parfor s = 1:size(subj,1)%change to "parfor" for parallel processing
    for r = 1:size(run,1)
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        temp = [];
        temp = load(coregister{s,r});
        spm_jobman('run',temp.matlabbatch);
    end
end





%% SEGMENTATION 

% define: segment --> NewSegmentation
%--------------------------------------------------------------------------
clear T1segment1 T1segment2 path*
for s = 1:size(subj,1)
    for r = 1
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        % delete c1sanlm* c2sanlm* c3sanlm* c4sanlm* c5sanlm* c6sanlm*
        matlabbatch = [];
        ngaus = [1 1 2 3 4 2];
        for tissue = 1:6
            matlabbatch{1}.spm.spatial.preproc.tissue(tissue).tpm    = {[TPMdir,'/TPM.nii,',num2str(tissue)]};
            matlabbatch{1}.spm.spatial.preproc.tissue(tissue).ngaus  = ngaus(tissue);
            matlabbatch{1}.spm.spatial.preproc.tissue(tissue).native = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(tissue).warped = [0 0];
        end
        matlabbatch{1}.spm.spatial.preproc.channel.vols     = cellstr(spm_select('FPlistRec',pwd,'^sanlm_s_.*_t1.nii$'));
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg  = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write    = [0 0];
        matlabbatch{1}.spm.spatial.preproc.warp.mrf         = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup     = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.reg         = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.preproc.warp.affreg      = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm        = 0;
        matlabbatch{1}.spm.spatial.preproc.warp.samp        = 3;
        matlabbatch{1}.spm.spatial.preproc.warp.write       = [1 1];
        save T1segment1 matlabbatch
        T1segment1{s,r} = [pwd,'/T1segment1.mat'];
        
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        matlabbatch = [];
        ngaus = [1 1 2 3 4 2];
        for tissue = 1:6
            matlabbatch{1}.spm.spatial.preproc.tissue(tissue).tpm    = {[TPMdir,'/TPM.nii,',num2str(tissue)]};
            matlabbatch{1}.spm.spatial.preproc.tissue(tissue).ngaus  = ngaus(tissue);
            matlabbatch{1}.spm.spatial.preproc.tissue(tissue).native = [1 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(tissue).warped = [0 0];
        end
        matlabbatch{1}.spm.spatial.preproc.channel.vols     = cellstr(spm_select('FPlistRec',pwd,'^sanlm_rs_.*_t1.nii$'));
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg  = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write    = [0 0];
        matlabbatch{1}.spm.spatial.preproc.warp.mrf         = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup     = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.reg         = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.preproc.warp.affreg      = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm        = 0;
        matlabbatch{1}.spm.spatial.preproc.warp.samp        = 3;
        matlabbatch{1}.spm.spatial.preproc.warp.write       = [1 1];
        save T1segment2 matlabbatch
        T1segment2{s,r} = [pwd,'/T1segment2.mat'];
        
        for tissue = 1:6
            path1{s,tissue} = [RUNdir{s,r},'/rc',num2str(tissue),'sanlm_s_',subj{s,1},'_t1.nii'];
            path2{s,tissue} = [DARTELdir{s,r},'/rc',num2str(tissue),'sanlm_s_',subj{s,1},'_t1.nii'];
        end
        for tissue = 1:6
            path3{s,tissue} = [RUNdir{s,r},'/rc',num2str(tissue),'sanlm_rs_',subj{s,1},'_t1.nii'];
            path4{s,tissue} = [DARTELdir{s,r},'/rc',num2str(tissue),'sanlm_rs_',subj{s,1},'_t1.nii'];
        end
    end
end

% run: segment
parfor s = 1:size(subj,1)%change to "parfor" for parallel processing
    for r = 1
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        if r == 1
            temp = [];
            temp = load(T1segment1{s,r});
            spm_jobman('run',temp.matlabbatch);
            
            temp = [];
            temp = load(T1segment2{s,r});
            spm_jobman('run',temp.matlabbatch);
            
            for tissue = 1:6
                copyfile(path1{s,tissue},path2{s,tissue});
            end
            
            for tissue = 1:6
                copyfile(path3{s,tissue},path4{s,tissue});
            end
        end
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CAT12: segmentation + VBM pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(VBMdir)
delete([VBMdir,'/report/catreport_sanlm_rs_',subj{end,1},'_t1.pdf'])
cat12('expert')

% copy files to VBM folder
for s = 1:size(subj,1)
    copyfile(...
        [STUDYdir,'/',subj{s,1},'/',run{1},'/sanlm_s_',char(subj{s,1}),'_t1.nii'],...
        [VBMdir,'/sanlm_s_',char(subj{s,1}),'_t1.nii']);
    copyfile(...
        [STUDYdir,'/',subj{s,1},'/',run{1},'/sanlm_rs_',char(subj{s,1}),'_t1.nii'],...
        [VBMdir,'/sanlm_rs_',char(subj{s,1}),'_t1.nii']);
end

% setup segmentation with CAT12
ima = {};
matlabbatch = [];
count = 0;
for s = 1:size(subj,1)
    % include only subjects for which y_*.nii does not exist
    if ~exist([NORMdir,'/mwp2sanlm_rs_',subj{s,1},'_t1'],'file')
        count = count + 1;
        ima{count,1} = [VBMdir,'/sanlm_rs_',subj{s,1},'_t1.nii'];
    end
end
matlabbatch{1}.spm.tools.cat.estwrite.data         = ima;
matlabbatch{1}.spm.tools.cat.estwrite.nproc        = 10;
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm     = {[SPMdir,'/tpm/TPM.nii']};
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg  = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.opts.accstr  = 0.5;

matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.APP              = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.NCstr            = -Inf;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASstr           = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.gcutstr          = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.cleanupstr       = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.WMHC             = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.SLC              = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.restypes.best    = [0.5 0.3];
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.dartel.darteltpm = {[CATdir,'/templates_volumes/Template_1_IXI555_MNI152.nii']};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox                           = 1.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtres                = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.scale_cortex          = 0.7;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.add_parahipp          = 0.1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.close_parahipp        = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.ignoreErrors            = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.verb                    = 2;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.print                   = 2;

matlabbatch{1}.spm.tools.cat.estwrite.output.surface                            = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40             = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra              = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers            = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr               = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.aal                = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.mori               = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.anatomy            = 1;

matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native      = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.warped      = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod         = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel      = 3;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native      = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.warped      = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod         = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel      = 3;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native     = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped     = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod        = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel     = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native     = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped     = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod        = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel     = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native      = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped      = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod         = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel      = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native   = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.dartel   = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.native   = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped   = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel   = 3;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.native    = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped    = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.dartel    = 3;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.native     = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped     = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel     = 3;
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps          = [1 1];
save vbm_cat12 matlabbatch
spm_jobman('run', matlabbatch);

% only continue when last segmentation file is written
while exist([VBMdir,'/report/catreport_sanlm_rs_',subj{end,1},'_t1.pdf'],'file') == 0
    pause(10)
end

%re-initialize parallel pool
try delete(gcp('nocreate')); end
parpool(nwork)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. step pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define: normalize
%--------------------------------------------------------------------------
clear fnormalise*
for s = 1:size(subj,1)
    for r = 1:size(run,1)
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        matlabbatch = [];
        matlabbatch{1}.spm.spatial.normalise.write.subj.def        = cellstr([NORMdir,'/y_sanlm_rs_',subj{s,1},'_t1.nii']);
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample   = cellstr(spm_select('FPlistRec',pwd,'^rf.*.nii$'));
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb     = bound;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox    = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 7;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w2';
        save fnormalise1 matlabbatch
        fnormalise1{s,r} = [pwd,'/fnormalise1.mat'];
    end
end

% run: normalize
parfor s = 1:size(subj,1)%change to parfor for parallel processing
    for r = 1:size(run,1)
        cd(RUNdir{s,r})
        temp = [];
        temp = load(fnormalise1{s,r});
        spm_jobman('run',temp.matlabbatch);
    end
end


% define: smooth
%--------------------------------------------------------------------------
clear smooth*
for s = 1:size(subj,1)
    for r = 1:size(run,1)
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        matlabbatch = [];
        matlabbatch{1}.spm.spatial.smooth.data   = cellstr(spm_select('FPlistRec',pwd,'^w2rf.*.nii$'));
        matlabbatch{1}.spm.spatial.smooth.fwhm   = [8 8 8];
        matlabbatch{1}.spm.spatial.smooth.dtype  = 0;
        matlabbatch{1}.spm.spatial.smooth.im     = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's8';
        save smooth8 matlabbatch
        smooth8{s,r} = [pwd,'/smooth8.mat'];
        
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        matlabbatch = [];
        matlabbatch{1}.spm.spatial.smooth.data   = cellstr(spm_select('FPlistRec',pwd,'^w2rf.*.nii$'));
        matlabbatch{1}.spm.spatial.smooth.fwhm   = [6 6 6];
        matlabbatch{1}.spm.spatial.smooth.dtype  = 0;
        matlabbatch{1}.spm.spatial.smooth.im     = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's6';
        save smooth6 matlabbatch
        smooth6{s,r} = [pwd,'/smooth6.mat'];
    end
end

% run: smooth
parfor s = 1:size(subj,1)%change to parfor for parallel processing
    for r = 1:size(run,1)
        cd(RUNdir{s,r})
        temp = [];
        temp = load(smooth8{s,r});
        spm_jobman('run',temp.matlabbatch);
        
        temp = [];
        temp = load(smooth6{s,r});
        spm_jobman('run',temp.matlabbatch);
    end
end


% define: T1 normalize
%--------------------------------------------------------------------------
clear T1nornmalise* SEGnornmalise
for s = 1:size(subj,1)
    for r = 1
        % normalize T1
        %------------------------------------------------------------------
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        matlabbatch = [];
        matlabbatch{1}.spm.spatial.normalise.write.subj.def        = cellstr([NORMdir,'/y_sanlm_rs_',subj{s,1},'_t1.nii']);
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample   = cellstr(spm_select('FPlistRec',pwd,'^sanlm_s_.*_t1.nii$'));
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb     = bound;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox    = [1 1 1];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 7;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
        save T1nornmalise1 matlabbatch
        T1nornmalise1{s,r} = [pwd,'/T1nornmalise1.mat'];
        
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        matlabbatch = [];
        matlabbatch{1}.spm.spatial.normalise.write.subj.def        = cellstr([NORMdir,'/y_sanlm_rs_',subj{s,1},'_t1.nii']);
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample   = cellstr(spm_select('FPlistRec',pwd,'^sanlm_rs_.*_t1.nii$'));
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb     = bound;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox    = [1 1 1];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 7;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
        save T1nornmalise2 matlabbatch
        T1nornmalise2{s,r} = [pwd,'/T1nornmalise2.mat'];
        
        % normalise segmented tissues
        %--------------------------------------------------------------
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        copyfile([NORMdir,'/p1sanlm_rs_',subj{s,1},'_t1.nii'],RUNdir{s,r});
        copyfile([NORMdir,'/p2sanlm_rs_',subj{s,1},'_t1.nii'],RUNdir{s,r});
        matlabbatch = [];
        matlabbatch{1}.spm.spatial.normalise.write.subj.def        = cellstr([NORMdir,'/y_sanlm_rs_',subj{s,1},'_t1.nii']);
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample   = cellstr(spm_select('FPlistRec',pwd,'^p.*sanlm_rs_*'));
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb     = bound;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox    = [1 1 1];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 7;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
        save SEGnornmalise matlabbatch
        SEGnornmalise{s,r} = [pwd,'/SEGnornmalise.mat'];
    end
end


% run: T1 normalize
parfor s = 1:size(subj,1)%change to parfor for parallel processing
    for r = 1
        cd(RUNdir{s,r})
        temp = [];
        temp = load(T1nornmalise1{s,r});
        spm_jobman('run',temp.matlabbatch);
        
        temp = [];
        temp = load(T1nornmalise2{s,r});
        spm_jobman('run',temp.matlabbatch);
        
        temp = [];
        temp = load(SEGnornmalise{s,r});
        spm_jobman('run',temp.matlabbatch);
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. step pre-processing (no parallel loop necessary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop subjects: create brain mask from white + grey matter segmentation
%--------------------------------------------------------------------------
for s = 1:size(subj,1)
    for r = 1
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        
        % single masks
        for cs = 1:2
            matlabbatch = [];
            tissues = cellstr(spm_select('FPlistRec',pwd,'^wp.*sanlm_s*'));
            matlabbatch{1}.spm.util.imcalc.input(1,1)     = tissues(cs);
            matlabbatch{1}.spm.util.imcalc.output         = [RUNdir{s,r},'/c',subj{s,1},'_',num2str(cs),'.nii'];
            matlabbatch{1}.spm.util.imcalc.outdir         = {RUNdir{s,r}};
            matlabbatch{1}.spm.util.imcalc.expression     = 'i1>0.03';
            matlabbatch{1}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = -7;
            matlabbatch{1}.spm.util.imcalc.options.dtype  = 4;
            spm_jobman('run',matlabbatch);
        end
        
        % grey + white mask
        matlabbatch = [];
        matlabbatch{1}.spm.util.imcalc.input(1,1)     = cellstr([RUNdir{s,r},'/c',subj{s,1},'_1.nii']);
        matlabbatch{1}.spm.util.imcalc.input(2,1)     = cellstr([RUNdir{s,r},'/c',subj{s,1},'_2.nii']);
        matlabbatch{1}.spm.util.imcalc.output         = [RUNdir{s,r},'/c_',subj{s,1},'_greywhite.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir         = {RUNdir{s,r}};
        matlabbatch{1}.spm.util.imcalc.expression     = '(i1+i2)>0';
        matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 8;
        matlabbatch{1}.spm.util.imcalc.options.dtyp   = 4;
        spm_jobman('run',matlabbatch);
        
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        matlabbatch = [];
        matlabbatch{1}.spm.spatial.smooth.data   = cellstr(spm_select('FPlistRec',pwd,['^c_',subj{s,1},'_greywhite.nii$']));
        matlabbatch{1}.spm.spatial.smooth.fwhm   = [0.5 0.5 0.5];
        matlabbatch{1}.spm.spatial.smooth.dtype  = 0;
        matlabbatch{1}.spm.spatial.smooth.im     = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        save smooth matlabbatch
        spm_jobman('run',matlabbatch);
        
        cd(RUNdir{s,r}); disp(subj{s,1}); disp(run{r})
        matlabbatch = [];
        matlabbatch{1}.spm.util.imcalc.input(1,1)     = cellstr(spm_select('FPlistRec',pwd,['^sc_',subj{s,1},'_greywhite.nii$']));
        matlabbatch{1}.spm.util.imcalc.output         = [RUNdir{s,r},'/bsc_',subj{s,1},'_greywhite.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir         = {RUNdir{s,r}};
        matlabbatch{1}.spm.util.imcalc.expression     = '(i1)>0';
        matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 8;
        matlabbatch{1}.spm.util.imcalc.options.dtyp   = 4;
        spm_jobman('run',matlabbatch);
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. CAT12 (XCAT): mean anatomical images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(T1dir)
clear matlabbatch expression
for s = 1:size(subj,1)
    matlabbatch{1}.spm.util.imcalc.input(s,1) = ...
        cellstr([STUDYdir,'/',subj{s,1},'/',char(run{1}),'/wsanlm_rs_',subj{s,1},'_t1.nii']);
end
expression = 'i1';

if size(subj,1) == 1
    expression = expression;
else
    for s = 2:size(subj,1)
        expression = [expression,'+i',num2str(s)];
    end
end
expression = ['((',expression,')/',num2str(length(subj)),')>0.5'];
matlabbatch{1}.spm.util.imcalc.expression     = expression;
matlabbatch{1}.spm.util.imcalc.output         = [T1dir,'/xcat.av.mask.grey.nii'];
matlabbatch{1}.spm.util.imcalc.outdir         = {T1dir};
matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 8;
matlabbatch{1}.spm.util.imcalc.options.dtyp   = 4;
spm_jobman('run',matlabbatch);

copyfile('xcat.av.mask.grey.nii',...
    [FACTdir,'/xcat.av.mask.grey.nii']);
copyfile('xcat.av.mask.grey.nii',...
    [FACTdir,'/xcat.av.mask.grey.nii']);

% mean greywhite
cd(T1dir)
clear matlabbatch expression
for s = 1:size(subj,1)
    matlabbatch{1}.spm.util.imcalc.input(s,1) = ...
        cellstr([STUDYdir,'/',subj{s,1},'/',char(run{1}),'/c_',subj{s,1},'_greywhite.nii']);
end
expression = 'i1';
if size(subj,1) == 1
    expression = expression;
else
    for s = 2:size(subj,1)
        expression = [expression,'+i',num2str(s)];
    end
end
expression = ['(',expression,')/',num2str(length(subj))];
matlabbatch{1}.spm.util.imcalc.expression     = expression;
matlabbatch{1}.spm.util.imcalc.output         = [T1dir,'/xcat.av.mean.greywhite.nii'];
matlabbatch{1}.spm.util.imcalc.outdir         = {T1dir};
matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 8;
matlabbatch{1}.spm.util.imcalc.options.dtyp   = 4;
spm_jobman('run',matlabbatch);

copyfile('xcat.av.mean.greywhite.nii',...
    [FACTdir,'/xcat.av.mean.greywhite.nii']);
copyfile('xcat.av.mean.greywhite.nii',...
    [FACTdir,'/xcat.av.mean.greywhite.nii']);

% mask greywhite
cd(T1dir)
clear matlabbatch expression
for s = 1:size(subj,1)
    matlabbatch{1}.spm.util.imcalc.input(s,1) = ...
        cellstr([STUDYdir,'/',subj{s,1},'/',char(run{1}),'/c',subj{s,1},'_1.nii']);
end
expression = 'i1';
if size(subj,1) == 1
    expression = expression;
else
    for s = 2:size(subj,1)
        expression = [expression,'+i',num2str(s)];
    end
end
expression = ['((',expression,')/',num2str(length(subj)),')>0.5'];
matlabbatch{1}.spm.util.imcalc.expression     = expression;
matlabbatch{1}.spm.util.imcalc.output         = [T1dir,'/xcat.av.mask.grey.nii'];
matlabbatch{1}.spm.util.imcalc.outdir         = {T1dir};
matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 8;
matlabbatch{1}.spm.util.imcalc.options.dtyp   = 4;
spm_jobman('run',matlabbatch);

copyfile('xcat.av.mask.grey.nii',...
    [FACTdir,'/xcat.av.mask.grey.nii']);
copyfile('xcat.av.mask.grey.nii',...
    [FACTdir,'/xcat.av.mask.grey.nii']);


