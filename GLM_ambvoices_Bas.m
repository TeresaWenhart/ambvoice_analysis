% CREATE GLM
% 
%__________________________________________________________________________
% (C) Teresa Wenhart, 2020


clear; close all; clc


% Paths
%--------------------------------------------------------------------------
ExpPath    = '/home/d.uzh.ch/twenha';
DataPath   = fullfile(ExpPath,'data','MRI','1_nifti_prepared');
behavPath  = fullfile(ExpPath,'data', 'behavioral');
MaskPath   = '/home/d.uzh.ch/twenha/Documents/MATLAB/spm12/toolbox/FieldMap/brainmask.nii,1';
ResultPath = fullfile(ExpPath,'models','first_level');
[~,~]      = mkdir(ResultPath);


% subject list
%--------------------------------------------------------------------------
sList = {...
    'av06a';...
    'av09b';...
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
%     'av33b';...
%     'av34b';...
   % 'av37b';...unable to open Bas behavioral file
    %'av38b';...
%     'av39b';...
%     'av42a';...
%     'av45b';...
%     'av46a';...
    };


% MRI parameters
%--------------------------------------------------------------------------
% 
nslices = 28;
TR      = 1.6;
stimduration = 0.6; 


% launch SPM (add toolbox to path if necessary)
%--------------------------------------------------------------------------
if exist('spm')~=2, addToolboxes(1); end
spm('defaults','fmri');
spm_jobman('initcfg')


%--------------------------------------------------------------------------
% create GLM for Baseline
%--------------------------------------------------------------------------
% 
runs = {'Bas';'Voc';'Tim';'Gen';'Ins'};
epis = {'baseline', 'voc', 'tim', 'gen', 'ins'};



for sIDX = 1:(numel(sList))
    r=1
    clear matlabbatch TD

    %----------------------------------------------------------------------
    % paths
    %----------------------------------------------------------------------
    dPath      = fullfile(DataPath,sList{sIDX});
    T1Path     = fullfile(dPath,'t1');
    ModelPath  = fullfile(ResultPath,sList{sIDX},epis{r});
    %[~,~]      = rmdir(ResultPath,'s');
    [~,~]      = mkdir(ModelPath);
    
    EPIPathList = fullfile(dPath,epis);
    cd(dPath)
    %----------------------------------------------------------------------
    % fMRI model specification: Baseline
    %----------------------------------------------------------------------
    % 
    
    clear EPIflist Motionfn fmri_spec
    
    fmri_spec.dir            = cellstr(ModelPath);
    fmri_spec.timing.units   = 'secs';
    fmri_spec.timing.RT      = TR;
    fmri_spec.timing.fmri_t  = 16;
    fmri_spec.timing.fmri_t0 = 8;
    
    EPIflist{1,r} = cellstr(spm_select('FPList',EPIPathList{r},['^s8w2rf.*\.nii$']));
    Motionfn{1,r} = cellstr(spm_select('FPList',EPIPathList{r},['^rp.*\.txt$']));
    
            
    %----------------------------------------------------------------------
    %Specify sessions
    %----------------------------------------------------------------------
    %one session is one task & run (runs =
    %{'Bas';'Voc';'Tim';'Gen';'Ins'};)
    fmri_spec.sess(1).scans = EPIflist{1,r}(1:end);
    
    % collect event timings
    %----------------------------------------------------------------------
    
    clear behav
    behav.data.Bas = readtable(fullfile(behavPath,sList{sIDX}, ['AV_Bas_' sList{sIDX} '.csv']));
    behav.data.Voc = readtable(fullfile(behavPath,sList{sIDX}, ['AV_Voc_' sList{sIDX} '.csv']));
    behav.data.Tim = readtable(fullfile(behavPath,sList{sIDX}, ['AV_Tim_' sList{sIDX} '.csv']));
    behav.data.Gen = readtable(fullfile(behavPath,sList{sIDX}, ['AV_Gen_' sList{sIDX} '.csv']));
    behav.data.Ins = readtable(fullfile(behavPath,sList{sIDX}, ['AV_Ins_' sList{sIDX} '.csv']));
 
    Ons.raw{1} = behav.data.Bas(:,'onset');
    Ons.raw{2} = behav.data.Voc(:,'onset');
    Ons.raw{3} = behav.data.Tim(:,'onset');
    Ons.raw{4} = behav.data.Gen(:,'onset');
    Ons.raw{5} = behav.data.Ins(:,'onset');
    
    % get trial indices to assign ambiguity as parametric
    % modulator below
    Tr{1} = behav.data.Bas(:,'morphrate1');
    Tr{2} = behav.data.Voc(:,'morphrate1');
    Tr{3} = behav.data.Tim(:,'morphrate1');
    Tr{4} = behav.data.Gen(:,'morphrate1');
    Tr{5} = behav.data.Ins(:,'morphrate1');
    
    %sound conditions for each run
    %BASELINE
    condition_Bas = {...
            'fem';...
            'mal';...
            'cla';...
            'bcl';...
            'vlc';...
            'vla';...
            };
    %Voice
            condition_Voc = {...
            'ins_1';...
            'voc_1';...
            'ins_08';...
            'voc_08';...
            'ins_06';...
            'voc_06';...
            };
    %Timbre
            condition_Tim = {...
            'vlc_1';...
            'cla_1';...
            'vlc_08';...
            'cla_08';...
            'vlc_06';...
            'cla_06';...
            };    
    %Gender
    condition_Gen = {...
            'mal_1';...
            'fem_1';...
            'mal_08';...
            'fem_08';...
            'mal_06';...
            'fem_06';...
            };
    %Instrument
        condition_Ins = {...
            'ins_1';...
            'voc_1';...
            'ins_08';...
            'voc_08';...
            'ins_06';...
            'voc_06';...
            };
    condition{1}.sound = condition_Bas
    condition{2}.sound = condition_Voc
    condition{3}.sound = condition_Tim
    condition{4}.sound = condition_Gen
    condition{5}.sound = condition_Ins
    
    % get onset times for secific sound conditions
    % BASELINE
    idXmal = contains(behav.data.Bas.Category1, 'Mal')
    Ons.run{1}.sound{1} = Ons.raw{1}(idXmal,:)
    
    idXfem = contains(behav.data.Bas.Category1, 'Fem')
    Ons.run{1}.sound{2} = Ons.raw{1}(idXfem,:)
    
    idXcla = contains(behav.data.Bas.Category1, 'Cla') 
    Ons.run{1}.sound{3} = Ons.raw{1}(idXcla,:)
    
    idXbcl = contains(behav.data.Bas.Category1, 'Bcl') 
    Ons.run{1}.sound{4} = Ons.raw{1}(idXbcl,:)
    
    idXvlc = contains(behav.data.Bas.Category1, 'Vlc') 
    Ons.run{1}.sound{5} = Ons.raw{1}(idXvlc,:)
    
    idXvla = contains(behav.data.Bas.Category1, 'Vla') 
    Ons.run{1}.sound{6} = Ons.raw{1}(idXvla,:)
    
     % VOICE
    idXins_1 = [behav.data.Voc.morphrate1] == 1
    Ons.run{2}.sound{1} = Ons.raw{2}(idXins_1,:)
    
    idXvoc_1 = [behav.data.Voc.morphrate1] == 0
    Ons.run{2}.sound{2} = Ons.raw{2}(idXvoc_1,:)
    
    idXins_08 = [behav.data.Voc.morphrate1] == 0.8
    Ons.run{2}.sound{3} = Ons.raw{2}(idXins_08,:)
    
    idXvoc_08 = [behav.data.Voc.morphrate1] == 0.2
    Ons.run{2}.sound{4} = Ons.raw{2}(idXvoc_08,:)
    
    idXins_06 = [behav.data.Voc.morphrate1] == 0.6
    Ons.run{2}.sound{5} = Ons.raw{2}(idXins_06,:)
    
    idXvoc_06 = [behav.data.Voc.morphrate1] == 0.4
    Ons.run{2}.sound{6} = Ons.raw{2}(idXvoc_06,:)
    % TIMBRE
    
    idXvlc_1 = [behav.data.Tim.morphrate1] == 0
    Ons.run{3}.sound{1} = Ons.raw{3}(idXvlc_1,:)
    
    idXcla_1 = [behav.data.Tim.morphrate1] == 1
    Ons.run{3}.sound{2} = Ons.raw{3}(idXcla_1,:)
    
    idXvlc_08 = [behav.data.Tim.morphrate1] == 0.2
    Ons.run{3}.sound{3} = Ons.raw{3}(idXvlc_08,:)
    
    idXcla_08 = [behav.data.Tim.morphrate1] == 0.8
    Ons.run{3}.sound{4} = Ons.raw{3}(idXcla_08,:)
    
    idXvlc_06 = [behav.data.Tim.morphrate1] == 0.4
    Ons.run{3}.sound{5} = Ons.raw{3}(idXvlc_06,:)
    
    idXcla_06 = [behav.data.Tim.morphrate1] == 0.6
    Ons.run{3}.sound{6} = Ons.raw{3}(idXcla_06,:)
    
    
    % GENDER
    idXmal_1 = [behav.data.Gen.morphrate1] == 0
    Ons.run{4}.sound{1} = Ons.raw{4}(idXmal_1,:)
    
    idXfem_1 = [behav.data.Gen.morphrate1] == 1
    Ons.run{4}.sound{2} = Ons.raw{4}(idXfem_1,:)
    
    idXmal_08 = [behav.data.Gen.morphrate1] == 0.2
    Ons.run{4}.sound{3} = Ons.raw{4}(idXmal_08,:)
    
    idXfem_08 = [behav.data.Gen.morphrate1] == 0.8
    Ons.run{4}.sound{4} = Ons.raw{4}(idXfem_08,:)
    
    idXmal_06 = [behav.data.Gen.morphrate1] == 0.4
    Ons.run{4}.sound{5} = Ons.raw{4}(idXmal_06,:)
    
    idXfem_06 = [behav.data.Gen.morphrate1] == 0.6
    Ons.run{4}.sound{6} = Ons.raw{4}(idXfem_06,:)
    
    % INSTRUMENT
    
    idXins_1 = [behav.data.Ins.morphrate1] == 1
    Ons.run{5}.sound{1} = Ons.raw{5}(idXins_1,:)
    
    idXvoc_1 = [behav.data.Ins.morphrate1] == 0
    Ons.run{5}.sound{2} = Ons.raw{5}(idXvoc_1,:)
    
    idXins_08 = [behav.data.Ins.morphrate1] == 0.8
    Ons.run{5}.sound{3} = Ons.raw{5}(idXins_08,:)
    
    idXvoc_08 = [behav.data.Ins.morphrate1] == 0.2
    Ons.run{5}.sound{4} = Ons.raw{5}(idXvoc_08,:)
    
    idXins_06 = [behav.data.Ins.morphrate1] == 0.6
    Ons.run{5}.sound{5} = Ons.raw{5}(idXins_06,:)
    
    idXvoc_06 = [behav.data.Ins.morphrate1] == 0.4
    Ons.run{5}.sound{6} = Ons.raw{5}(idXvoc_06,:)
    
    
    % specify model
    %----------------------------------------------------------------------
    
  
    for cIDX = 1:6
        % conditions
        fmri_spec.sess(r).cond(cIDX).name     = condition{r}.sound{cIDX};
        fmri_spec.sess(r).cond(cIDX).onset    = table2array(Ons.run{r}.sound{cIDX});
        fmri_spec.sess(r).cond(cIDX).duration = stimduration;
        fmri_spec.sess(r).cond(cIDX).tmod     = 0;
        
        % parametric modulators
%         fmri_spec.sess(r).cond(cIDX).pmod(1).name = 'soundcategory';
%         fmri_spec.sess(r).cond(cIDX).pmod(1).param = Pmodsoundcategory{cIDX};
%         fmri_spec.sess(r).cond(cIDX).pmod(1).poly = 1;
%         fmri_spec.sess(r).cond(cIDX).pmod(2).name = 'ambiguity';
%         fmri_spec.sess(r).cond(cIDX).pmod(2).param = Pmodambiguity{cIDX};
%         fmri_spec.sess(r).cond(cIDX).pmod(2).poly = 2;
%         
        fmri_spec.sess(r).cond(cIDX).orth     = 1;
    end
    fmri_spec.sess(r).multi     = {''};
    fmri_spec.sess(r).regress   = struct([]);
    fmri_spec.sess(r).hpf       = 128;
    fmri_spec.sess(r).multi_reg = Motionfn{1,r};
%     mov_reg = load(Motionfn{1,r});
%     mov_reg = mov_reg(1:end,:);
%     dlmwrite(['mov_reg_short_',num2str(r),'.txt'],mov_reg,'delimiter','\t');
%     fmri_spec.sess(r).multi_reg = cellstr(spm_select('FPlistRec',pwd,['mov_reg_short_',num2str(r),'.txt']));
     
   
 
    fmri_spec.fact             = struct([]);
    fmri_spec.bases.hrf.derivs = [0 0];
    fmri_spec.volt             = 1;
    fmri_spec.global           = 'None';
    fmri_spec.mthresh          = 0.8;
    fmri_spec.mask             = {''};
    fmri_spec.cvi              = 'AR(1)';

    matlabbatch{1}.spm.stats.fmri_spec = fmri_spec;

    %----------------------------------------------------------------------
    % Model Estimation
    %----------------------------------------------------------------------
    fmri_est.method.Classical = 1;
    fmri_est.spmmat           = {[pwd,'/SPM.mat']};
    
    clear fmri_est
    fmri_est.spmmat(1)        = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    fmri_est.write_residuals  = 0;
    fmri_est.method.Classical = 1;
    
    matlabbatch{2}.spm.stats.fmri_est = fmri_est;
    
    
    %----------------------------------------------------------------------
    % Contrast Manager
    %----------------------------------------------------------------------
    cname = [];
        cons = [];
        
        % baseline
        cname{1}  = 'fem';
        cname{2}  = 'mal';
        cname{3}  = 'cla';
        cname{4}  = 'bcl';
        cname{5}  = 'vlc';
        cname{6}  = 'vla';
        
        
        cons{1}   = [1 0 0 0 0 0 zeros(1,6)];
        cons{2}   = [0 1 0 0 0 0 zeros(1,6)];
        cons{3}   = [0 0 1 0 0 0 zeros(1,6)];
        cons{4}   = [0 0 0 1 0 0 zeros(1,6)];
        cons{5}   = [0 0 0 0 1 0 zeros(1,6)];
        cons{6}   = [0 0 0 0 0 1 zeros(1,6)];
        
        cname{7}  = 'voc';
        cname{8}  = 'nvoc';
        cname{9}  = 'string';
        cname{10}  = 'ww';
        
        cons{7}   = [1 1 0 0 0 0 zeros(1,6)];
        cons{8}   = [0 0 1 1 1 1 zeros(1,6)];
        cons{9}   = [0 0 1 1 0 0 zeros(1,6)];
        cons{10}  = [0 0 0 0 1 1 zeros(1,6)];
        
   
        
        % contrasts
        
        cname{11}  = 'fem_mal';
        cname{12}  = 'mal_fem';
        cname{13} = 'str_ww';
        cname{14} = 'ww_str';
        cname{15} = 'voc_nvoc';
        cname{16} = 'nvoc_voc';
        
        cons{11}  = [1 -1 0 0 0 0 zeros(1,6)];
        cons{12}  = [-1 1 0 0 0 0 zeros(1,6)];
        cons{13}  = [0 0 1 -1 0 0 zeros(1,6)];
        cons{14}  = [0 0 -1 1 0 0 zeros(1,6)];
        cons{15}  = [1 1 -0.5 -0.5 -0.5 -0.5  zeros(1,6)];
        cons{16}  = [-1 -1 0.5 0.5 0.5 0.5  zeros(1,6)];
        
        for j = 1:size(cons,2)
            con.consess{j}.tcon.name    = cname{j};
            con.consess{j}.tcon.convec  = cons{j};
            con.consess{j}.tcon.sessrep = 'none';
        end
        nmain = numel(con.consess);
        % spm_jobman('run',matlabbatch);
        
        % make F contrasts
        %------------------------------------------------------------------
        cname = [];
        cons = [];
        cname{1} = 'eoi';
        cons{1}  = [eye(6) zeros(6,6)];
        
        for j = 1:size(cons,2)
            con.consess{j+nmain}.fcon.name = cname{j};
            con.consess{j+nmain}.fcon.convec = cons{j};
            con.consess{j+nmain}.fcon.sessrep = 'none';
        end
        con.spmmat = {[ModelPath,'/SPM.mat']};
        con.delete = 0;
        
        
        matlabbatch{3}.spm.stats.con = con;
        %spm_jobman('run',matlabbatch);
        
    
    
    %----------------------------------------------------------------------
    % Results Report
    %----------------------------------------------------------------------
    
    mSize = numel(matlabbatch);
    clear results
    for cIDX = 1:2:(numel(con.consess)*2)
        results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        if isfield(con.consess{(cIDX+1)/2},'tcon')
            results.conspec.titlestr = con.consess{(cIDX+1)/2}.tcon.name;
        elseif isfield(con.consess{(cIDX+1)/2},'fcon')
            results.conspec.titlestr = con.consess{(cIDX+1)/2}.fcon.name;
        else
            disp('ERROR')
        end
        results.conspec.contrasts = (cIDX+1)/2;
        results.conspec.threshdesc = 'none';
        results.conspec.thresh = 0.001;
        results.conspec.extent = 0;
        results.conspec.conjunction = 1;
        results.conspec.mask.image.name = cellstr(MaskPath);
        results.conspec.mask.image.mtype = 0;
        results.units = 1;
        results.print = 'pdf';
        results.write.none = 1;
        
        matlabbatch{cIDX+mSize}.spm.stats.results = results;
    end
    
    
    %----------------------------------------------------------------------
    % Print Results
    %----------------------------------------------------------------------
    
    clear print
    for cIDX = 2:2:(numel(con.consess)*2)
        print.fname = fullfile(ResultPath,sprintf('c%i_%s',cIDX/2,sList{sIDX}));
        print.fig.fighandle = NaN;
        print.opts = 'png';
        
        matlabbatch{cIDX+mSize}.spm.util.print = print;
    end
    
    %======================================================================
    % RUN matlabbatch
    %======================================================================
    
    save(fullfile(ExpPath, 'Analysis', 'Batch', ...
        ['batch_GLMEval_' sList{sIDX} '_new.mat']), ...
        'matlabbatch')
    
    spm_jobman('run',matlabbatch);
    
end

