% Enter paths of MRI experiment (pTop) and rawdata (pRaw, subfolder including
% subfolders of participants). The script loops through all participants,
% converts the rec files to nifti and moves the converted files (converted into raw-participant folders) into the pOut (nifti folder in pTop, one fore each subject)
%
% CAVE: for conversion you first need to follow the instructions to
% install pati (within pati_linux) folder, i.e. install active Perl etc.
% This might cause problem with newer versions of active Perl (do not
% require XML anymore (and ppm not possible); There is a new version of the
% rec2nifti perl-file, (24/04/2019), which you need to add to the
% pati_linux folder.
%__________________________________________________________________________
% (C) Teresa Wenhart, Cognitive and Affective Neuroscience, Zurich, 2020
%

clear
clc

%directories
pTop = 'C:\Users\teresa\Documents\projects\ambiguous_voices\FMRI_study\data\MRI';
pRaw = 'C:\Users\teresa\Documents\projects\ambiguous_voices\FMRI_study\data\MRI\raw';

%loop through participants.
D = dir(pRaw); %subfolders of participants
for k = 3:length(D) % avoid using the first ones
    subj = D(k).name; %
    
    pSubj = fullfile(pRaw, subj);
    pOut = fullfile(pTop,'nifti');
    pOut = fullfile(pOut, subj);
    pRec = fullfile(pTop, 'rec');
    pRec = fullfile(pOut, subj);


    rFiles = dir(fullfile(pSubj,'*.par'));
    rFiles = {rFiles.name};

    tvaevent_convert_rec(pSubj) %WARNING OUTPUIT INTO RAW DIR

    %move nii

    conds = {'baseline','ins','tim','voc','gen','t1','tva','b0'}; %Seperate out scan things _you can make this better!


    niis = dir(fullfile(pSubj,'*.nii'));
    niis = {niis.name};

    logs = dir(fullfile(pSubj,'*.log'));
    logs = {logs.name};

    recs = dir(fullfile(pSubj,'*.rec'));
    recs = {recs.name};

    for iC = 1:length(conds);
        pCond = fullfile(pOut,conds{iC});
        %make outdir
        mkdir(pCond);
        idxNii = find(contains(niis,conds{iC})); %index of niis for iC condition.

        condNiis = niis(idxNii); %relevant Niis
        for iN = 1:length(condNiis);
            %move files
            movefile(fullfile(pSubj,condNiis{iN}),pCond);
        end

        %move physio files
        pLogOut = fullfile(pSubj,'physio');
        mkdir(pLogOut);  
        idxLog = find(contains(logs,conds{iC}));

        condLog = logs(idxLog); %logs to be moved.

        for iL = 1:length(condLog);
            %move files
            movefile(fullfile(pSubj,condLog{iL}),pLogOut);
        end  
        
        idxRec = find(contains(recs,conds{iC}));
        condRec = recs(idxRec);
        for iR = 1:length(condRec);
            movefile(fullfile(pSubj,condRec{iR}),pRec);
        end
        
    end
end
