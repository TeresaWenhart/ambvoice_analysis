%script to rename MRI files according to Sascha Fr?hholz standard
%file prefixes (replace "p2"?):
%   "subj" = ID of subject e.g. av01 (taken from participants subfolder
%   *b0map* --> b0_subj +...
%   *epi* --> f_subj +...
%   *survey* --> f_subj +...
%   *t1* /*anat* --> s_subj + ...
%__________________________________________________________________________
% (C) Teresa Wenhart, Cognitive and Affective Neuroscience, Zurich, 2020
%
%
%
%
clear
clc

%directories
pTop = 'C:\Users\teresa\Documents\projects\ambiguous_voices\FMRI_study\data\MRI\nifti';

%loop through participants.
D = dir(pTop); %subfolders of participants
for k = 3:length(D) %exclude empty directories ., ..
    subj = D(k).name; %
    pSubj = fullfile(pTop, subj);
    conds = dir(pSubj); %usually: conds = {'baseline','ins','tim','voc','gen','t1','tva','b0', 'physio'};
    for iC = 4:length(conds); %exclude empty directories ., ..
        pCond = fullfile(pSubj,conds(iC).name);
        % find specific files
        b0s = dir(fullfile(pCond,'*b0map*'));
        b0s = {b0s.name};
        survs = dir(fullfile(pCond,'*survey*'));
        survs = {survs.name};
        epis = dir(fullfile(pCond,'*epi*'));
        epis = {epis.name};
        anats = dir(fullfile(pCond,'*anat*'));
        anats = {anats.name};
        % add prefix to filenames
        nb0s = strcat('b0_',subj, '_', b0s);
        nepis = strcat('f_',subj, '_', epis);
        nsurvs = strcat('f_',subj,'_', survs);
        nanats = strcat('s_',subj,'_', anats);
        nanatsOUT = strcat('s_',subj,'_t1.nii');
        %rename
         
        for j = 1: length(b0s);
            movefile(fullfile(pCond,b0s{j}),fullfile(pCond, nb0s{j}));
            end
        for e = 1: length(epis);
            movefile(fullfile(pCond,epis{e}),fullfile(pCond, nepis{e}));
            end
        for s = 1: length(survs);
            movefile(fullfile(pCond,survs{s}),fullfile(pCond, nsurvs{s}));
            end
        for a = 1: length(anats);
            movefile(fullfile(pCond,anats{a}),fullfile(pCond, nanatsOUT));
            end
    end
end

