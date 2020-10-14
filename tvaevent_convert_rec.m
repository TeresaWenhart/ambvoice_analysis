function tvaevent_convert_rec(datapath)

% enter path to .rec files from Philips scanner to convert them to analyze
% or nifti format
%
% example:
% datapath = 'D:\data\rawMRfiles\';
% convert_rec(datapath)

% substrings of file names to identify EPI, B0 and T1 files

% if all your rec files are in the same folder, these wild cards can be
% used to find them
% EPI, T1 and B0 should be identified separately
EPI_name = '*epi*.rec';  % this should be unique to your EPI files
B0_name  = '*b0*.rec';  % this should be unique to your field map files
T1_name  = '*t1*.rec';   % this should be unique to your anatomical files

% convert EPI

% find all files in the data path containing the string defined above
filename_EPI = dir(fullfile(datapath,EPI_name));
filename_EPI = {filename_EPI.name}';

for epi_idx = 1:numel(filename_EPI)
    a=1;
% -s -vol will create 3D nifty for each volume instead of one big 4D nifti
% file
% the ampersand & at the end of the command will allow to execute multiple
% conversions in parallel (it's just faster then waiting for each one to
% finish)
command = ['C:\Perl64\bin\perl.exe C:\pati_linux\rec2nifti.pl -s -vol -f ' fullfile(datapath,filename_EPI{epi_idx}) ' &'];

% execute the command (does exactly the same as running the line in Windows cmd)
[status,cmdout]= system(command,'-echo');

pause(15)

end

% convert B0

filename_B0 = dir(fullfile(datapath,B0_name));
filename_B0 = {filename_B0.name}';

for epi_idx = 1:numel(filename_B0)

command = ['C:\Perl64\bin\perl.exe C:\pati_linux\rec2nifti.pl -split -f ' fullfile(datapath,filename_B0{epi_idx}) ' &'];

[status,cmdout]= system(command,'-echo');

end

% convert T1

filename_T1 = dir(fullfile(datapath,T1_name));
filename_T1 = {filename_T1.name}';

for epi_idx = 1:numel(filename_T1)

% no -s or -vol for T1, since it's one big 3D file
command = ['C:\Perl64\bin\perl.exe C:\pati_linux\rec2nifti.pl ' fullfile(datapath,filename_T1{epi_idx}) ' &'];

[status,cmdout]= system(command,'-echo');

end
