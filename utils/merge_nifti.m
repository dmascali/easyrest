function merge_nifti(folder,wildcard,output_name,silent_mode)
%MERGE_NIFTI(FOLDER,WILDECARD,OUTPUT_NAME)
%Merge, using fslmerge, all nifti (or nii.gz) defined by WILDCARD in FOLDER 
%and in each subfolders. Merged file has OUTPUT_NAME. The actual merged
%file list will be output as a .txt file with the same output_name.
%   example:
%   merge_nifti('/home/Desktop/','epi*.nii','merged_epi')
%
%MERGE_NIFTI(.,SILENT_MODE)
%Use SILENT_MODE flag to suppres stdo.
%   example:
%   merge_nifti('/home/Desktop/','epi*.nii','merged_epi',1)
%_______________________________________________________________________
% Copyright (C) 2017 danielemascali@gmail.com

if nargin < 3
    help merge_nifti
    return
end

if nargin == 4
    silent_mode = 1;
else
    silent_mode = 0;
end

list = conn_dir([folder,'/',wildcard]);
total = size(list,1);
if total == 0
    fprintf('\nNo file found.\n');
    return
end
% print on file
fid = fopen([output_name,'.txt'],'w');
fprintf(fid,'Total nifti file merged: %d.',total);
for l = 1:total
    fprintf(fid,'\n%s',c_b(list(l,:)));
end
fclose(fid);
%------------------------------------

% print on screen if desired
if ~silent_mode
    for l = 1:total
        fprintf('\n%s',c_b(list(l,:)));
    end
    fprintf('\nThese files are going to be merged (n = %d).',total);
end
%------------------------------------


merged_str = [];
for l = 1:size(list,1)
    merged_str = [merged_str,c_b(list(l,:)),' '];
end

system(['fslmerge -t ',output_name,' ',merged_str]);


if ~silent_mode
    fprintf('Completed. Bye.\n');
end

return
end


function [str] = c_b(str)    %c_b stands for cut blanks
%remove blank space at end of the string (issue due to conn_dir)
str(str==' ') = '';
return
end