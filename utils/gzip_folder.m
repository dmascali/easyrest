function gzip_folder(folder,wildcard,mode,script_mode)
%GZIP_FOLDER(FOLDER,WILDECARD,MODE)
%compress or decompress files defined by WILDCARD in FOLDER and in each
%subfolders. 
%MODE=1 COMPRESS
%MODE=0 DECOMPRESS
%   example:
%   gzip_folder('/home/Desktop/','*.nii',1)
%   gzip_folder('/home/Desktop/','*.nii.gz',0)
%
%_______________________________________________________________________
% Copyright (C) 2016 danielemascali@gmail.com

if nargin < 3
    help gzip_folder
    return
end

if nargin < 4
    script_mode = 0;
end

if mode == 1
    str = 'Compressing files...';
    str_sys = '';
elseif mode == 0
    tmp = wildcard(end-1:end);
    if  ~strcmp(tmp,'gz')
        error('Use as a wildcard for decompression only filetypes *.gz');
    end
    str = 'Decompressing files...';
    str_sys = ' -d';
else
    error('Modality bad defined. Use 1 for compress and 2 for decompress.');
end

list = conn_dir([folder,'/',wildcard]);
total = size(list,1);
if total == 0
    fprintf('\nNo file found.\n');
    return
end
for l = 1:total
    fprintf('\n%s',c_b(list(l,:)));
end
fprintf('\nThese files are going to be compressed/decompressed.');
if ~script_mode
    [opzione] = choose_yn ('Press a key to continue (type ''n'' to exit)','y');
    if opzione == 'n'
        fprintf('Aborted. Bye.\n');
        return
    end
end
if mode == 1 %computing filesize before compression
    dimension = [];
    for l = 1:total
        tmp = dir(c_b(list(l,:)));
        dimension = [dimension,tmp.bytes];
    end
    total_dimension = sum(dimension);
end
    
if ~script_mode
    h = waitbar(0,str);
end
for l = 1:total
    %[~,~] = system(['gzip ',c_b(list(l,:)),str_sys]);
     system(['gzip -f ',c_b(list(l,:)),str_sys])
    if ~script_mode
        waitbar(l/total);
    end
end

if mode == 1 %computing filesize after compression
    dimension_after = [];
    for l = 1:total
        tmp = dir([c_b(list(l,:)),'.gz']);
        dimension_after = [dimension_after,tmp.bytes];
    end
    total_dimension_after = sum(dimension_after);
    fprintf('\nSpace occupied before compression: \t%.2f GB',total_dimension/10^9);
    fprintf('\nSpace occupied after  compression: \t%.2f GB',total_dimension_after/10^9);
    fprintf('\nFree space gained:  \t\t\t%.2f GB\n',(total_dimension - total_dimension_after)/10^9);
end    

if ~script_mode
    close(h);
end



fprintf('Completed. Bye.\n');

return
end

function [str] = c_b(str)    %c_b stands for cut blanks
%remove blank space at end of the string (issue due to conn_dir)
str(str==' ') = '';
return
end

function [opzione]= choose_yn (stringa,default)
%choose yes or no
%Federico Giove

disp(' ')
opzione=' ';
while ( ~strcmp(opzione,'y') && ~strcmp(opzione,'n') && ~strcmp(opzione,'Y') && ~strcmp(opzione,'N')) || length(opzione)>1
    opzione=input ([stringa,' (y/n, default ',default,')? '],'s');
    if isempty(opzione)
        opzione=default;
    end
    if ( ~strcmp(opzione,'y') && ~strcmp(opzione,'n') && ~strcmp(opzione,'Y') && ~strcmp(opzione,'N')) || length(opzione)>1
        disp('Reply with "y" or "n"!')
        disp(' ')
    end
end
if opzione=='Y'
    opzione='y';
end
if opzione=='N'
    opzione='n';
end
return
end
  
  