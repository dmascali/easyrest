function [filenamesout,filestructsout]=conn_dir(pathname,filter)
% CONN_DIR RECURSIVE DIR
% e.g. filenames=CONN_DIR('*.m');
%      filenames=CONN_DIR('c:/*.nii');
%

persistent filenames filestructs dodisp;
if nargin==1,% init
    [pathname,name,ext]=fileparts(pathname);
    if isempty(pathname),pathname=pwd;end
    filter=[name,ext];
    filenames=[];filestructs=[];
    dodisp=~nargout;
end
filterrest=filter;
[filternow,filterrest]=strtok(filterrest,';');
while ~isempty(filternow),
    filename=fullfile(pathname,fliplr(deblank(fliplr(deblank(filternow)))));
    dir0=dir(filename);
    [names,idx]=sortrows(strvcat(dir0(:).name));
    for n1=1:length(dir0),
        if ~dir0(idx(n1)).isdir,
            txt=fullfile(pathname,dir0(idx(n1)).name);
            if size(filenames,1)<1e5, % Change this value to increase the maximum number of filenames displayed
                filenames=strvcat(filenames,txt);
                filestructs=cat(1,filestructs,dir0(idx(n1)));
            else, return; end
            if dodisp,disp(txt);end
        end
    end
    [filternow,filterrest]=strtok(filterrest,';');
end
dir0=dir(pathname);
[names,idx]=sortrows(strvcat(dir0(:).name));
for n1=1:length(dir0),
    if dir0(idx(n1)).isdir && ~strcmp(dir0(idx(n1)).name,'.') && ~strcmp(dir0(idx(n1)).name,'..'),
        conn_dir(fullfile(pathname,dir0(idx(n1)).name),filter);
    end
end
if dodisp,filenamesout=[];filestructsout=[];else,filenamesout=filenames;filestructsout=filestructs;end

return
end