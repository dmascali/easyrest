function labels = atlas_label_ROImode(atlas,mni2atlas_selector,varargin)
%ATLAS_LABEL_ROIMODE returns labels obtained via mni2atlas in ROI mode (see
% mni2atls for info). MNI2ATLAS_SELECTOR is a vector containing indexes 
% corrisponding to the FSL atlases (see mni2atls for info).
% The function returns the first not empty label found among the selected
% FSL atlases.
%
% Additional char input variables are:
% "SHOW_PERCENT"   [1/0, def 0] shows the frequency of that label. 
% "PRESERVE_HOLES" [1/0, def 0] In case of missing indexes in ATLAS, preserves
%                  the holes. NB: it must be set to 0 to be compatible with 
%                  ER output. 
% Usage example:
%    atlas_label_ROImode('300_ROIs_parcel.nii',[2 3],'show_percent',1,'preserve_holes',0)
%__________________________________________________________________________
%SYSTEM REQUIREMENTS
%   1) SPM
%   2) NifTI and ANALYZE tool (version > 2012-10-12)
%   3) fsl atlases
%__________________________________________________________________________
%   Version 1.0
%   Copyright (C) 2018   danielemascali@gmail.com

% TODO: instead of extracting the first not empty label, compare more
% atlases and pick the greatest. 

if nargin == 0
    help atlas_label_ROImode
    return
end

%--------------VARARGIN----------------------
params  =  {'show_percent','preserve_holes'};
defparms = {            0,               0'};
legalvalues{1} = [0 1];
legalvalues{2} = [0 1];
[show_percent,preserve_holes] = parse_varargin(params,defparms,legalvalues,varargin);
%-------------------------------------------

hdr = spm_vol(atlas);
img = spm_read_vols(hdr);

[~,name,~] = fileparts(atlas);

ind_img=setdiff(unique(img(:)),0);
if max(ind_img)~= length(ind_img)
    if preserve_holes
        warning('There is/are hole(s) (ie, missing idexes) in the atlases. ATTENTION: holes will be preserved (this is NOT how ER deals with holes!)');
    else
        warning('There is/are hole(s) (ie, missing idexes) in the atlases. They will be removed (as done in ER). BEAWRE: sequential indexing will not be reliable. Consider to use the option "preseve_holes".');
    end
end

number_of_altases = length(mni2atlas_selector);

if preserve_holes
    indx_rois = 1:max(ind_img);
    labels = cell(max(ind_img),1);
else
    indx_rois = 1:length(ind_img);
    labels = cell(length(ind_img),1);
end


for l = indx_rois
    tmp = img;
    if preserve_holes
        tmp(tmp ~= l) = 0;
    else
        tmp(tmp ~= ind_img(l)) = 0;
    end
    tmp(tmp > 0 ) = 1; if preserve_holes && sum(tmp(:)) == 0; labels{l} = 'Missing ROI'; continue; end;
    
    info.label = [];
    count = 1;
    while isempty(info.label) && count <= number_of_altases
        selector_index = mni2atlas_selector(count);
        info = mni2atlas(tmp,selector_index);
        count = count +1;
    end
    if isempty(info.label)
        labels{l} = 'Label not found';
        continue
    end
    
    if show_percent
        labels{l} = info.label{1};
    else
        indx = find(info.label{1} == '%');
        labels{l} = info.label{1}(indx+2:end);
    end
   
end

% print labels
fid = fopen([name,'.txt'],'w');
for l = indx_rois
    fprintf(fid,'%s\n',deblank(labels{l}));
end
fclose(fid);


return
end