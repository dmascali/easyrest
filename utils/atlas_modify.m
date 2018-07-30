function atlas_modify(atlas,selector,output_name,binary_flag)
%ATLAS_MODIFY(ATLAS,SELECTOR,OUTPUT_NAME)
% Given a 3D-brain ATLAS, ATLAS_MODIFY creates another atlas with only the
% lables indicated by SELECTOR (an integer vector). The name of the new 
% atlas is given by OUTPUT_NAME.
% The function also creates a new .txt label file.
%
%ATLAS_MODIFY(ATLAS,SELECTOR,OUTPUT_NAME,BINARY_FLAG)
% If BINARY_FLAG is provided (ie: 1) only a binary mask (ie, a ROI) with 
% the selected lables will be in output (useful to combine several labels
% in a single mask).
%
% Examples:
%   For constructing a new atlas:
%       atlas_modify('atlas.nii',[2 3 4 7 10],'new_atlas.nii')
%
%   For constructing a single mask:
%       atlas_modify('atlas.nii',[2 3 4 7 10],'visual_roi.nii',1)
%_______________________________________________________________________
% Copyright (C) 2015 danielemascali@gmail.com


if nargin == 4
   binary_mode = true;
else
   binary_mode = false;
end
 

hdr = spm_vol(atlas);
img = spm_read_vols(hdr);
try
    fid = fopen([atlas(1:end-3),'txt']);
    line = fgets(fid);
    txt_flag = 1;
catch
    warning(sprintf('No txt file with labels found. \nSkipping creation of updated label file.'));
    txt_flag = 0;
end

    

if txt_flag
    roi_names = [];
    while ischar(line)    
        roi_names = strvcat(roi_names,line); 
        line = fgets(fid);
    end
    fclose(fid);
end

new_img = zeros(size(img));


count = 0;
for l = selector
    count = count +1;
    if binary_mode
        new_img (img == l) = 1;
    else
        new_img (img == l) = count;
    end
end

if txt_flag
    fid = fopen([output_name(1:end-3),'txt'],'w');
    for l = selector
        fprintf(fid,'%s\n',deblank(roi_names(l,:)));
    end
    fclose(fid);
end

hdr.pinfo(1) = 1; hdr.pinfo(2) = 0;
hdr.fname = output_name;
hdr.private.dat.fname = output_name;
spm_write_vol(hdr,new_img);


return
end
