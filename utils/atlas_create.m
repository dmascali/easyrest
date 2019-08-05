function atlas_create(roi_list_file,radius,template,output_name,extra_flag)
%ATLAS_CREATE ( ROI_LIST_FILE, RADIUS, TEMPLATE, OUTPUT_NAME )
%produces an atlas volume (3D nifti) from a list of ROI coordinates in MNI
%space and a mat file with (Euclidian) distances among ROIs. 
%Spherical ROIs are created. 
%If you need to create an atlas from a list of nifti files (instead of 
%list of coordinates), this functionality is embedded in easyrest 
%(use ROItoROI with ER.others.RtoR.roi_folder_or_file = 'folder_path')
% 
%Input variables:
% ROI_LIST_FILE = txt file with MNI coordinates on each row.
% RADIUS        = desired radius in mm for each ROI.
% TEMPLATE      = a 3D image (in MNI) that matches the desire final resolution.
% OUTPUT_NAME   = string name of the atlas and distance file.
%_______________________________________________________________________
% Copyright (C) 2016 danielemascali@gmail.com

%NOTE: october 2017: checked for left right flipping: NO problems.

if nargin < 4
    help create_atlas
    error('Not enough input.');
end

if nargin == 5
    call_from_ER = 1;
else
    call_from_ER = 0;
end

%get rois coordinates
fid = fopen(roi_list_file,'r');
line = fgets(fid);
count = 0;
while ischar(line)    
    count = count +1;
    roi_mni(count,:) = str2num(line); 
    line = fgets(fid);
end
fclose(fid);

% get distances among rois
d = pdist(roi_mni);
d= squareform(d);

if call_from_ER == 0
    h = waitbar(0,'Creating ROIs, please wait...');
end
% create rois
for l = 1:size(roi_mni,1)
    c = sphere_eq_string(roi_mni(l,:));
    if l < 10
        out = ['_temp_00',num2str(l),'.nii'];
    elseif l < 100
        out = ['_temp_0',num2str(l),'.nii'];
    elseif l < 1000
        out = ['_temp_',num2str(l),'.nii'];
    end
    if call_from_ER == 0
        system(['3dcalc -a ',template,' -LPI -expr ''step(',num2str(radius*radius),'-(x',c{1},')*(x',c{1},')-(y',c{2},')*(y',c{2},')-(z',c{3},')*(z',c{3},'))'' -prefix ',out,' -datum short']);
    else
         [~,~]= system(['3dcalc -a ',template,' -LPI -expr ''step(',num2str(radius*radius),'-(x',c{1},')*(x',c{1},')-(y',c{2},')*(y',c{2},')-(z',c{3},')*(z',c{3},'))'' -prefix ',out,' -datum short']);
    end
    if l == 1
        hdr = spm_vol(out);
        img = spm_read_vols(hdr);
        if sum(img(:)) == 0
            system('rm _temp_*.nii');
            error('Empty roi at the first iteration');
        end
    else
        tmp = spm_read_vols(spm_vol(out));
        if sum(tmp(:)) == 0
            system('rm _temp_*.nii');
            error(['Empty roi at interation number: ',num2str(l)]);
        end
        img = img + l*tmp;
    end
    if max(img(:)) > l
        close(h);
        system('rm _temp_*.nii');
        error(['Overlapping rois at interation number: ',num2str(l)]);
    end
    if call_from_ER == 0
        waitbar(l/length(roi_mni));
    end
end
if call_from_ER == 0
    close(h);
end
[~,name,~] = fileparts(output_name);
output_name = [name,'.nii'];
hdr.fname = output_name;
hdr.private.dat.fname = output_name;
spm_write_vol(hdr,img);
system('rm _temp_*.nii');
save([output_name(1:end-4),'_distance.mat'],'d','roi_mni');

return
end


function [c_str] =  sphere_eq_string(c)

c_str = cell(3,1);
for l = 1:3
    if c(l) >= 0
        c_str{l} = ['-',num2str(c(l))];
    else
        c_str{l} = ['+',num2str(abs(c(l)))];
    end
end
return
end
