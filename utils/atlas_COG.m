function [atlas_COG] = atlas_COG(atlas)
% Find center of gravity of atlas and print to a text file
%
%   Remarks:
%       If the true center of gravity is outside the region, the voxel 
%       inside the region that is closest to this point is used as the 
%       actual center of gravity. 
%
%       MNI coordinates are reported. 
%
%   adjusted from NBS toolbox

% checked left right flipping.


hdr = spm_vol(atlas);
img = spm_read_vols(hdr);

%find origin
mat = hdr.mat;
origin = mat(1:3,4); 
vox_size = [hdr.mat(1,1);hdr.mat(2,2);hdr.mat(3,3)];
origin = abs(origin./vox_size)';

ind_img=setdiff(unique(img(:)),0);

coor=zeros(length(ind_img),3);
for i=1:length(ind_img)
    [x,y,z]=ind2sub(size(img),find(ind_img(i)==img));
    [val,ind_min]=min(sqrt((mean(x)-x).^2+(mean(y)-y).^2+(mean(z)-z).^2)); 
    coor(i,:)=[x(ind_min(1)),y(ind_min(1)),z(ind_min(1))]; 
end
%Map voxel coordinates to MNI coordinates
coor=(coor-repmat(origin,length(ind_img),1));
coor = bsxfun(@times,coor,vox_size');


%Write to a text file
[~,name,~] = fileparts(atlas);
dlmwrite([name,'_COG.txt'],coor,'delimiter',' ','precision','%.3f');
if nargout > 0
    atlas_COG = coor;
end

return
end

