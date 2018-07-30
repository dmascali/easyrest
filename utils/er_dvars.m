function er_dvars(input,output_name,method,mask)
% method 'dvars' 'srms' see 3dTto1D for additional options and info
% if mask is empty then use autmoask

if isempty(mask)
    mask_str = ' -automask';
else
    mask_str = [' -mask ',mask];
end

system(['3dTto1D -input ',input,mask_str,' -method ',method,' -prefix ',output_name]);


return
end