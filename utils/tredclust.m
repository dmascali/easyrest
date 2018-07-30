function [struct_info] = tredclust(clip1,clip2,volume,input_name,output_name)
% Clusterize statistical maps (t,f,z-score) specifing clips values and 
% minimun cluster volume. If no output name is specified no map will 
% be printed. Cluster info will be always reported in struct_info.
%
% If mni2atlas is in the matlab path, labels information will be printed in
% struct_info and also in txt files. According to mni2atlas functionality,
% two types of anatomical labels will be printed 1) VECTOR modality, i.e.,
% anatomical label of the peak value of the cluster; 2) ROI modality, i.e.,
% anatomical composition of the cluster (the frequency of a label in the
% roi). Roi modality is available only if the output_name is specified
% (i.e., if the clusters are printed as nii file)
%
% danielemascali@uniroma1.it

%update: connection radius now is 1.5 to be consistent with spm results
% dm 29/01/2015
%update: spm for loading nii (before was load_nii)
% dm 30/01/2015
%update: clear tmp variable before running script 
% dm 05/11/2015
%update: several update, including intergration with mni2atlas.m
% dm 05/27/2015 (?)
%update: now an ALLbinary mask is also produced in output
% dm 06 feb 2017

print_labels = 1;       %was 1
if exist('mni2atlas.m')~=2
    warning off backtrace
    warning('Mni2atlas.m not found. Anatomical labels will be not printed');
    warning on backtrace
    print_labels = 0;
end

if nargin == 5
   produce_roi = true;
else
   produce_roi = false;
   if print_labels
        warning off backtrace
        warning('Anatomical labels will be printed only for peak values. If you want the anatomical composition of the entire cluster consider to print cluster maps recalling this function with the addiontal input: output_name.');
        warning on backtrace
   end
end
 
if exist(input_name)~=2
    error('Could not find input image.');
end


struct_info.name = [];
struct_info.cluster_number = [];
struct_info.mni_position= [];
struct_info.volume = [];
struct_info.peak_value = [];

[~,~]= system('rm _temp.nii');

[~,stdo] = system (['3dclust -savemask _temp.nii -orient LPI -2clip ',num2str(clip1),' ',num2str(clip2),' -dxyz=1 1.5 ',num2str(volume),' ',input_name]);


if ~isempty(dir('_temp.nii'))
       
    struct_info.name = input_name;

    %--------Parsing stdo--------------------------------------------
    
    top_limit = strfind(stdo,'.nii');
    down_limit = strfind(stdo,'#------');
    
    if down_limit(end) < top_limit(end)
        
        parsed = stdo(1,(top_limit(end) + 4):end);   %pare che se trova un solo cluster non inserisce #------
        
    else
    
        parsed = stdo(1,(top_limit(end) + 4):(down_limit(end) -1));
    
    end
    
    tab_info = str2num(parsed);
    
    tab_dim = size(tab_info);   %1 -> row, 2 -> colum
    
    if tab_dim(2)~=16
        fprintf('\n\n----------------------------\nERROR in reading std output!\n------------------------------\n\n');
    end
    
    struct_info.mni_position = tab_info(1:tab_dim(1),14:16);
    struct_info.volume = tab_info(1:tab_dim(1),1);
    struct_info.peak_value = tab_info(1:tab_dim(1),13);
    struct_info.cluster_number = tab_dim(1);
    
    fprintf(['\nCluster found for ',input_name(1:end-4),':\t\t',num2str(tab_dim(1)),'\n']);
    
    %----------------------------------------------------------------
    
    if produce_roi
        nii = spm_read_vols(spm_vol('_temp.nii'));
        clust_number = max(nii(:));
        [~,~] = system(['rm ',output_name,'_0*nii ',output_name,'_1*nii ',output_name,'_2*nii ',output_name,'_3*nii ',output_name,'_4*nii ',output_name,'_5*nii ',output_name,'_6*nii ',output_name,'_7*nii ',output_name,'_8*nii ',output_name,'_9*nii ',output_name,'_ALL.nii ',output_name,'_ALLbinary.nii']);
        for i=1:clust_number
            if i<=9
                [~,~] = system(['3dcalc -a _temp.nii -expr ''within(a,',num2str(i),',',num2str(i),')'' -prefix ',output_name,'_0',num2str(i),'.nii -overwrite']);
            else
                [~,~] = system(['3dcalc -a _temp.nii -expr ''within(a,',num2str(i),',',num2str(i),')'' -prefix ',output_name,'_',num2str(i),'.nii -overwrite']);
            end
        end
        [~,~] = system(['3dcopy _temp.nii ',output_name,'_ALL.nii']);
        [~,~] = system(['3dcalc -a _temp.nii -prefix ',output_name,'_ALLbinary.nii -expr ''step(a)'' -overwrite']);
    end    
   
    system('rm _temp.nii');

else
    fprintf('> No cluster found.\n');
    struct_info.name = input_name;
    struct_info.cluster_number = 0;
end
    

if print_labels;
    %vectro modality
    if ~isempty(struct_info.mni_position)
        temp = [];
        for j=1:size(struct_info.mni_position,1)
            tmp = mni2atlas(struct_info.mni_position(j,:));
            temp =[temp;tmp];
        end
        struct_info.label_vec = temp;
    else
        struct_info.label_vec = [];
    end
    print_labels_do(struct_info,[input_name(1:end-4),'_VEC'],clip1,clip2,volume);
    
    %roi modality
    if produce_roi
        list_rois = dir([output_name,'_*.nii']);
        temp = [];
        for j=1:(length(list_rois)-1)   %elimino _ALL
            tmp = mni2atlas(list_rois(j).name);
            temp =[temp;tmp];
        end
        struct_info.label_roi = temp;
        print_labels_do(struct_info,[input_name(1:end-4),'_ROI'],clip1,clip2,volume);
    end
end
    
return
end


function print_labels_do(a,output_name,clip1,clip2,volume)
% print on file

if strcmp(output_name(end),'C')
    %modality vect
    vect_mod = 1;
    label = a.label_vec;
else
    vect_mod = 0;
    label = a.label_roi;
end

fid = fopen([output_name,'.txt'],'w');
fprintf(fid,['Map name:\t',a.name,'\n']);
fprintf(fid,['Clip1:\t',num2str(clip1),'\n']);
fprintf(fid,['Clip2:\t',num2str(clip2),'\n']);
fprintf(fid,['MinVol:\t',num2str(volume),'\n']);
fprintf(fid,['Cluster found:\t',num2str(a.cluster_number),'\n']);
if vect_mod
    fprintf(fid,['Mni2atlas modalitiy:\t VECTOR. Following antamical labels refer to peak value.\n']);
else
    fprintf(fid,['Mni2atlas modalitiy:\t ROI. Following antamical labels refer to the average composition of the cluster.\n']);
end
fprintf(fid,'For more details see the help function of mni2atlas.m\n\n');

for j=1:a.cluster_number
    fprintf(fid,'---------------------------------------\n');
    fprintf(fid,'%d) %d\t%2.2f\t%d %d %d ',j,a.volume(j),a.peak_value(j),a.mni_position(j,:));

    for k=1:size(label,2)

        if ~isempty(label(j,k).label)

                fprintf(fid,'\n%s',label(j,k).name);

                for l=1:length(label(j,k).label)

                     fprintf(fid,'\n\t%s',label(j,k).label{l});

                end

                fprintf(fid,'\n');

        end


    end

end

fprintf(fid,'\n');

fclose(fid);

return
end


