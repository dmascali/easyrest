function smooth_after_processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       editable fields:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list_measures{1} = 'ReHo';  % put the name of the first level folders to be smoothed

FWHM = 4;

%type = '3dBlurInMask';   % type of smoothing
type = '3dBlurToFWHM';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -global opt
[fname_UNSM, pathname_UNSM] = uigetfile('*.mat', 'Select the .mat of the unsmoothed project');
[fname_SM, pathname_SM] = uigetfile('*.mat', 'Select the .mat of the final project (smoothed)');

unsm = load ([pathname_UNSM,'/',fname_UNSM]);
sm = load ([pathname_SM,'/',fname_SM]);
wd = pwd;

count = 0;
for r = 1:length(list_measures)
    pathname_UNSM = [unsm.opt.folders.firstlevel,'/',list_measures{r}];
    if ~exist(pathname_UNSM,'dir')
        continue
    end
    cd (pathname_UNSM);
    list_nii = dir('*.nii');
    if isempty(list_nii)
        warning(['Measure ',list_measures{r},' has empty directory']);
        continue
    end
    count = count +1;
    % remove eventually _s_*.nii files (previous execution)
    system(['rm ',pathname_UNSM,'/_s_*.nii']);
    
    % find LFF series if 3dBlurToFWHM is chosen----------------------------
    if strcmpi(type,'3dBlurToFWHM')
        lff_dir = unsm.opt.folders.preprocessing;
        list_lff = dir([lff_dir,'/*_LFF.nii']);
    end
    %----------------------------------------------------------------------
    
    % now work on the target (smoothed) project----------------------------
    pathname_SM = [sm.opt.folders.firstlevel,'/',list_measures{r}];
    if exist(pathname_SM,'dir')
        temp_list = dir([pathname_SM,'/*.nii']);
        if ~isempty(temp_list)
            choose = choose_yn('Destination directory not empty. Next step is to remove these data, shall I continue?','y');
            if choose == 'n'
                return
            else
                system(['rm ',pathname_SM,'/*.nii']);
            end
        end
    else
        mkdir(sm.opt.folders.firstlevel,list_measures{r})
    end
    %----------------------------------------------------------------------
    h = waitbar(0,['Smoothing ',list_measures{r}]);
    for l = 1:length(list_nii)
        switch type
            case {'3dBlurInMask'}
                system(['3dBlurInMask -input ',list_nii(l).name,' -prefix _s_',list_nii(l).name,' -FWHM ',num2str(FWHM),' -mask ',unsm.opt.AUX.bmask_path{1,1}]);
            case {'3dBlurToFWHM'}
                found = 0;
                for z = 1:length(list_lff)
                    if ~isempty(strfind(list_nii(l).name,list_lff(z).name(1:end-8)));
                        found = 1;
                       break
                    end
                end; if found == 0; error('Can''t find LFF.nii');end;
                system(['3dBlurToFWHM -input ',list_nii(l).name,' -blurmaster ',lff_dir,'/',list_lff(z).name,' -prefix _s_',list_nii(l).name,' -FWHM ',num2str(FWHM),' -mask ',unsm.opt.AUX.bmask_path{1,1}]);
        end
        system(['3dcopy _s_',list_nii(l).name,' ',pathname_SM,'/',list_nii(l).name]);
        waitbar(l/length(list_nii),h);
    end
    close (h);
end

cd (wd);

return
end