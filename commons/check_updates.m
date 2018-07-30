function [exit] = check_updates(version)
exit = 0;

fprintf('\nChecking for updates...');
url_version_history = 'https://dl.dropboxusercontent.com/u/21956743/easyrest_released/easyrest_version.html';
url_for_download = 'https://dl.dropboxusercontent.com/u/21956743/easyrest_released/ER_toolbox.zip';


% Get the latest version available:
%[NewVersion,status] = urlread(url_version_history,'Timeout',10);   %
%timeout non è presente nelle versioni più vecchie di matlab
[NewVersion,status] = urlread(url_version_history); 

% Check if latest version is newer than this version:
if status~=0 && str2double(version)<str2double(NewVersion(1:3))
    fprintf('a new version has been released!\n');
    fprintf('\nUpdate History:\n');
    fprintf('%s',NewVersion(4:end));
    fprintf('\nYou can download the latest version from:');
    fprintf('\n%s\n',url_for_download);
    [opzione] = scelta_yn('Do you want to continue without updating?','n');
    if opzione == 'n'
        exit = 1;
        return
    end
elseif status == 0
    fprintf('connection to the repository failed.\n--Check your internet connection! \n--Otherwise, if you are correctly connected report this bug to danielemascali@gmail.com\n');
else
    fprintf('you are running the latest version!\n');
end

return
end