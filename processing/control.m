function control(status,string)
global opt

if status ~= 0 
    opt.AUX.log = strvcat(opt.AUX.log,string);
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf(opt.AUX.fid,'%s\n',string');
    error(string)
end

if ~isempty(strfind(string,'ERROR'))
    opt.AUX.log = strvcat(opt.AUX.log,string);
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf(opt.AUX.fid,'%s\n',string');
    error(string)
end

opt.AUX.log = strvcat(opt.AUX.log,string);
fprintf(opt.AUX.fid,'%s\n',string');

return
end