function denoising_analysis_plot_2ndl
global opt


% sessions comparisons for each group. requires s > 1
if opt.session_number > 1
    for g = 1:opt.group_number
        data = create_data(opt.session_number,g,'group');
        title_str = ['Session comp.:\n ',opt.AUX.group_names_prefix{g}];
        if opt.group_number == 1
            output_name = ['session_comparison'];
        else
            output_name = ['session_comparison_',opt.AUX.group_names_prefix{g}(1:end-1)];
        end
        for l = 1:length(opt.AUX.session_names_prefix)
            sample_names{l} = remove_u (opt.AUX.session_names_prefix{l});
        end

        da_2ndl_plot(data,title_str,sample_names,output_name,1,0)
    end
end

% group comparison. requires g > 1. 
sample_names = [];
if opt.group_number > 1
    for j = 1:opt.session_number
        data = create_data(opt.group_number,j,'session');
        title_str = ['Group comp.:\n ',opt.AUX.session_names_prefix{j}];
        if opt.session_number == 1
            output_name = ['group_comparison'];
        else
            output_name = ['group_comparison',opt.AUX.session_names_prefix{j}];
        end
        for l = 1:length(opt.AUX.group_names_prefix)
            sample_names{l} = remove_u (opt.AUX.group_names_prefix{l});
        end
        da_2ndl_plot(data,title_str,sample_names,output_name,0,0)
    end
    %todo
end

return
end



function data = create_data(samp_numb,indx,str)         % sample are the number of entity that you want to compare 

global opt 
data.fd_TRnorm = create_data_aux(samp_numb,indx,opt.DA_extracted.fd_TRnorm,str);
data.rms = create_data_aux(samp_numb,indx,opt.DA_extracted.rms,str);
data.rms_dvars = create_data_aux(samp_numb,indx,opt.DA_extracted.rms_dvars,str);
data.rms_sd = create_data_aux(samp_numb,indx,opt.DA_extracted.rms_sd,str);
for r = 1:size(opt.DA_extracted.vTv,4)
    data.vTv{r} = create_data_aux(samp_numb,indx,opt.DA_extracted.vTv,str,r);
    data.vTv_mean{r} = create_data_aux(samp_numb,indx,opt.DA_extracted.vTv_mean,str,r);
    data.var_exp{r}= create_data_aux(samp_numb,indx,opt.DA_extracted.var_exp,str,r);
    for s = 1:size(opt.DA_extracted.specifcity,5)
        data.specificity{r,s}= create_data_aux(samp_numb,indx,opt.DA_extracted.specifcity,str,r,s);
    end
end
return
end

function Y = create_data_aux(samp_numb,indx,y,str,reg_n,sp)
global opt
switch str
    case 'group'
        select_group = 1;
    case 'session'
        select_group = 0;
end

if select_group
    for s = 1:samp_numb
        tmp = [];
        for k = 1:opt.subject_number(indx)
            if nargin == 5
                tmp = [tmp;y{indx,s,k,reg_n}];
            elseif nargin == 6
                tmp = [tmp;y{indx,s,k,reg_n,sp}];
            else
                tmp = [tmp;y{indx,s,k}];
            end
        end
        Y{s} = tmp;
    end
else
    for s = 1:samp_numb
        tmp = [];
        for k = 1:opt.subject_number(s)
            if nargin == 5
                tmp = [tmp;y{s,indx,k,reg_n}];
            elseif nargin == 6
                tmp = [tmp;y{s,indx,k,reg_n,sp}];
            else
                tmp = [tmp;y{s,indx,k}];
            end
        end
        Y{s} = tmp;
    end
end


return
end

% function data = create_data_gc(samp_numb,gindx)
% global opt 
% data.fd = create_data_aux(samp_numb,gindx,opt.DA_extracted.fd);
% data.rms = create_data_aux(samp_numb,gindx,opt.DA_extracted.rms);
% data.rms_dvars = create_data_aux(samp_numb,gindx,opt.DA_extracted.rms_dvars);
% data.rms_sd = create_data_aux(samp_numb,gindx,opt.DA_extracted.rms_sd);
% for r = 1:size(opt.DA_extracted.vTv,4)
%     data.vTv{r} = create_data_aux(samp_numb,gindx,opt.DA_extracted.vTv,r);
%     data.vTv_mean{r} = create_data_aux(samp_numb,gindx,opt.DA_extracted.vTv_mean,r);
%     data.var_exp{r}= create_data_aux(samp_numb,gindx,opt.DA_extracted.var_exp,r);
% end
% return
% end


function [str] = remove_u(str)    %r_u stands for replace underscore
%replace blank space with '_' 
str(str=='_') = '';
return
end

function [str] = add_u(str)    %r_u stands for replace underscore
%replace blank space with '_' 
str(str==' ') = '_';
return
end