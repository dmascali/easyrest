function varargout = parse_varargin(params,defparams,legalvalues,var_arg)
%PARSE_VARARGIN(PARAMS,DEFPARAMS,LEGALPARAMS,VAR_ARG)
%parses VAR_ARG and returns a variable number of outputs (equal to
%the number of PARAMS) with the specified or default values (DEFPARMS).
%
%Performs several checks: 1) returns error if property names are not valid 
%(ie, are not among PARAMS). 2) [optinal] checks validity of property values
%(requires LEGALVALUES). 
%
%LEGALVALUES can be empty (no check #2)
%
%Usage example:
%
%   varargin = {'tail','right','delta',42,'zscore','true'};
% 
%   params  = {'tail','zscore','delta'};
%   defparms = {'both',   'off',     3 };
%
%   legalvalues{1} = {'both','right','left'};
%   legalvalues{2} = {'true','false','off','on'};
%   legalvalues{3} = [];
%
%   [tail,zscore_flag,delta] = parse_varargin(params,defparms,legalvalues,varargin)
%__________________________________________________________________________
%
%   Version 1.0
%   Copyright (C) 2018   danielemascali@gmail.com

if nargin == 0
    help parse_varargin
    return
end

n_params = length(params);

%check valid property names
l = 1;
while l <= size(var_arg,2)
    check_property_name = sum(strcmpi(params, var_arg{l}));
    if ~check_property_name 
        %print valid property names
        str = [];
        for m = 1:n_params
            str = [str,' "',params{m},'"'];
        end
        error(sprintf(['There is no "',var_arg{l},'" property.\nValid property names are:',str,'.\n']));
    end
    l = l +2;
end

for l = 1:n_params 
    inputExist = find(strcmpi(var_arg, params{l}));
    if inputExist
        %-------------check validity if desired-----------
        if ~isempty(legalvalues) && ~isempty(legalvalues{l}) 
            if iscell(legalvalues{l}) 
                valid_check = sum(strcmpi(legalvalues{l}, var_arg{inputExist+1}));
            else % numbers must not be in cells
                valid_check = sum(legalvalues{l} == var_arg{inputExist+1});
            end
            if valid_check == 0
                if iscell(legalvalues{l}) 
                    str = ['\nLegal values are:'];
                    for m = 1:length(legalvalues{l})
                        str = [str,' "',legalvalues{l}{m},'"'];
                    end
                    str = [str,'.\n'];
                else
                    str = ['\n'];
                end
                error(sprintf(['Invalid value for parameter: "',var_arg{inputExist},'".',str]));
            end
        end
        %-------------------------------------------------
        varargout{l} = var_arg{inputExist+1};
    else
        varargout{l} = defparams{l};
    end
end

return
end