function out = char2logical(inp)
%convert char varargin to logical
switch inp
    case {'on','true'}
        out = 1;
    case {'off','false'}
        out = 0;
end
out = logical(out);
return
end

