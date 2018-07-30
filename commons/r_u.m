function [str] = r_u(str)    %r_u stands for replace underscore
%replace '_' with blank space 
str(str=='_') = ' ';
return
end