function [str] = c_b(str)    %c_b stands for cut blanks
%remove blank space at end of the string (issue due to conn_dir)
str(str==' ') = '';
return
end