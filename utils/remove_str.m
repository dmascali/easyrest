function str = remove_str(str,toremove)
%removes the pattern TOREMOVE from the STR
indx = strfind(str,toremove);

str(indx:(indx+length(toremove)-1)) = [];

return
end