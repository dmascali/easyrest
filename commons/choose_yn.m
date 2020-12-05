function [opzione]= choose_yn (stringa,default)
%choose yes or no
%Federico Giove

disp(' ')
opzione=' ';
while ( ~strcmp(opzione,'y') && ~strcmp(opzione,'n') && ~strcmp(opzione,'Y') && ~strcmp(opzione,'N')) || length(opzione)>1
    opzione=input ([stringa,' (y/n, default ',default,')? '],'s');
    if isempty(opzione)
        opzione=default;
    end
    if ( ~strcmp(opzione,'y') && ~strcmp(opzione,'n') && ~strcmp(opzione,'Y') && ~strcmp(opzione,'N')) || length(opzione)>1
        disp('Reply with "y" or "n"!')
        disp(' ')
    end
end
if opzione=='Y'
    opzione='y';
end
if opzione=='N'
    opzione='n';
end
return
end
  