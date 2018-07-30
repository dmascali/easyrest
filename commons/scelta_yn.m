function [opzione]= scelta_yn (stringa,default)
%scelta_yn.m Sceglie tra y e n
%Versione 1.2
%Federico Giove, ultima modifica 30/09/2002
%Versione 1.3 : ora se si inserisce una stringa con length > 1 non
% da più errore. daniele mascali 18/02/2015

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
  