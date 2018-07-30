function [Y] = BandPassOrt(N,TR,F1,F2,invert)
% similar to AFNI function 1dBport
% the constant component is not produced


%define frequency grid
deltaf = 1/(N*TR);
nyquist = 1/(2*TR);
freq = deltaf:deltaf:nyquist;


index_selector = ( F1 <= freq & freq <= F2);

if invert
    index_selector = not(index_selector);
end

Nreg = sum(index_selector);

if Nreg == 0
    Y = [];
    disp('No frequency found.');
    return
end

index = find(index_selector);

x = 0:1:(N-1);
t = x.*TR;

Y = [];
for l = 1:Nreg
    Y = [Y,cos(2*pi*freq(index(l))*t)'];
    Y = [Y,sin(2*pi*freq(index(l))*t)'];
end

% The sine of the nyquist frequency must be removed 
if freq(index(end)) == nyquist
    Y = Y(:,1:end-1);
end

%figure; plot(Y);

return
end
