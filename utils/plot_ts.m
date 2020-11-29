function plot_ts(ts,ts2,varargin)


%--------------VARARGIN----------------------
params  =  {   'ordering','constant','xlim','ylim','linewidth','DataCleared','filter','PolOrder','FullOrt'};
defParms = { 'sequential',     'off',    [],    [],         [],      'false',      [],        1      'off'};
legalValues{1} = {'sequential','random','stda','stdd'};
legalValues{2} = {'on','off'};
legalValues{3} = [];
legalValues{4} = [];
legalValues{5} = {'on','off'};
legalValues{6} = {'true','false'};
legalValues{7} = [];
legalValues{8} = [-1 0 1 2];
legalValues{9} = {'on','off'};
[confounds,firstmean,deri,squares,TvarNormalise,DataCleared,freq,PolOrder,FullOrt] = ParseVarargin(params,defParms,legalValues,varargin,1);
% --------------------------------------------

if nargin == 1
    ts2 = [];
end

figure; 

n = size(ts,2);


if isempty(ts2)
    for l = 1:n
        plot(ts(:,l));
        pause
    end
else
    for l = 1:n
        plot(ts(:,l)); hold on;
        pause
        plot(ts2(:,l)); 
        pause
        hold off;
    end

end