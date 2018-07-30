function denoising_analysis_plot_specificity_timeseries(DA,s)

global opt

width = 12.5;     % Width in inches
height =7;    % Height in inches
%alw = 0.75;    % AxesLineWidth
%fsz = 11;      % Fontsize
lw = 1;      % LineWidth
lw_leg = 3;
lw_leg2 = 2;
msz = 9;       % MarkerSize
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);
fntunit = 'points'; % default 'points','inches','centimeters','normalized','pixels'
text_.fnts = 12;
text_.fntw = 'bold';
FntsLeg = 9;

figure('visible','off');
clear h h2 h3 hl;

colors = opt.DA.colors;

colors_ts = [255 51 51;
             0 128 255;
             0 0 0]./255;

sp_r = 11;
sp_c = 5;

subplot(sp_r,sp_c, [(1)],'Visible','off')
%ax = axes('Position',[0 0 1 1],'Visible','off');
str = {['Proj. name: ',opt.folders.prject_name];
          ['Subj name: ',r_u(DA.subj_name)];
          ['Mean FD:              ',num2str(DA.rp_based.fd_mean),' mm'];
          ['Mean FD/tr:          ',num2str(DA.rp_based.fd_TRnorm_mean),' mm'];
          ['RMS movement:  ',num2str(DA.rp_based.rms_mean),' mm'];
          ['Network:       ',opt.DA.roi_specificity(s).name];
          };
text(0,0.7,str,'interpreter','none','FontUnits',fntunit,'FontSize',text_.fnts,'FontWeight',text_.fntw);  %was 0,0.7

% print the legend of the regression
subplot(sp_r,sp_c, [(2)],'Visible','off')
hold on;
h{1} = plot([0 0],[0 0],'Color',colors(1,:),'LineWidth', lw_leg);
set(gca,'visible','off');
leg_str{1,1} = DA.reg(2).name_contract;
legg = h{1};
for r = 3:length(DA.reg);
    h{r-1} =  plot([0 0],[0 0],'Color',colors(r-1,:),'LineWidth', lw_leg);
    set(gca,'visible','off');
    leg_str{1,r-1} = DA.reg(r).name_contract;
    legg = [legg,h{r-1}];
end
hl = legend(legg,leg_str);
set(hl, 'Position',[0.35, 0.65, 0.01, 0.07],'box','off','FontSize',FntsLeg);

% print the legend of the timeseries
clear legg leg_str h
subplot(sp_r,sp_c, [(sp_c +1)],'Visible','off')
hold on;
h{1} = plot([0 0],[0 0],'Color',colors_ts(1,:),'LineWidth', lw_leg2);
set(gca,'visible','off');
h{2} =  plot([0 0],[0 0],'Color',colors_ts(2,:),'LineWidth', lw_leg2);
set(gca,'visible','off');
h{3} =  plot([0 0],[0 0],'--','Color',colors_ts(3,:),'LineWidth', lw_leg2);
set(gca,'visible','off');
legg = [h{1},h{2},h{3}];
if length(DA.reg(2).specificity(s).names.all) == 4
    leg_str{1,1} = DA.reg(2).specificity(s).names.all{1}; leg_str{1,2} = DA.reg(2).specificity(s).names.all{2}; leg_str{1,3} = ['<',DA.reg(2).specificity(s).names.all{3},',',DA.reg(2).specificity(s).names.all{4},'>'];
else
    leg_str{1,1} = DA.reg(2).specificity(s).names.all{1}; leg_str{1,2} = DA.reg(2).specificity(s).names.all{2}; leg_str{1,3} = [DA.reg(2).specificity(s).names.all{3}];
end
hl = legend(legg,leg_str);
set(hl, 'Position',[0.56, 0, 0.01, 0.07],'Orientation','horizontal','box','off','FontSize',FntsLeg,'interpreter','none');


%.-------------------------------------------------
% print z-Fisher stats
%.-------------------------------------------------
% ---- specificity -----------------
media = [];
for r = 2:length(DA.reg);
    media = [media;DA.reg(r).specificity(s).zf.s];
end
subplot(sp_r,sp_c, [(2*sp_c +1) (3*sp_c +1) (4*sp_c +1)])
hold on
for r = 1:length(media)
    bar([r],media(r),'FaceColor',colors(r,:),'EdgeColor','none');
end
set(gca,'XTick',[]);
ylabel('Specificity');
limy = get(gca,'ylim');
limx = get(gca,'xlim');
ylim([limy(1) 1]);
plot([limx(1) limx(2)],[media(1) media(1)],':k');
title ('Specificity');
box on

% ---- z-fisher T-----------------
media = [];
for r = 2:length(DA.reg);
    media = [media;DA.reg(r).specificity(s).zf.t];
end
subplot(sp_r,sp_c, [(6*sp_c +1) (7*sp_c +1)])
hold on
for r = 1:length(media)
    bar([r],media(r),'FaceColor',colors(r,:),'EdgeColor','none');
end
set(gca,'XTick',[]);
ylabel('z-fisher');
title(DA.reg(2).specificity(s).names.t,'interpreter','none');
limy = get(gca,'ylim');
limx = get(gca,'xlim');
ylim([limy(1) 1.5]);
plot([limx(1) limx(2)],[media(1) media(1)],':k');
box on

% ---- z-fisher R-----------------
media = [];
for r = 2:length(DA.reg);
    media = [media;DA.reg(r).specificity(s).zf.r];
end
subplot(sp_r,sp_c, [(9*sp_c +1) (10*sp_c +1)])
hold on
for r = 1:length(media)
    bar([r],media(r),'FaceColor',colors(r,:),'EdgeColor','none');
end
set(gca,'XTick',[]);
ylabel('z-fisher');
title(DA.reg(2).specificity(s).names.r,'interpreter','none');
limy = get(gca,'ylim');
limx = get(gca,'xlim');
ylim([limy(1) 1.5]);
plot([limx(1) limx(2)],[media(1) media(1)],':k');
box on

%.-------------------------------------------------
% print z-score stats
%.-------------------------------------------------

% ---- z-score T-----------------
media = [];
for r = 2:length(DA.reg);
    media = [media;DA.reg(r).specificity(s).zs.t];
end
subplot(sp_r,sp_c, [(6*sp_c +2) (7*sp_c +2)])
hold on
for r = 1:length(media)
    bar([r],media(r),'FaceColor',colors(r,:),'EdgeColor','none');
end
set(gca,'XTick',[]);
ylabel('z-score');
title(DA.reg(2).specificity(s).names.t,'interpreter','none');
limx = get(gca,'xlim');
plot([limx(1) limx(2)],[media(1) media(1)],':k');
box on

% ---- z-score R-----------------
media = [];
for r = 2:length(DA.reg);
    media = [media;DA.reg(r).specificity(s).zs.r];
end
subplot(sp_r,sp_c, [(9*sp_c +2) (10*sp_c +2)])
hold on
for r = 1:length(media)
    bar([r],media(r),'FaceColor',colors(r,:),'EdgeColor','none');
    string_name{r} = DA.reg(r+1).name_contract;
end
set(gca,'XTick',[]);
% set(gca,'XTick',1:length(string_name));
% set(gca,'XTickLabel',string_name);
ylabel('z-score');
title(DA.reg(2).specificity(s).names.r,'interpreter','none');
limx = get(gca,'xlim');
plot([limx(1) limx(2)],[media(1) media(1)],':k');
box on



%.-------------------------------------------------
% print FD
%.-------------------------------------------------

v = DA.rp_based;
x = 1:size(v.rp,1);
N = length(x);
subplot(sp_r,sp_c, [3 4 5 ])
plot(x,v.fd_TRnorm,'r',x,v.tra_abs,':b',x,v.rot_abs,':k','linewidth',lw);
title('FD/tr and absolute movements');
set(gca,'XTick',[]);
%set(gca,'YTick',[]);
xlim([0,N]);
ylim(opt.DA.fd_ylim);
ylabel('mm');
%--------------------------------------


%.-------------------------------------------------
% print TIMESERIES
%.-------------------------------------------------
count = 1;
for r = 2: length(DA.reg);
    y = DA.reg(r).specificity(s).tmseries;
    subplot(sp_r,sp_c, [((count)*sp_c +3) ((count)*sp_c +4) ((count)*sp_c +5) ((count+1)*sp_c +3) ((count+1)*sp_c +4) ((count+1)*sp_c +5)])
    %set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
    if size(y,2) == 4
        ref = mean([y(:,3),y(:,4)],2);
    else
        ref = y(:,3);
    end
    plot(x,y(:,1),'Color',colors_ts(1,:));hold on;
    plot(x,y(:,2),'Color',colors_ts(2,:));
    plot(x,ref,'--','Color',colors_ts(3,:));
    %(text(0,1.2, DA.reg(r).name_contract,'FontUnits',fntunit,'FontSize',reg_text.fnts,'FontWeight',reg_text.fntw,'Color',colors(r-1,:));
    set(gca,'XColor',colors(r-1,:),'YColor',colors(r-1,:),'LineWidth',2);
    set(gca,'YTick',[]);
    ylabel('%BOLDx10')
    xlim([0,N]);
    if r == length(DA.reg);
        xlabel('Volumes');        
    else
        set(gca,'XTick',[]);
    end
    count = count + 2;
end
xlabel('Volumes');
%--------------------------------------

%.-------------------------------------------------
% Save 
%.-------------------------------------------------
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

% Save the file as PNG
print([opt.folders.denoising_results_spe,'/ER_DA_',opt.folders.prject_name,'_',DA.subj_name,'_SPE_',opt.DA.roi_specificity(s).name],'-dpng','-r300');



return
end

