function denoising_analysis_plot(DA)

global opt

width = 12;     % Width in inches
height =7;    % Height in inches
%alw = 0.75;    % AxesLineWidth
%fsz = 11;      % Fontsize
lw = 1;      % LineWidth
msz = 9;       % MarkerSize
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);
fntunit = 'points'; % default 'points','inches','centimeters','normalized','pixels'
text_.fnts = 12;
text_.fntw = 'bold';
FntsLeg = 8;

figure('visible','off');
clear h h2 h3;

colors = opt.DA.colors;

sp_r = 11;
sp_c = 5;
%.-------------------------------------------------
% print all rp and intensity_based measures
%.-------------------------------------------------

subplot(sp_r,sp_c, [(1)],'Visible','off')
%ax = axes('Position',[0 0 1 1],'Visible','off');
str = {['Proj. name: ',opt.folders.prject_name];
          ['Subj name: ',r_u(DA.subj_name)];
          ['Mean FD:              ',num2str(DA.rp_based.fd_mean),' mm'];
          ['Mean FD/tr:          ',num2str(DA.rp_based.fd_TRnorm_mean),' mm'];
          ['RMS movement:  ',num2str(DA.rp_based.rms_mean),' mm'];
          };
text(0,0.75,str,'interpreter','none','FontUnits',fntunit,'FontSize',text_.fnts,'FontWeight',text_.fntw);  %was 0,0.7

subplot(sp_r,sp_c, [(sp_c +1) (sp_c +2) (2*sp_c +1) (2*sp_c +2)])
%set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
v = DA.rp_based;
x = 1:size(v.rp,1);
N = length(x);
plot(x,v.tra(:,1),'b',x,v.tra(:,2),'g',x,v.tra(:,3),'r','linewidth',lw);
ylabel('mm');
hl = legend('x', 'y', 'z');
%set(hl, 'Position',[0.065, 0.75, 0.01, 0.07],'box','off');
set(hl, 'Orientation','horizontal','box','off','FontSize',FntsLeg,'Location','best');
set(gca,'XTick',[]);
xlim([0,N]);

subplot(sp_r,sp_c, [(3*sp_c +1) (3*sp_c +2) (4*sp_c +1) (4*sp_c +2)])
%set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
v = DA.rp_based;
x = 1:size(v.rp,1);
plot(x,v.rot(:,1),'b',x,v.rot(:,2),'g',x,v.rot(:,3),'r','linewidth',lw);
ylabel(opt.AUX.rp_rot_unit);
hl = legend('pitch','roll','yaw');
%set(hl, 'Position',[0.065, 0.60, 0.01, 0.07],'box','off');
set(hl, 'Orientation','horizontal','box','off','FontSize',FntsLeg,'Location','best');
set(gca,'XTick',[]);
xlim([0,N]);


subplot(sp_r,sp_c, [(5*sp_c +1) (5*sp_c +2) (6*sp_c +1) (6*sp_c +2)])
plot(x,v.fd_TRnorm,'r',x,v.tra_abs,':b',x,v.rot_abs,':k','linewidth',lw);
hl = legend('FD/tr', '|tra|', '|rot|');
%set(hl, 'Position',[0.065, 0.45, 0.01, 0.07],'box','off');
set(hl, 'Orientation','horizontal','box','off','FontSize',FntsLeg,'Location','best');
ylabel('mm');
set(gca,'XTick',[]);
xlim([0,N]);

subplot(sp_r,sp_c, [(7*sp_c +1) (7*sp_c +2) (8*sp_c +1) (8*sp_c +2)])
V = DA.intensity_based;
plot(x,V.DVARS,'linewidth',lw);
hl = legend('DVARS');
%set(hl, 'Position',[0.06, 0.42, 0.01, 0.07],'box','off');
set(hl,'box','off');
%set(gca,'YTick',[]);
ylabel('|\Delta%BOLDx10|');
set(gca,'XTick',[]);
xlim([0,N]);

subplot(sp_r,sp_c, [(9*sp_c +1) (9*sp_c +2) (10*sp_c +1) (10*sp_c +2)])
V = DA.intensity_based;
plot(x,V.SD,'g',x,V.GS,'k','linewidth',lw);
hl = legend('SD','GS','Orientation','horizontal','Location','best');
%set(hl, 'Position',[0.06, 0.22, 0.01, 0.07],'box','off');
set(hl,'box','off');
%set(gca,'YTick',[]);
ylabel('%BOLDx10');
xlabel('Volumes');
xlim([0,N]);




%.-------------------------------------------------
% print regressions
%.-------------------------------------------------

subplot(sp_r,sp_c, [3 4 ])

plot(x,v.fd_TRnorm,'r',x,v.tra_abs,':b',x,v.rot_abs,':k','linewidth',lw);
title('FD/tr and absolute movements');
set(gca,'XTick',[]);
%set(gca,'YTick',[]);
xlim([0,N]);
ylim(opt.DA.fd_ylim);
ylabel('mm');

count = 1;
for r = 2: length(DA.reg);
    
    subplot(sp_r,sp_c, [((count)*sp_c +3) ((count)*sp_c +4) ((count+1)*sp_c +3) ((count+1)*sp_c +4)])
    %set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
    Y = zscore(DA.reg(r).data)';
    imagesc(Y,[-1 1]); colormap gray;
    %(text(0,1.2, DA.reg(r).name_contract,'FontUnits',fntunit,'FontSize',reg_text.fnts,'FontWeight',reg_text.fntw,'Color',colors(r-1,:));
    set(gca,'XColor',colors(r-1,:),'YColor',colors(r-1,:),'LineWidth',2);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    ylabel('GM Voxels')
    count = count + 2;
end
xlabel('Volumes');




%.-------------------------------------------------
% print correlations
%.-------------------------------------------------

subplot(sp_r,sp_c, [5 (sp_c +5) (2*sp_c +5) ])


hold on
X = DA.reg(2).vTv_sampled.x;
h{1} = plot(X,DA.reg(2).vTv_sampled.y,'Color',colors(1,:),'LineWidth', lw);
leg_str{1,1} = DA.reg(2).name_contract;

legg = h{1};
for r = 3:length(DA.reg);
    h{r-1} =  plot(X,DA.reg(r).vTv_sampled.y,'Color',colors(r-1,:),'LineWidth', lw);
    leg_str{1,r-1} = DA.reg(r).name_contract;
    legg = [legg,h{r-1}];
end
hl = legend(legg,leg_str);
set(hl, 'Position',[0.66, 0.03, 0.01, 0.05],'Orientation','horizontal','box','off','FontSize',FntsLeg);
lim = get(gca,'ylim');
plot([0 0],[0 lim(2)],':k');
xlim([-0.7 0.7]);
set(gca,'YTick',[]);
xlabel ('Correlation (r)');
ylabel ('Normalized Count');
title('Voxel to Voxel correlation');
box on


% %.-------------------------------------------------
% % print correlations stats vTv
% %.-------------------------------------------------
% 
% media = [];
% stdv = [];
% for r = 2:length(DA.reg);
%     media = [media;DA.reg(r).vTv_sampled.mean];
%     stdv = [stdv;DA.reg(r).vTv_sampled.std];
% end
% subplot(sp_r,sp_c, [(4*sp_c +5) (5*sp_c +5)])
% h2=errorbar(media, stdv,'k'); hold on
% for r = 1:length(media)
%     plot([r],media(r),'x','MarkerEdgeColor',colors(r,:),'MarkerSize',msz);
% end
% lim = get(gca,'xlim');
% plot([lim(1) lim(2)],[0 0],':k');
% title('Mean correlation');
% ylabel('Correlation (r)');
% set(gca,'XTick',[]);
% set(h2,'LineStyle','none');
% set(h2,'LineWidth',1);

%.-------------------------------------------------
% print tSNR
%.-------------------------------------------------

media = [];
stdv = [];
for r = 2:length(DA.reg);
    media = [media;DA.reg(r).tSNR.mean];
    stdv = [stdv;DA.reg(r).tSNR.std];
end
subplot(sp_r,sp_c, [(4*sp_c +5) (5*sp_c +5)]); hold on
for r = 1:length(media);
    bar([r],media(r),'FaceColor',colors(r,:),'EdgeColor','none');
end
h3=errorbar(media, stdv,'k'); 
set(h3,'LineStyle','none');
set(h3,'LineWidth',1);
set(gca,'XTick',[]);
ylabel('tSNR');
limx = get(gca,'xlim');
plot([limx(1) limx(2)],[media(1) media(1)],':k');
box on; grid on;




%.-------------------------------------------------
% print variance explained.
%.-------------------------------------------------
media = [];
stdv = [];
for r = 2:length(DA.reg);
    media = [media;DA.reg(r).explained_var.mean];
    stdv = [stdv;DA.reg(r).explained_var.std];
end
subplot(sp_r,sp_c, [(6*sp_c +5) (7*sp_c +5)]); hold on
for r = 1:length(media);
    bar([r],media(r),'FaceColor',colors(r,:),'EdgeColor','none');
end
h3=errorbar(media, stdv,'k'); 
set(h3,'LineStyle','none');
set(h3,'LineWidth',1);
set(gca,'XTick',[]);
ylabel('Explained Variance');
ylim([0 100]);
limx = get(gca,'xlim');
plot([limx(1) limx(2)],[media(1) media(1)],':k');
box on; grid on;

%.-------------------------------------------------
% print specificity.
%.-------------------------------------------------
% in this plot let's consider only the first specificity (DMN)
media = [];
for r = 2:length(DA.reg);
    media = [media;DA.reg(r).specificity(1).zf.s];   
end
subplot(sp_r,sp_c, [(9*sp_c +5) (10*sp_c +5)])
hold on
for r = 1:length(media)
    bar([r],media(r),'FaceColor',colors(r,:),'EdgeColor','none');
end
limy = get(gca,'ylim');
limx = get(gca,'xlim');
ylim([limy(1) 1]);
plot([limx(1) limx(2)],[media(1) media(1)],':k');
title('FC PCC:MPFC vs PCC:VIS');
ylabel('Specificity');
set(gca,'XTick',[]);
box on


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
print([opt.folders.denoising_results,'/ER_DA_',opt.folders.prject_name,'_',DA.subj_name],'-dpng','-r300');


%remember to restore defpos

return
end
