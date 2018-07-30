function denoising_analysis_plot_extra(DA,unit_r)

global opt

width = 12;     % Width in inches
height =8;    % Height in inches
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

if unit_r %perason otherwise zscore
    matrix_r_lim = [-1,1];
    plot_r_lim = [-1 1];
    flag_r = 'r';   %for output name
else
    matrix_r_lim = [-4,4];
    plot_r_lim = [-9 9];
    flag_r = 'zs';   %for output name
end


sp_r = 5;
sp_c = 6;
%.-------------------------------------------------
% print all rp and intensity_based measures
%.-------------------------------------------------

subplot(sp_r,sp_c, [(5)],'Visible','off')
%ax = axes('Position',[0 0 1 1],'Visible','off');
str = {['Proj. name: ',opt.folders.prject_name];
          ['Subj name: ',r_u(DA.subj_name)];
          ['Mean FD:              ',num2str(DA.rp_based.fd_mean),' mm'];
          ['Mean FD/tr:          ',num2str(DA.rp_based.fd_TRnorm_mean),' mm'];
          ['RMS movement:  ',num2str(DA.rp_based.rms_mean),' mm'];
          };
text(0,0.75,str,'interpreter','none','FontUnits',fntunit,'FontSize',text_.fnts,'FontWeight',text_.fntw);  %was 0,0.7




%.-------------------------------------------------
% print roi to roi matrix
%.-------------------------------------------------
count = 0;
for r = 2: length(DA.reg);
    if unit_r
        y = DA.reg(r).roitoroi.r.matrix;
    else
        y = DA.reg(r).roitoroi.zs.matrix;
    end
    subplot(sp_r,sp_c, ((count)*sp_c +1));
    imagesc(y,[matrix_r_lim(1) matrix_r_lim(2)]);
	
    pos = get(gca,'pos');
    cb = colorbar('Westoutside');
    set(gca,'pos',pos);
    if unit_r
        ylabel(cb,'Pearson r');
    else
        ylabel(cb,'z-score');
    end
        
	
    pbaspect([1 1 1]);
    set(gca,'XColor',colors(r-1,:),'YColor',colors(r-1,:),'LineWidth',2);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    count = count +1;
    box on
end

%.-------------------------------------------------
% print roi to roi r vs distance
%.-------------------------------------------------
count = 0;
for r = 2: length(DA.reg);
    if unit_r
        y1 = DA.reg(r).roitoroi.r.vect;
        y2 = DA.reg(r).roitoroi.r.vect_avg;
    else
        y1 = DA.reg(r).roitoroi.zs.vect;
        y2 = DA.reg(r).roitoroi.zs.vect_avg;
    end
    subplot(sp_r,sp_c, [((count)*sp_c +2) ((count)*sp_c +3)]);
    hs = scatter(DA.reg(r).roitoroi.d_vect,y1,'Marker','.','MarkerEdgeColor','k'); hold on;
    plot([0 180],[0 0],'-b');
    plot(DA.reg(r).roitoroi.d_vect_avg,y2,'-r','Linewidth',4);
    set(gca,'XColor',colors(r-1,:),'YColor',colors(r-1,:),'LineWidth',2);
    uistack(hs,'bottom');
    if r == length(DA.reg)
        xlabel('mm');
    else        
        set(gca,'XTick',[]);
    end
    if r == 2
        title('ROI-to-ROI FC vs distance');
    end
    %set(gca,'YTick',[]);
    if unit_r
        ylabel('Pearson r');
    else
        ylabel('z-score');
    end
    ylim([plot_r_lim(1) plot_r_lim(2)]);
    count = count +1;
    xlim([0 180]);
    box on;
    
end

%.-------------------------------------------------
% print roi to roi matrix DELTA R
%.-------------------------------------------------
count = 1;
for r = 3: length(DA.reg);
    if unit_r
        y = DA.reg(r).roitoroi.r.delta_matrix;
    else
        y = DA.reg(r).roitoroi.zs.delta_matrix;
    end
    subplot(sp_r,sp_c, ((count)*sp_c +4));
    imagesc(y,[matrix_r_lim(1) matrix_r_lim(2)]);
    if r == 3
        pos = get(gca,'pos');
        cb = colorbar('northoutside');
        set(gca,'pos',pos);
         if unit_r
            title(cb,'Pearson r');
         else
            title(cb,'z-score');
         end
    end
    pbaspect([1 1 1]);
    set(gca,'XColor',colors(r-1,:),'YColor',colors(r-1,:),'LineWidth',2);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    count = count +1;
    box on
end

%.-------------------------------------------------
% print roi to roi DELTAr vs distance
%.-------------------------------------------------
count = 1;
for r = 3: length(DA.reg);
    if unit_r
        y1 = DA.reg(r).roitoroi.r.delta_vect;
        y2 = DA.reg(r).roitoroi.r.delta_vect_avg;
    else
        y1 = DA.reg(r).roitoroi.zs.delta_vect;
        y2 = DA.reg(r).roitoroi.zs.delta_vect_avg;
    end
    subplot(sp_r,sp_c, [((count)*sp_c +5) ((count)*sp_c +6)]);
    scatter(DA.reg(r).roitoroi.d_vect,y1,'Marker','.','MarkerEdgeColor','k'); hold on;
    plot([0 180],[0 0],'-b');
    plot(DA.reg(r).roitoroi.d_vect_avg,y2,'-w','Linewidth',2);
    set(gca,'XColor',colors(r-1,:),'YColor',colors(r-1,:),'LineWidth',2);
    if r == length(DA.reg)
        xlabel('mm');
    else        
        set(gca,'XTick',[]);
    end
    if r == 3
        title('Delta ROI-to-ROI FC vs distance');
    end
    %set(gca,'YTick',[]);
    if unit_r
        ylabel('Pearson r');
    else
        ylabel('z-score');
    end
    ylim([plot_r_lim(1) plot_r_lim(2)]);
    xlim([0 180]);
    count = count +1;
    box on
end

% %.-------------------------------------------------
% % print BOLD:FD correlation (as violin plot)
% %.-------------------------------------------------
% y= [];
% for r = 2:length(DA.reg);
%     if unit_r
%         y = [y,DA.reg(r).bold_motion.fd.corr'];
%     else
%         y = [y,DA.reg(r).bold_motion.fd.zf'];
%     end
% end
% subplot(sp_r,sp_c, [(sp_c +1) (sp_c +2)]);
% plot([0 (size(y,2)+1)],[0 0],'k'); hold on
% violin(y,'facecolor',colors,'edgecolor','none','xlabel',[],'medc',[],'facealpha',0.7);
% set(gca,'XTick',[]);
% title('BOLD:FD correlation');
% box on
% 
% %.-------------------------------------------------
% % print BOLD:FD correlation (as violin plot)
% %.-------------------------------------------------
% y= [];
% for r = 2: length(DA.reg);
%     if unit_r
%         y = [y,DA.reg(r).bold_motion.dvars.corr'];
%     else
%         y = [y,DA.reg(r).bold_motion.dvars.zf'];
%     end
% end
% subplot(sp_r,sp_c, [(sp_c*2 +1) (sp_c*2 +2)]);
% plot([0 (size(y,2)+1)],[0 0],'k'); hold on
% violin(y,'facecolor',colors,'edgecolor','none','xlabel',[],'medc',[],'facealpha',0.7);
% set(gca,'XTick',[]);
% title('BOLD:DVARS correlation');
% box on


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
print([opt.folders.denoising_results_rtr,'/ER_DA_',opt.folders.prject_name,'_',DA.subj_name,'_extra_',flag_r],'-dpng','-r300');


%remember to restore defpos

return
end
