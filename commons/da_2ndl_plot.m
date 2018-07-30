function da_2ndl_plot(data,title_str,sample_names,output_name,paired, flag_mod,reg_names_, proj_name_)
global opt
global leg_glob_flag resize_plot
leg_glob_flag = 0;

if flag_mod == 0 %lunched from easyrest
    visible = 'off';
    reg_names = opt.DA_extracted.name_contract;
    proj_name = opt.folders.prject_name;
else    % lunched from easyrest_rois_plot
    visible = 'on';
    reg_names = reg_names_;
    proj_name = proj_name_;
end
    


width = 10;     % Width in inches
height =7.5;    % Height in inches
%alw = 0.75;    % AxesLineWidth
%fsz = 11;      % Fontsize
lw = 1;      % LineWidth
msz = 9;       % MarkerSize
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);
fntunit = 'points'; % default 'points','inches','centimeters','normalized','pixels'
text_.fnts = 10;
text_.fntw = 'bold';
FntsLeg = 7;
fntsz = 8;
fntw = 'normal';

figure('visible',visible);

sp_r = 5;
sp_c = 6;

colors = [0 0 0; ...
          204 0 0; ...
          76 153 0; ...
          0 153 153; ...
          76 0 153];
colors = colors./255;

%.-------------------------------------------------
% print all rp and intensity_based measures
%.-------------------------------------------------
subplot(sp_r,sp_c, [(1)],'Visible','off')
hold on
legg = [];
for r = 1:length(data.vTv)
    h{r} =  plot([0],[0],'Color',colors(r,:),'Linewidth',4,'visible','off');
    %leg_str{1,r} = reg_names{r};
    %legg = [legg,h{r}];
end
%hl = legend(legg,leg_str);
%set(hl,'box','off');
%set(hl, 'Position',[0.55, 0.03, 0.01, 0.05],'Orientation','horizontal','box','off','FontSize',FntsLeg);
%ax = axes('Position',[0 0 1 1],'Visible','off');
hold off

str = {'Proj. name:'
    ['  ',proj_name];
        title_str;       };
text(-0.7,1.2,str,'interpreter','none','FontUnits',fntunit,'FontSize',text_.fnts,'FontWeight',text_.fntw,'Unit','Normalize');
resize_plot = 1;
isto_plot(data.fd_TRnorm,[sp_r sp_c (sp_c +1)],'FD/tr',sample_names,'mm',paired,fntsz,fntw);
isto_plot(data.rms,[sp_r sp_c (2*sp_c +1)],'rms RP',sample_names,'mm',paired,fntsz,fntw);
isto_plot(data.rms_dvars,[sp_r sp_c (3*sp_c +1)],'rms DVARS',sample_names,'|\Delta%BOLDx10|',paired,fntsz,fntw);
isto_plot(data.rms_sd,[sp_r sp_c (4*sp_c +1)],'rms SD',sample_names,'%BOLDx10',paired,fntsz,fntw);
resize_plot = 0;
%---------------------------------------
% plot histograms
%---------------------------------------
for r = 1:length(data.vTv)
    if r == length(data.vTv)
        leg_glob_flag = 1;
    end
    plot_histograms(data.vTv{r},[sp_r sp_c], [((r-1)*sp_c +2) ((r-1)*sp_c +3)],sample_names)
     set(gca,'XColor',colors(r,:),'YColor',colors(r,:),'LineWidth',1);

     if r == 1
         title('Voxel to Voxel Correlation','FontSize',fntsz,'FontWeight',fntw);
     end
      if r == length(data.vTv)
         xlabel ('Correlation (r)','FontSize',fntsz,'FontWeight',fntw);
      else
          set(gca,'Xtick',[]);
      end
      ylabel ('Normalized Count','FontSize',fntsz,'FontWeight',fntw);
    xlim([-0.7 0.7]);
    ylim([0  0.08]);
end

%---------------------------------------
% plot mean vTv
%---------------------------------------
for r = 1:length(data.vTv)
    
    isto_plot(data.vTv_mean{r},[sp_r sp_c ((r-1)*sp_c +4)],[], sample_names,'Correlation (r)',paired,fntsz,fntw);
    set(gca,'XColor',colors(r,:),'YColor',colors(r,:),'LineWidth',1);
     if r == 1
         title('Mean correlation','FontSize',fntsz,'FontWeight',fntw);
     end
end

%---------------------------------------
% plot explained var
%---------------------------------------
for r = 1:length(data.vTv)
    
    isto_plot(data.var_exp{r},[sp_r sp_c ((r-1)*sp_c +5)],[],sample_names,[],paired,fntsz,fntw);
    set(gca,'XColor',colors(r,:),'YColor',colors(r,:),'LineWidth',1);
     if r == 1
         title('Explained Variance','FontSize',fntsz,'FontWeight',fntw);
     end
    ylim([0 100])
end

%---------------------------------------
% plot Specificity
%---------------------------------------
%here we plot just the first specificity
for r = 1:length(data.vTv)
    isto_plot(data.specificity{r,1},[sp_r sp_c ((r-1)*sp_c +6)],[],sample_names,[],paired,fntsz,fntw);
    set(gca,'XColor',colors(r,:),'YColor',colors(r,:),'LineWidth',1);
    if r == 1
         title('Specificity','FontSize',fntsz,'FontWeight',fntw);
    end
    ylim([-1 1])
end


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
print(['ER_DA_',proj_name,'_',output_name],'-dpng','-r300');



return
end

function plot_histograms(Y,subplot_xy,subplot_n,sample_names)

global leg_glob_flag

subplot(subplot_xy(1),subplot_xy(2),subplot_n)
hold on

sample_number = length(Y);
x = -0.99:0.01:0.99;
   
stats = cell(1,sample_number);
colors = hsv(sample_number);
colors_m = abs(colors - (60/255));
hs = cell(sample_number,1);
hp = hs;
leg_indx = [];
for s = 1:sample_number
    stats{s}(1,:) = mean(Y{s},1);
    stats{s}(2,:) = std(Y{s},[],1)./sqrt(size(Y{s},1)-1);
    hs{s} = shadedplot(x', (stats{s}(1,:)-stats{s}(2,:)), (stats{s}(1,:)+stats{s}(2,:)), colors(s,:));
    hp{s} = plot(x',stats{s}(1,:),'linewidth',0.8);
    set(hp{s},'Color',colors_m(s,:))
    leg_indx = [leg_indx, hp{s}];
end
if leg_glob_flag
    legend(leg_indx,sample_names,'box','off')
end

%         %s1 = shadedplot(freq, (media_REST-sd_REST), (media_REST+sd_REST), [153/255 255/255 162/255]);
%         s2 = shadedplot(freq, (media_L-sd_L), (media_L+sd_L), [255/255 102/255 51/255]);
%         s3 = shadedplot(freq, (media_R-sd_R), (media_R+sd_R), [153/255 153/255 255/255]);
%title([temp_rois(i).name],'fontsize',10,'FontWeight','bold');
%xlabel ('Correlation (r)');

%title('Voxel to Voxel Correlation');
box on
grid off
hold off


return
end

function isto_plot(dati,subplot_xyn,title_str,sample_names,ylab,flag_paired,fntsz,fntw)
global resize_plot
z = subplot(subplot_xyn(1),subplot_xyn(2),subplot_xyn(3));
if resize_plot
    pos = get(z,'Position');
    pos(3) = pos(3) - (pos(3)./100)*20;
    set(z,'Position',pos);
end
%pbaspect([1 1.618 1]);
hold on

sample_number = length(dati);
   
media = zeros(sample_number,1);
sem = media;
for s = 1: sample_number
    media(s) = er_nanmean(dati{s});
    sem(s) = er_nanstd(dati{s})/sqrt(sum(~isnan(dati{s}))-1);
end

flag_type_plot = 1;
do_stat = 1;

switch flag_type_plot
    case 2  %box plot
        temp_data = [];
        temp_g = [];
        for s = 1:sample_number
            temp_data = [temp_data;dati{s}];
            temp_g = [temp_g; s*ones(length(dati{s}),1)];
        end
        boxplot(temp_data,temp_g)
%         if flag_plot_poits
%             if flag_paired; line_x = cell(2,1); line_y = line_x; count_s = 0;end
%             for s = 1:sample_number
%                 if sample(s).same_group_subject
%                     colors = hsv(sample(s).merged_subj(1));
%                     colors = repmat(colors,[sample(s).merging_number,1]);
%                 else
%                     colors = hsv(sample(s).tot_subj);
%                 end
%                 a = 1:sample(s).merging_number; % for separing group
%                 a = flipdim((zscore(a)./(sample(s).merging_number+3)),2);
%                 if flag_paired; line_x_ = []; line_y_ = [];end
%                 cl_indx = 0;
%                 for l = 1:sample(s).merging_number
%                         tmp = dati{s}([sample(s).merged_subj_index{l}]); %temporarily store data in variable "tmp"
%                         x = repmat(s,1,length(tmp)); %the x axis location
%                         x = (x+(rand(size(x))-0.5)*0.2) -a(l); %add a little random "jitter" to aid visibility
%                         for j=1:sample(s).merged_subj(l)
%                             cl_indx = cl_indx +1;
%                             plot(x(j),tmp(j),markers{l},'Color',colors(cl_indx,:),'MarkerFaceColor',colors(cl_indx,:),'MarkerSize',5);
%                         end
%                         if flag_paired 
%                             line_x_ = [line_x_,x];
%                             line_y_ = [line_y_,tmp'];
%                         end
%                 end
%                 if flag_paired;
%                     line_x{s} = line_x_;
%                     line_y{s} = line_y_;
%                     count_s = count_s +1;
%                 end
%             end
%             if flag_paired && count_s == 2;
%                 x = [line_x{1};line_x{2}];
%                 y = [line_y{1};line_y{2}];
%                 %plot(x,y,'--k','MarkerFaceColor','none')
%                 plot(x,y,':','Color',[0.3,0.3,0.3],'MarkerFaceColor','none')
%             end
%         end
    case 1  % bar plot
        %bar_col = bone(sample_number);
        for z =1: sample_number
            %bar([z],media(z),'FaceColor',bar_col(z,:)); %'EdgeColor',[0.3,0.3,0.3]); 'EdgeColor','none');
            bar([z],media(z),'FaceColor',[0.4,0.4,0.4],'EdgeColor','none'); %'EdgeColor',[0.3,0.3,0.3]); 'EdgeColor','none');
        end
%         if flag_plot_poits
%             for s = 1:sample_number
%                     if sample(s).same_group_subject
%                         colors = hsv(sample(s).merged_subj(1));
%                         colors = repmat(colors,[sample(s).merging_number,1]);
%                     else
%                         colors = hsv(sample(s).tot_subj);
%                     end
%                     a = 1:sample(s).merging_number; % for separing group
%                     a = flipdim((zscore(a)./(sample(s).merging_number+3)),2);
%                     cl_indx = 0;
%                     for l = 1:sample(s).merging_number
%                             tmp = dati{s}([sample(s).merged_subj_index{l}]); %temporarily store data in variable "tmp"
%                             x = repmat(s,1,length(tmp)); %the x axis location
%                             x = (x+(rand(size(x))-0.5)*0.2) -a(l); %add a little random "jitter" to aid visibility
%                             for j=1:sample(s).merged_subj(l)
%                                 cl_indx = cl_indx +1;
%                                 plot(x(j),tmp(j),markers{l},'Color',colors(cl_indx,:),'MarkerFaceColor',colors(cl_indx,:),'MarkerSize',4)
%                             end
%                     end
% 
%             end
%         end
        h=errorbar(media, sem,'k');
        set(h,'LineStyle','none');
        set(h,'LineWidth',1);
    case 3 % line plot
%         if flag_plot_poits
%             for s = 1:sample_number
%                     if sample(s).same_group_subject
%                         colors = hsv(sample(s).merged_subj(1));
%                         colors = repmat(colors,[sample(s).merging_number,1]);
%                     else
%                         colors = hsv(sample(s).tot_subj);
%                     end
%                     a = 1:sample(s).merging_number; % for separing group
%                     a = flipdim((zscore(a)./(sample(s).merging_number+3)),2);
%                     cl_indx = 0;
%                     for l = 1:sample(s).merging_number
%                             tmp = dati{s}([sample(s).merged_subj_index{l}]); %temporarily store data in variable "tmp"
%                             x = repmat(s,1,length(tmp)); %the x axis location
%                             x = (x+(rand(size(x))-0.5)*0.2) -a(l); %add a little random "jitter" to aid visibility
%                             for j=1:sample(s).merged_subj(l)
%                                 cl_indx = cl_indx +1;
%                                 plot(x(j),tmp(j),markers{l},'Color',colors(cl_indx,:),'MarkerFaceColor',colors(cl_indx,:),'MarkerSize',4)
%                             end
%                     end
%             end
%         end
        h=errorbar(media, sem,'k');
        set(h,'LineWidth',1);
end

%-------------testing for significance-----------------------
if do_stat

    if sample_number == 1
        MM = 0;
        [~,p,~,stat] = ttest(dati{1},MM,0.05,'both');
        
        combs = [1];
    else
        combs = combs_no_rep(sample_number,2);
        if flag_paired == 1
            for l = 1:size(combs,1)   
                [~,p(l),~,stat(l)] = ttest(dati{combs(l,1)},dati{combs(l,2)},0.05,'both');  
            end
        else
            for l = 1:size(combs,1)   
                [~,p(l),~,stat(l)] = ttest2(dati{combs(l,1)},dati{combs(l,2)},0.05,'both');  
                
            end

        end
    end

    for l = 1:size(combs,1)
        sigstar_cell{l} = combs(l,:);
    end
    hsig = sigstar(sigstar_cell,p);
    set(hsig(:,1), 'LineWidth',0.8);
    set(hsig(:,2), 'FontSize',11);
    %set(hsig(:,2), 'FontWeight','bold');

end
%------------------------------------------------------------

%set(gca,'YTick',[0:0.1:0.4])
%set(gca,'YTickLabel',['0';' ';'0.1';' ';'0.2';' ';'0.3';' ';'0.4'])

%a = ylim;
%text(-1.1,a(2)+0.01/a(2),['\fontsize{12px}',letter]);
%text(xpos,1.070,['\fontsize{12px}',letter],'Units', 'Normalized', 'VerticalAlignment', 'Top')


set(gca,'xtick',[1:sample_number]);
set(gca,'XTickLabel',sample_names);
if not(isempty(ylab))
    ylabel(ylab,'fontsize',fntsz,'FontWeight',fntw);
end
if not(isempty(title_str))
    title(title_str,'fontsize',fntsz,'FontWeight',fntw);
end

%     set(gca,'FontSize',9);     
%set(gca,'FontWeight','bold'); 
%     set(gca, 'LineWidth'   , 0.7     );


% if plot_stim.do
%     y_lim = ylim;
%     wei_y = y_lim(2) - y_lim(1);
%     for n = 1:plot_stim.n_stim
%        for nn = 1:length(plot_stim.POS{n})
%            pos = [plot_stim.POS{n}(nn) y_lim(1) plot_stim.WEI{n}(nn) wei_y];
%            r = rectangle('Position',pos,'FaceColor', plot_stim.COL{n},'LineStyle','none');
%            uistack(r,'bottom')
%        end
%    end
% 
% end

box on
hold off
xlim([0 (sample_number+1)]);
    
return
end