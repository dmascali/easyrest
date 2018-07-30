function ER_MEvsSE_plot_specificity(data,sample_names,paired,sample)

width = 7;     % Width in inches
height =3;    % Height in inches
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

figure('visible','on');

sp_r = size(data.specificity,1);
sp_c = size(data.specificity,2);
for r = 1:5
    for s = 1:size(data.specificity,2)
      isto_plot(data.specificity{r,s},sample,[sp_r sp_c (sp_c*(r-1)+s)],'t-SNR',sample_names,'',paired,fntsz,fntw);
    end
end

return
end


function isto_plot(dati,sample,subplot_xyn,title_str,sample_names,ylab,flag_paired,fntsz,fntw)
global resize_plot
z = subplot(subplot_xyn(1),subplot_xyn(2),subplot_xyn(3));
if resize_plot
    pos = get(z,'Position');
    pos(3) = pos(3) - (pos(3)./100)*20;
    set(z,'Position',pos);
end
%pbaspect([1 1.618 1]);
hold on

markers = {'o','s','d','+','*','>','<','^'};  

sample_number = length(dati);
   
media = zeros(sample_number,1);
sem = media;
for s = 1: sample_number
    media(s) = er_nanmean(dati{s});
    sem(s) = er_nanstd(dati{s})/sqrt(sum(~isnan(dati{s}))-1);
end

flag_type_plot = 1;
flag_plot_poits = 1;
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
        if flag_plot_poits
            for s = 1:sample_number
                    if sample(s).same_group_subject
                        colors = hsv(sample(s).merged_subj(1));
                        colors = repmat(colors,[sample(s).merging_number,1]);
                    else
                        colors = hsv(sample(s).tot_subj);
                    end
                    a = 1:sample(s).merging_number; % for separing group
                    a = flipdim((zscore(a)./(sample(s).merging_number+3)),2);
                    cl_indx = 0;
                    for l = 1:sample(s).merging_number
                            tmp = dati{s}([sample(s).merged_subj_index{l}]); %temporarily store data in variable "tmp"
                            x = repmat(s,1,length(tmp)); %the x axis location
                            x = (x+(rand(size(x))-0.5)*0.2) -a(l); %add a little random "jitter" to aid visibility
                            for j=1:sample(s).merged_subj(l)
                                cl_indx = cl_indx +1;
                                plot(x(j),tmp(j),markers{l},'Color',colors(cl_indx,:),'MarkerFaceColor',colors(cl_indx,:),'MarkerSize',4)
                            end
                    end

            end
        end
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