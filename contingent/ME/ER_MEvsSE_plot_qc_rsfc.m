function ER_MEvsSE_plot_qc_rsfc(data,sample_names)

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

sp_r = 5;
sp_c = length(sample_names);

resolution = 6; %mm
l = min(data.d_vect):resolution:max(data.d_vect);

for r = 1:5
    for s = 1:length(sample_names)
        
        Y = data.roitoroi.zs.vect{r}{s}';
        Y = zscore(Y);
        X = zscore(data.fd{s});
        c = X\Y;
        subplot(sp_r,sp_c, [((r-1)*sp_c+ s)]);
        scatter(data.d_vect,c,'Marker','.','MarkerEdgeColor','k'); hold on;
        
        
        
        %test random
        Xr = X(randperm(length(X)));
        cr = Xr\Y;
        scatter(data.d_vect,cr,'Marker','.','MarkerEdgeColor','g');
        
        plot([0 180],[0 0],'-k','Linewidth',1.5);
        
        zf = atanh(c);
        zf_r = atanh(cr);
        % average
        C_avg = [];
        for z = 2:(length(l))     
            C_avg(z-1) = nanmean(zf(data.d_vect<= l(z) & data.d_vect > l(z-1)));
        end
        
        C_avg_r = [];
        for z = 2:(length(l))     
            C_avg_r(z-1) = nanmean(zf_r(data.d_vect<= l(z) & data.d_vect > l(z-1)));
        end
        C_avg_d = l(1:end-1) + diff(l)/2;
        plot(C_avg_d,C_avg,'-r','Linewidth',2.5);
        plot(C_avg_d,C_avg_r,'-w','Linewidth',2);
        
        %plot(DA.reg(r).roitoroi.d_vect_avg,y2,'-r','Linewidth',4);
        %set(gca,'XColor',colors(r-1,:),'YColor',colors(r-1,:),'LineWidth',2);
        %uistack(hs,'bottom');
        %if r == length(DA.reg)
            xlabel('mm');
        %else        
         %   set(gca,'XTick',[]);
        %end
        %if r == 2
            title('ROI-to-ROI FC vs distance');
        %end
        %set(gca,'YTick',[]);
        %if unit_r
%             ylabel('Pearson r');
%         else
            ylabel('z-score');
        %end
        %ylim([matrix_r_lim(1) matrix_r_lim(2)]);
       
        xlim([0 180]);
        box on;
        
        % z-score transf
    end
end
    
return
end