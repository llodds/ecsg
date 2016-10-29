% This script is used to generate:
%  (1) receiver data trace comparisons;
%  (2) source-receiver setup on a composite grid.

close all;
clear all;

pic_dir='/Users/muhongzhou/Documents/PhD/reports-doc/R20-ecsg-tacc/tests/test1-cmp-eci-ii/matlab-results/';

test_list = {'test1', 'test2', 'test3'};
grid_list = {'grid1', 'grid2', 'grid3'};
for test_i=1:3
    test_num = test_list{test_i};
    %% read in data
    %==========================================
    %==========================================
    % 1. read in receiver trace data from ECI folder
    cd /Users/muhongzhou/Documents/PhD/reports-doc/R20-ecsg-tacc/tests/test1-cmp-eci-ii/Fig/tests/test1-cmp-eci-ii/eci/
    rec_out = sprintf('rec_eci_%s.hh',test_num);
    dims=rsf_dim(rec_out); rec_ECI = zeros(dims'); rsf_read(rec_ECI,rec_out);
    n=dims(2); %number of traces
    % read in receiver trace geometry data
    rec_trace = sprintf('geo_eci_%s.hh',test_num);
    dims2=rsf_dim(rec_trace); rec_geo = zeros(dims2'); rsf_read(rec_geo, rec_trace);
    
%     % 2. read in receiver trace data from UNI folder
%     cd /Users/muhongzhou/Documents/PhD/reports-doc/R20-ecsg-tacc/tests/test1-cmp-eci-ii/Fig/tests/test1-cmp-eci-ii/uni/
%     rec_out = sprintf('rec_uni_%s.hh',test_num);
%     rec_UNI = zeros(dims'); rsf_read(rec_UNI,rec_out);
    
    % 3. read in receiver trace data from II folder
    cd /Users/muhongzhou/Documents/PhD/reports-doc/R20-ecsg-tacc/tests/test1-cmp-eci-ii/Fig/tests/test1-cmp-eci-ii/ii/
    rec_out = sprintf('rec_ii_%s.hh',test_num);
    rec_II = zeros(dims'); rsf_read(rec_II,rec_out);
    
    
    %% find RMS error
    %==========================================
    %==========================================
%     rms_err_ECI = zeros(n,1);
%     rms_err_II = zeros(n,1);
%     for i=1:n
%         rms_err_ECI(i) = rms(rec_ECI(:,i)-rec_UNI(:,i))/rms(rec_UNI(:,i));
%         rms_err_II(i) = rms(rec_II(:,i)-rec_UNI(:,i))/rms(rec_UNI(:,i));
%     end
%     disp(rms_err_ECI)
%     disp(rms_err_II)
    
    
    %% compare traces
    %==========================================
    %==========================================
    gcf = figure();
    line_width = 0.8;
    font_size = 9;
    
    %read in dt
    dt = rsf_par(rec_out,'d1','f',0.5);
    time=dt*((1:dims(1))-1);
    for i=1:n
        ax = subplot(n,1,i);
        hold on;
        %plot(time,rec_UNI(:,i),'k','LineWidth',line_width+0.4);
        plot(time,rec_II(:,i),'b','LineWidth',line_width);
        plot(time,rec_ECI(:,i),'r','LineWidth',line_width);
        %plot(time,abs(rec_II(:,i)-rec_UNI(:,i)),'b-.','LineWidth',line_width);
        %plot(time,abs(rec_ECI(:,i)-rec_UNI(:,i)),'r-.','LineWidth',line_width);
        
        str=sprintf('rec%d (%.0fm, %.0fm, %.0fm)', i, rec_geo(1,i), rec_geo(2,i), rec_geo(3,i));
        text(420,ax.YLim(2)*0.7,str, 'FontWeight','bold','FontSize',12);
        
        xlabel('Time [ms]', 'FontSize', font_size);
        xlim([0,dt*dims(1)]);
    end
%     l3=legend('UNI(h)','II(H=2h)','ECI(H=2h)',...
%         'abs(II-UNI)', 'abs(ECI-UNI)',...
%         'Orientation','horizontal',...
%         'Location','best');
    l3.FontSize = 12;
    l3.Position = [0.13,0.00,0.78,0.08];
    
    set(gcf, 'PaperPosition', [0 0.5 9 9.5]);
    set(gcf, 'PaperSize', [9 10]);
    pic_dir_tmp=strcat(pic_dir,test_num);
    saveas(gcf, pic_dir_tmp, 'pdf')
    
    
    %%  generate source receiver setup
    lx=1600;
    lz1=300;
    lz2=500;
    
    lz1_tmp=300-30;
    gcf2=figure();
    fill([0,lx,lx,0],[0,0,lz1_tmp,lz1_tmp],'y');
    hold on;
    fill([0,lx,lx,0],[lz1_tmp,lz1_tmp,lz1+lz2,lz1+lz2],'c');
    font_size=30;
    hold on;
    xlim([0,lx]);
    ylim([0,lz1+lz2]);
    set(gca,'Ydir','reverse')
    color1=[0.4 0.6 1.0];
    color2=[1.0 0.4 0.6];
    %   1. fine grid
    fg_step=20;
    for y=0:fg_step:lz1
        A=[0,lx];
        B=[y,y];
        plot(A,B,'Color',color1)
    end
    for x=0:fg_step:lx
        A=[x,x];
        B=[0,lz1];
        plot(A,B,'Color',color1)
    end
    %   2. coarse grid
    for y=lz1:2*fg_step:lz2+lz1
        A=[0,lx];
        B=[y,y];
        plot(A,B,'Color',color2)
    end
    for x=0:2*fg_step:lx
        A=[x,x];
        B=[lz1,lz1+lz2];
        plot(A,B,'Color',color2)
    end
    
    % plot source and receivers
    marker_size=17;
    plot(rec_geo(2,1),rec_geo(3,1),'rp','markers',marker_size,'MarkerFaceColor','r');
    plot(rec_geo(1,1),rec_geo(3,1),'b^','markers',marker_size,'MarkerFaceColor','b');
    plot(rec_geo(1,2),rec_geo(3,1),'b^','markers',marker_size,'MarkerFaceColor','b');
    plot(rec_geo(1,3),rec_geo(3,1),'b^','markers',marker_size,'MarkerFaceColor','b');
    plot(rec_geo(1,4),rec_geo(3,1),'b^','markers',marker_size,'MarkerFaceColor','b');
    
    str=sprintf(' src (x=%.0fm,y=800m,z=%.0fm)',rec_geo(2,i),rec_geo(3,i));
    
    xlabel('x');
    ylabel('z')
    set(gcf2, 'PaperPosition', [0 0 16 8]);
    set(gcf2, 'PaperSize', [16 8]);
    set(gca, 'FontSize', 25)
    title(str,'FontSize', 35);
    
    pic_dir_tmp=strcat(pic_dir,grid_list{test_i});
    saveas(gcf2, pic_dir_tmp, 'pdf')
    
end



