close all;
clear all;

%   1. read in receiver trace data from ECI folder
cd /Users/muhongzhou/Documents/PhD/reports-doc/R20-ecsg-tacc/tests/test1-cmp-eci-ii/eci
dims=rsf_dim('rec_out.hh');
rec_ECI = zeros(dims'); 
rsf_read(rec_ECI,'rec_out.hh');
n=dims(2); %number of traces
%read in receiver trace geometry data
dims2=rsf_dim('rec_trace.hh');
rec_geo = zeros(dims2');
rsf_read(rec_geo, 'rec_trace.hh');


%   2. read in receiver trace data from UNI folder
cd /Users/muhongzhou/Documents/PhD/reports-doc/R20-ecsg-tacc/tests/test1-cmp-eci-ii/uni
rec_UNI = zeros(dims'); 
rsf_read(rec_UNI,'rec_out.hh');

%   3. read in receiver trace data from II folder
cd /Users/muhongzhou/Documents/PhD/reports-doc/R20-ecsg-tacc/tests/test1-cmp-eci-ii/ii
rec_II = zeros(dims'); 
rsf_read(rec_II,'rec_out.hh');

%%  find RMS error
rms_err_ECI = zeros(n,1);
rms_err_II = zeros(n,1);
for i=1:n
    rms_err_ECI(i) = rms(rec_ECI(:,i)-rec_UNI(:,i))/rms(rec_UNI(:,i));
    rms_err_II(i) = rms(rec_II(:,i)-rec_UNI(:,i))/rms(rec_UNI(:,i));
end
disp(rms_err_ECI)
disp(rms_err_II)
        
%%   compare traces
gcf = figure();
line_width = 0.8;
font_size = 9;

%   read in dt 
dt = rsf_par('rec_out.hh','d1','f',0.5);
time=dt*((1:dims(1))-1);
for i=1:n
    subplot(n,1,i)
    hold on;
    plot(time,rec_UNI(:,i),'k','LineWidth',line_width+0.4);
    plot(time,rec_II(:,i),'b','LineWidth',line_width);
    plot(time,rec_ECI(:,i),'r','LineWidth',line_width);
    plot(time,abs(rec_II(:,i)-rec_UNI(:,i)),'b-.','LineWidth',line_width);
    plot(time,abs(rec_ECI(:,i)-rec_UNI(:,i)),'r-.','LineWidth',line_width);
    
    str=sprintf('rec# %d (%.0fm,%.0fm,%.0fm)', i, rec_geo(1,i), rec_geo(2,i), rec_geo(3,i));
    title(str);
    
    xlabel('Time [ms]', 'FontSize', font_size);
    xlim([0,dt*dims(1)]);
end
l3=legend('UNI(h)','II(H=2h)','ECI(H=2h)',...
    'abs(II-UNI)', 'abs(ECI-UNI)',...
    'Orientation','horizontal',...
    'Location','best');
l3.FontSize = 12;
l3.Position = [0.12,0.00,0.78,0.11];

set(gcf, 'PaperPosition', [0 0 9 11]); %Position the plot further to the left and down. Extend the plot to fill entire paper.
set(gcf, 'PaperSize', [9 11]); %Keep the same paper size
saveas(gcf, '/Users/muhongzhou/Documents/PhD/reports-doc/R20-ecsg-tacc/tests/test1-cmp-eci-ii/matlab-results/test-5-098-20', 'pdf')


