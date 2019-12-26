% zerogseq_04_vibration_analysis.m
% Analyze high frequency vibration data
%
% Should be run twice, first with:
% dataset = 'Ground' and
% dataset = 'Flight'
%
% Edit code to adjust this setting.
%

%% Start fresh
% Uncomment below two lines if using this as a standalone script
% clear all; close all; clc;
% addpath('./code');

% Uncomment one of the following options (must use capital letters)
% if using this as a standalone script.
% dataset = 'Ground';
% dataset = 'Flight';

% Graphics format for saving figures {file_extension, MATLAB print option}
% Type 'help print' for examples
figformat = {'eps' '-depsc'};
%figformat = {'pdf' '-dpdf'};

% LineWidth for figure plotting
lw = 0.5;
fontsize = 7;
fontname = 'Helvetica';

% Folder for outputs (e.g., figures)
outfolder = ['./analysis/vibration/' dataset];

%% Create output folder
if ~exist(outfolder,'dir'), mkdir(outfolder); end

% Load High Frequency vibration data in channel 08
load(sprintf('./data/SlamStick/%s/Ch08.mat',dataset));

% SlamStick MAT files store time as POSIX
% To convert, use following syntax

% Get the start time
t_start_local = datetime(ADC(1,1),'ConvertFrom','posixtime','TimeZone','America/New_York');
t_start_utc = datetime(ADC(1,1),'ConvertFrom','posixtime','TimeZone','UTC');

% Convert to elapsed time
t_elapsed_s = ADC(1,:)-ADC(1,1);

% Get sampling frequency and period
Ts = mean(diff(t_elapsed_s));
Fs = 1/Ts;

% save variables (data is stored as g-level equivalent vibration)
g_x = ADC(2,:);
g_y = ADC(3,:);
g_z = ADC(4,:);

% Compute equivalent g-level scalar
g = sqrt(g_x.^2 + g_y.^2 + g_z.^2);

% Compute PSD using Welch's method (pwelch.m) with default parameters.
% Units: g^2 / Hz. "By default, X is divided into the longest possible 
% sections, to get as close to but not exceeding 8 segments with 50% 
% overlap."
[psd,f]=pwelch(g,[],[],[],Fs);

%% Generate the PSD figure

% Filename for this figure
fn = 'zerogseq_04_fig_vibration_psd';

% Typical plot limits for PSD for aircraft vibrations
ylim = [0 1E-3];
xlim = [0 Fs/2];

% Make the figure
figure; set(gcf,'color',[1 1 1]);

% compute cumulative sum of power spectrum with frequency
cpsd = cumsum(psd)/sum(psd);

% First plot
subplot(2,1,1);
plot(f,psd,'linewidth',lw);
set(gca,'ylim',ylim,'xlim',xlim);
xlabel('Frequency (Hz)');
ylabel('vibration PSD (g^2 / Hz)');

% Plot cumulative sum of power spectrum with frequency
subplot(2,1,2); 
plot(f,cpsd,'linewidth',lw); hold on;
set(gca,'xscale','log');
xlabel('Frequency (Hz)'); ylabel('PSD Normalized Cum. Sum');
set(gca,'xlim',xlim);

% Print
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

% Show inset plot of spectrum power near DC

fig=figure; set(gcf,'color',[1 1 1]); 
plot(f,20*log10(psd),'linewidth',lw); 
set(gca,'xlim',[0 5],'ylim',[-200 100]); 
xlabel('Frequency (Hz)'); ylabel('PSD (dB)');
fig.PaperUnits = 'inches'; w = 3; h = 2;
fig.PaperPosition = [(8.5-w)/2 (11-h)/2 w h];

% Print
print(fullfile(outfolder,[fn '.inset.' figformat{1}]),figformat{2});

%% Design filter
d = designfilt('highpassiir', 'StopbandFrequency', 5, 'PassbandFrequency', 10, 'StopbandAttenuation', 60, 'PassbandRipple', 1, 'SampleRate', Fs);

% show frequency response of filter
figure('color',[1 1 1]); set(gca,'fontname',fontname,'fontsize',fontsize);
[mag,frq] = freqz(d,4096,Fs);
plot(frq,abs(mag),'-b','LineWidth',lw);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
set(gca,'xlim',[0 20]);

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 3 2]);

fn = 'zerogseq_04_fig_rms_vibration_filter';
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

%% Calculate the filtered "G-level" (this is really vibration)
% Filtering g_level preserves frequency content; calculating g-level from
% filtered g_x, g_y, g_z, does not preserve high frequency peaks
g_f = filtfilt(d,g);

% Compute PSD using Welch's method (pwelch.m) with default parameters.
% Units: g^2 / Hz. "By default, X is divided into the longest possible 
% sections, to get as close to but not exceeding 8 segments with 50% 
% overlap."
[psd_f,f_f]=pwelch(g_f,[],[],[],Fs);

%% Generate the filtered PSD figure

% Filename for this figure
fn = 'zerogseq_04_fig_vibration_filt_psd';
% Typical plot limits for PSD for aircraft vibrations
ylim = [0 1E-3];
xlim = [0 Fs/2];

% Make the figure
figure; set(gcf,'color',[1 1 1]);

% compute cumulative sum of power spectrum with frequency
cpsd_f = cumsum(psd_f)/sum(psd_f);

% First plot
subplot(2,1,1);
plot(f_f,psd_f,'linewidth',lw);
set(gca,'ylim',ylim,'xlim',xlim);
xlabel('Frequency (Hz)');
ylabel('vibration PSD (g^2 / Hz)');

% Plot cumulative sum of power spectrum with frequency
subplot(2,1,2); 
plot(f_f,cpsd_f,'linewidth',lw); hold on;
set(gca,'xscale','log');
xlabel('Frequency (Hz)'); ylabel('PSD Normalized Cum. Sum');
set(gca,'xlim',xlim);

% Print
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

% Show inset plot of spectrum power near DC

fig=figure; set(gcf,'color',[1 1 1]); 
plot(f_f,20*log10(psd_f),'linewidth',lw); 
set(gca,'xlim',[0 5],'ylim',[-200 100]); 
xlabel('Frequency (Hz)'); ylabel('PSD (dB)');
fig.PaperUnits = 'inches'; w = 3; h = 2;
fig.PaperPosition = [(8.5-w)/2 (11-h)/2 w h];

% Print
print(fullfile(outfolder,[fn '.inset.' figformat{1}]),figformat{2});

% Print -- for main paper
fig = figure('color',[1 1 1]);
fig.PaperUnits = 'inches'; w = 3; h = 4;
fig.PaperPosition = [(8.5-w)/2 (11-h)/2 w h];
windowsize = floor(numel(f)/max(f)); % 1 Hz window
semilogy(f,smoothdata(psd,'movmean',windowsize),'linewidth',lw*2,'Color',[1 0 0]); hold on;
semilogy(f_f,smoothdata(psd_f,'movmean',windowsize),'linewidth',lw*2,'Color',[0 0 1]);
set(gca,'ylim',[1E-8 1E-3],'xlim',xlim);
xlabel('Frequency (Hz)');
ylabel('Vibration PSD (g^2 / Hz)');
print(fullfile(outfolder,[fn '.1kHz.' figformat{1}]),figformat{2});

% Print -- for main paper
xlim = [0.1 Fs/2];
fig = figure('color',[1 1 1]);
fig.PaperUnits = 'inches'; w = 7; h = 4;
fig.PaperPosition = [(8.5-w)/2 (11-h)/2 w h];
windowsize = floor(numel(f)/max(f)); % 1 Hz window
loglog(f,smoothdata(psd,'movmean',windowsize),'linewidth',lw*2,'Color',[1 0 0]); hold on;
loglog(f_f,smoothdata(psd_f,'movmean',windowsize),'linewidth',lw*2,'Color',[0 0 1]);
set(gca,'ylim',[1E-8 1E-3],'xlim',xlim);
xlabel('Frequency (Hz)');
ylabel('Vibration PSD (g^2 / Hz)');
print(fullfile(outfolder,[fn '.1kHz.loglog.' figformat{1}]),figformat{2});

%% calculate RMS as function of time
N_s = ceil(max(t_elapsed_s));
r = NaN(1,N_s);
for k=1:N_s
    bT = bitand(t_elapsed_s>=(k-1),t_elapsed_s<k);
    r(1,k) = rms(g_f(bT));
end

%% Plot vibration data, filtered vibration, and RMS vibration (1 s bin)

h = figure('color',[1 1 1]);
subplot(3,1,1);
if strcmp(dataset,'Flight')
    plot(t_elapsed_s,g_x-mean(g_x)+3); hold on;
    plot(t_elapsed_s,g_y-mean(g_y)+1.5); hold on;
    plot(t_elapsed_s,g_z-mean(g_z)+0); hold on;
    % Add custom 1g scale bar
    plot(6500*[1 1],[1 2],'-k','linewidth',1);
    set(gca,'ylim',[-2.5 5]);
    set(gca,'xlim',[0 6750]);
else
    plot(t_elapsed_s,g_x-mean(g_x)+20); hold on;
    plot(t_elapsed_s,g_y-mean(g_y)+10); hold on;
    plot(t_elapsed_s,g_z-mean(g_z)+0); hold on;
    plot(2400*[1 1],[10 20],'-k','linewidth',1);
end
xlabel('Elapsed time (s)');
ylabel('Vibration (g)');
set(gca,'ytick',[]);
set(gca,'yticklabel',{});

subplot(3,1,2);
plot(t_elapsed_s,g,'Linewidth',lw); hold on;
plot(t_elapsed_s,g_f,'Linewidth',lw); hold on;
xlabel('Elapsed time (s)');
ylabel('Vibration (g)');
if strcmp(dataset,'Flight') set(gca,'xlim',[0 6750]); end;

subplot(3,1,3);
plot(1:N_s,r(1,:)); hold on;
xlabel('Elapsed time (s)');
ylabel('RMS Vibration (g)');
if strcmp(dataset,'Flight') set(gca,'xlim',[0 6750]); end;

fn = 'zerogseq_04_fig_overview';
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

% add legends
subplot(3,1,1);
legend('X','Y','Z','scale');
subplot(3,1,2);
legend('Unfiltered','Filtered');

fn = 'zerogseq_04_fig_overview_with_legends';
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});


%% Plot the RMS vibration (g-level)
h = figure('color',[1 1 1]);
histogram(r,1000);
xlabel('RMS g_filt');
ylabel('Frequency (1s intervals)');

fn = 'zerogseq_04_fig_rms_vibration_histogram';
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

% Save RMS data
save(fullfile(outfolder,'zerogseq_04_rms_vibration.mat'),'r');
