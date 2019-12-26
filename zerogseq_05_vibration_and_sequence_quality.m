% zerogseq_05_vibration_and_sequence_quality.m
%
% Analyze the relationship between RMS vibration, G-level, and sequence phred
% quality score via plotting and regression.
%
% Run script twice, once with dataset = 'Flight' and once with dataset =
% 'Ground'. Edit script to set this option.
%

% Start fresh
% Uncomment below two lines if using this as a standalone script
% clear all; close all; clc;
% addpath('./code');

% Uncomment one of the following options (must use capital letters)
% if using this as a standalone script.
% dataset = 'Ground';
% dataset = 'Flight';

%% Genome length
genome_length = 48502;

%% Choose correct options

% FLIGHT 
if strcmpi(dataset,'Flight')
    vib_data = './analysis/vibration/Flight/zerogseq_04_rms_vibration.mat';
    g_data = './analysis/acceleration/G_filt.mat';
    read_data = './analysis/MinION/Flight/flight_run/Reads.mat';
    mux_read_data = './analysis/MinION/Flight/flight_mux/Reads.mat';
    base_data = './analysis/MinION/Flight/flight_run/Bases.mat';
    run_offset_s = 1061;
    mux_run_offset_s = 329;
    outfolder = './analysis/Combined/Flight';
    % Regression restricted to just before and just after parabolas to avoid 
    % complicating factors of g-level and vibration changes during aircraft 
    % descent
    t_max = 4000; 
    logfile = fullfile(outfolder,'zerogseq_05_vibration_sequence_quality_regression_Flight.txt');
    periods_file = './analysis/acceleration/periods.txt';
    xlim = [0 7000];
end

% GROUND
if strcmpi(dataset,'Ground')
    vib_data = './analysis/vibration/Ground/zerogseq_04_rms_vibration.mat';
    g_data = '';
    read_data = './analysis/MinION/Ground/ground_run/Reads.mat';
    mux_read_data = './analysis/MinION/Ground/ground_mux/Reads.mat';
    run_offset_s = 406;
    mux_run_offset_s = -23;
    outfolder = './analysis/Combined/Ground';
    % Regression across all data
    t_max = 2266; 
    logfile = fullfile(outfolder,'zerogseq_05_vibration_sequence_quality_regression_Ground.txt');
    periods_file = '';
    xlim = [0 2500];
end

% Check if output folder exists
if ~(exist(outfolder,'dir')==7), mkdir(outfolder); end

% Start logging
diary(logfile);

% Graphics format for saving figures {file_extension, MATLAB print option}
% Type 'help print' for examples
figformat = {'eps' '-depsc' '-painters'};
%figformat = {'pdf' '-dpdf'};

% LineWidth for figure plotting
lw = 0.5;
fontsize = 7;
fontname = 'Helvetica';

% Load vibration data
rms = load(vib_data);
rms = rms.r;
% Load read info
reads = load(read_data);
reads = reads.reads;
% load mux read info
muxreads = load(mux_read_data);
muxreads = muxreads.reads;

% Caluclate Tombo Read times in accelerometer elapsed time
fs = reads.sampling_rate_Hz;
t_start_elapsed_read_s = reads.start_time_samples./fs;
t_tombo_offset_s = reads.tombo_alignment_offset_samples./fs;
t_tombo_duration_s = reads.tombo_duration_samples./fs;
t_start_s = run_offset_s + t_start_elapsed_read_s + t_tombo_offset_s;
t_stop_s = t_start_s + t_tombo_duration_s;
% Now for mux
mux_fs = muxreads.sampling_rate_Hz;
mux_t_start_elapsed_read_s = muxreads.start_time_samples./mux_fs;
mux_t_tombo_offset_s = muxreads.tombo_alignment_offset_samples./mux_fs;
mux_t_tombo_duration_s = muxreads.tombo_duration_samples./mux_fs;
mux_t_start_s = mux_run_offset_s + mux_t_start_elapsed_read_s + mux_t_tombo_offset_s;
mux_t_stop_s = mux_t_start_s + mux_t_tombo_duration_s;

% Plot
figure('color',[1 1 1]);
set(gca,'fontsize',fontsize,'fontname',fontname);
plot(rms,'LineWidth',lw); hold on;
xlabel('Elapsed Time(s)');
ylabel('RMS vibration (g), 1s bin');
yyaxis right
for k=1:numel(t_start_s)
    % plot this read
    plot([t_start_s(k) t_stop_s(k)],reads.basecalled_quality_phred(k)*[1 1],'-r','LineWidth',lw);
end
ylabel('Phred Quality Score');

% Add mux reads to plot
yyaxis right
for k=1:numel(mux_t_start_s)
    % plot this read
    plot([mux_t_start_s(k) mux_t_stop_s(k)],muxreads.basecalled_quality_phred(k)*[1 1],'-k','LineWidth',lw);
end
set(gca,'xlim',xlim);

% Save figure
fn = 'zerogseq_05_fig_rms_vibration_and_sequence_quality';
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

%% Figure: Mars Lunar/Europa Focus

if strcmpi(dataset,'Flight')
    % Read in periods data
    periods = readtable(periods_file);
    % Read in filtered g-level data
    gdata = load('./analysis/acceleration/G_filt.mat');
    
    % Periods for first five parabolas are 5, 9, 13, 17, 21
    % Plot region period 4:22
    p_region = [2:24];
    t_region = [periods.time_start(p_region(1)) periods.time_stop(p_region(end))];
    
    % Identify reads in region
    bReadEndsInRegion = bitand(t_stop_s>=t_region(1),t_stop_s<=t_region(end));
    bReadStartsInRegion = bitand(t_start_s>=t_region(1),t_start_s<=t_region(end));
    bReadInRegion = bitor(bReadEndsInRegion,bReadStartsInRegion);
    
    id = find(bReadInRegion);
    
    % Identify g-level data in region
    bGlevel_region = bitand(gdata.t>=t_region(1),gdata.t<=t_region(2));
    t_gl = gdata.t(bGlevel_region);
    g_gl = gdata.g_filt(bGlevel_region);
    
    % Plot
    fig = figure('color',[1 1 1]);
    set(gca,'fontsize',fontsize,'fontname',fontname);
    % options
    ylim = [0 1.8];
    xlim = [min(t_gl) max(t_gl)];
    % Plot transition points
    tc = [0.9 0.9 0.9];
    for k=2:2:(numel(p_region)-2)
        % line at end of region
        t_k = periods.time_stop(p_region(k));
        t_k1 = periods.time_stop(p_region(k+1));
        xd = [t_k t_k t_k1 t_k1];
        yd = [ylim(1) ylim(2) ylim(2) ylim(1)];
        fill(xd,yd,tc,'edgecolor','none'); hold on;
    end
    % Plot g-level
    plot(t_gl,g_gl,'-k'); hold on;
    xlabel('Elapsed Time(s)');
    ylabel('Acceleration or Vibration (g)');
    % Plot reads
    yyaxis right
    for k=1:numel(id)
        % ID for this read
        id_k = id(k);
        % plot this read
        if reads.within_period(id_k)
            read_color = [1 0 0];
        else
            read_color = [0.8 0.8 0.8];
        end
        plot([t_start_s(id_k) t_stop_s(id_k)],reads.basecalled_quality_phred(id_k)*[1 1],'-','color',read_color,'LineWidth',lw); hold on;
    end
    ylabel('Phred Quality Score');
    yyaxis left
    set(gca,'xlim',xlim,'ylim',ylim);
    plot(1:numel(rms),rms,'LineWidth',lw,'Color',[0 0 0.5]); hold on;
    
    % Save figure
    fig.PaperPositionMode = 'manual';
    orient(fig,'Landscape');
    fig.PaperUnits = 'inches'; w = 11; h = 6;
    fig.PaperPosition = [(11-w)/2 (8.5-h)/2 w h];
    fn = 'zerogseq_05_fig_rms_vibration_and_sequence_quality_Mars_Lunar-Europa_Zero_Focus';
    print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});
end

%% Prepare data for median quality and RMS

% Calculate median read quality as a function of time
t = 1:max([ceil(max(t_stop_s)) numel(rms)]);
q_med = NaN(size(t));
for k=1:numel(t)
    bMatch = bitand(t_start_s<=t(k),t_stop_s>=t(k));
    q_k = reads.basecalled_quality_phred(bMatch);
    if isempty(q_k)
        q_med(k) = NaN;
    else
        q_med(k) = median(q_k);
    end
end

% Pad rms with NaN if necessary
% Already padded q_med above
% Now rms and q_med are equal length
if numel(t)>numel(rms)
    rms(numel(rms):numel(t))=NaN;
end

%% Plot median quality and RMS (1 s bins) together

figure('color',[1 1 1]); set(gca,'fontsize',fontsize,'fontname',fontname);
plot(t,rms,'LineWidth',lw); hold on;
xlabel('Elapsed Time(s)');
ylabel('RMS vibration (g), 1s bin');
yyaxis right
plot(t,q_med,'-r','LineWidth',lw);
ylabel('Median Phred Quality Score');

% Save figure
fn = 'zerogseq_05_fig_rms_vibration_and_median_sequence_quality';
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

% Figure useful for visualizing 3-dimensional relationship
figure('color',[1 1 1]);  set(gca,'fontsize',fontsize,'fontname',fontname);
%t_max = numel(rms);
t_min = min(t);
plot3(q_med,rms,t,'.k')
xlabel('Median Quality Score');
ylabel('RMS vibration (g)');
zlabel('Time (s)');
grid on;

% Save figure
fn = 'zerogseq_05_fig_rms_vibration_median_sequence_quality_and_time';
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

% Load g-level data if available
g_1s = NaN(1,numel(t));
if ~isempty(g_data), 
    gdata = load(g_data);
    % Start time is identical to vibration data so no timing offset is required
    % Calculate g_filt for each 1 second interval
    for k=1:numel(t)
        bPeriod = bitand(gdata.t>=k-1,gdata.t<k);
        g_1s(k) = mean(gdata.g_filt(bPeriod));
    end
    
    % Stepwise Linear Regression using Time, RMS vibration, G-level
    disp('Stepwise Linear Regression using Time, RMS vibration, G-Level');
    X = [(t_min:t_max)' rms(t_min:t_max)' g_1s(t_min:t_max)'];
    y = q_med(1:find(t==t_max))';
    % Stepwise Linear Regression
    % Usage: [b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(X,y)
    mdl1 = stepwiselm(X,y,'VarNames',{'Time' 'Vibration' 'g-level' 'Read Quality'},'Verbose',2)
else
    % Stepwise Linear Regression using Time, RMS vibration
    disp('Stepwise Linear Regression using Time, RMS vibration');
    X = [(t_min:t_max)' rms(t_min:t_max)'];
    y = q_med(1:find(t==t_max))';
    % Stepwise Linear Regression
    % Usage: [b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(X,y)
    mdl1 = stepwiselm(X,y,'VarNames',{'Time' 'Vibration' 'Read Quality'},'Verbose',2)
end

diary off;
