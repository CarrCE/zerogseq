% zerogseq_06_vibration_and_ionic_current_noise.m
%
% Analyze the relationship between RMS vibration and ionic current noise as
% measured by Tombo and stored in base tables.
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

%% Choose correct options

% FLIGHT 
if strcmpi(dataset,'Flight')
    vib_data = './analysis/vibration/Flight/zerogseq_04_rms_vibration.mat';
    read_data = './analysis/MinION/Flight/flight_run/Reads.mat';
    base_data = './analysis/MinION/Flight/flight_run/Bases.mat';
    run_offset_s = 1061;
    outfolder = './analysis/Combined/Flight';
    bPeriods = true;
    periods_file = './analysis/acceleration/periods.txt';
    % Regression restricted to just before and just after parabolas to 
    % avoid complicating factors of g-level and vibration changes during 
    % aircraft descent
    t_max = 4000; 
    g_data = '/Users/chrisc/Desktop/ZeroG-Tombo/Accel/analysis/G_filt.mat';
    logfile = fullfile(outfolder,'zerogseq_06_Flight_vibration_and_ionic_current_noise.txt');
    fig_prefix = 'zerogseq_06_Flight';
end

% GROUND
if strcmpi(dataset,'Ground')
    vib_data = './analysis/vibration/Ground/zerogseq_04_rms_vibration.mat';
    read_data = './analysis/MinION/Ground/ground_run/Reads.mat';
    base_data = './analysis/MinION/Ground/ground_run/Bases.mat';
    run_offset_s = 406;
    outfolder = './analysis/Combined/Ground';
    bPeriods = false;
    periods_file = '';
    t_max = 2266; % Regression across all data
    g_data = '';
    logfile = fullfile(outfolder,'zerogseq_06_Ground_vibration_and_ionic_current_noise.txt');
    fig_prefix = 'zerogseq_06_Ground';
end

% start loggin
diary(logfile);

% Check if output folder exists
if ~(exist(outfolder,'dir')==7), mkdir(outfolder); end

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
% Load base info
bases = load(base_data);

% Calculate Tombo Read times in accelerometer elapsed time
t_start_s = bases.t_start_s + run_offset_s;
t_stop_s = bases.t_stop_s + run_offset_s;

% Plot
figure('color',[1 1 1]);
set(gca,'fontsize',fontsize,'fontname',fontname);
plot(rms,'LineWidth',lw); hold on;
xlabel('Elapsed Time(s)');
ylabel('RMS vibration (g), 1s bin');
yyaxis right
h = scatter(t_start_s,bases.norm_std,'.k','MarkerFaceColor','k');
h.MarkerFaceAlpha = 0.01;
h.MarkerEdgeAlpha = 0.01;
ylabel('Ionic Current Noise (norm_{std})');

% Save figure
fn = [fig_prefix '_fig_rms_vibration_and_ionic_current_noise'];
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

% We observe a repetitive change in the ionic current noise that has a
% frequency of around 20 to 25 seconds. Effects we are considering:
% 1. Nanopore flicking (seems too short; typical flick timing is 10 min+)
% 2. EM noise due to aircraft (beat frequency?) or some kind of radar.
% 3. Other?

% Calculate median ionic current as a function of time
t = 1:max([ceil(max(t_stop_s)) numel(rms)]);
norm_std_med = NaN(size(t));
n_bases = NaN(size(t));
for k=1:numel(t)
    bMatch = bitand(t_start_s<=t(k),t_stop_s>=t(k));
    norm_std_k = bases.norm_std(bMatch);
    if isempty(norm_std_k)
        norm_std_med(k) = NaN;
        n_bases(k) = 0;
    else
        norm_std_med(k) = median(norm_std_k);
        n_bases(k) = numel(norm_std_k);
    end
end

% Look up period as a function of time if we are using a Periods table
if bPeriods
    periods = readtable(periods_file);
    time.period = zeros(size(t));
    time.transition = false(size(t));
    time.parabola = false(size(t));
    time.hypergravity = false(size(t));
    time.other = false(size(t));
    for k=1:numel(t)
        % get period for current time
        bMatch = bitand(t(k)>=periods.time_start,t(k)<=periods.time_stop);
        id = find(bMatch,1);
        if ~isempty(id)
            time.period(k) = periods.period(id);
            time.transition(k) = periods.transition(id);
            time.parabola(k) = periods.parabola(id);
            time.hypergravity(k) = periods.hypergravity(id);
            time.other(k) = ~bitor(periods.transition(id),bitor(periods.parabola(id),periods.hypergravity(id)));
        else
            % already assigned zero (period) or false (other fields)
        end
    end
end

if bPeriods
    % Timeline figure
    hf = figure('color',[1 1 1]); set(gca,'fontsize',fontsize,'fontname',fontname);
    sz = 3;
    h = plot(t(time.parabola),time.parabola(time.parabola)+0,'.b'); hold on;
    h.MarkerSize = sz;
    h = plot(t(time.transition),time.transition(time.transition)+0.02,'.m'); hold on;
    h.MarkerSize = sz;
    h = plot(t(time.hypergravity),time.hypergravity(time.hypergravity)+0.06,'.r'); hold on;
    h.MarkerSize = sz;
    h = plot(t(time.other),time.other(time.other)+.04,'.k'); hold on;
    h.MarkerSize = sz;
    xlabel('Elapsed time (s)');
    set(gca,'ylim',[0.5 1.5]);

    % Save figure
    fn = [fig_prefix '_fig_timeline'];
    print(hf,fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

    set(gca,'ylim',[0.94 1.58]);
    set(gca,'xlim',[1100 1550]);
    set(gca,'ytick',[]);
    set(gca,'ycolor',[1 1 1]);
    
    % Save figure
    fn = [fig_prefix '_fig_timeline_zoom'];
    print(hf, fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});
    % pause
end

% Pad rms with NaN if necessary
% Already padded norm_std_med above
% Now rms and norm_std_med are equal length
if numel(t)>numel(rms)
    rms(numel(rms):numel(t))=NaN;
end

figure('color',[1 1 1]); set(gca,'fontsize',fontsize,'fontname',fontname);
plot(t,rms,'LineWidth',lw); hold on;
xlabel('Elapsed Time(s)');
ylabel('RMS vibration (g), 1s bin');
yyaxis right
plot(t,n_bases,'-r','LineWidth',lw);
ylabel('Number of aligned bases');

% Save figure
fn = [fig_prefix '_fig_rms_vibration_and_number_of_aligned_bases'];
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

figure('color',[1 1 1]); set(gca,'fontsize',fontsize,'fontname',fontname);
plot(t,rms,'LineWidth',lw); hold on;
xlabel('Elapsed Time(s)');
ylabel('RMS vibration (g), 1s bin');
yyaxis right
plot(t,norm_std_med,'-r','LineWidth',lw);
ylabel('Median Ionic Current Noise (norm_{std})');

% Save figure
fn = [fig_prefix '_fig_rms_vibration_and_median_ionic_current_noise'];
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

fig=figure('color',[1 1 1]); set(gca,'fontsize',fontsize,'fontname',fontname);
if bPeriods
    % plot separately for transition, parabola, hypergravity
    al = 0.3;
    sz = 75;    
      
    h=scatter(rms(time.other),norm_std_med(time.other),'.k'); hold on;
    h.MarkerFaceAlpha = al;
    h.MarkerEdgeAlpha = al;
    h.SizeData = sz;
    h=scatter(rms(time.transition),norm_std_med(time.transition),'om'); hold on;
    h.MarkerFaceAlpha = al;
    h.MarkerEdgeAlpha = al;
    h.SizeData = sz/4;
    h=scatter(rms(time.hypergravity),norm_std_med(time.hypergravity),'.r'); hold on;
    h.MarkerFaceAlpha = al;
    h.MarkerEdgeAlpha = al;
    h.SizeData = sz;
    h=scatter(rms(time.parabola),norm_std_med(time.parabola),'.b'); hold on;
    h.MarkerFaceAlpha = al;
    h.MarkerEdgeAlpha = al;
    h.SizeData = sz;
    legend('Other','Transition','Hypergravity','Parabola');
else
    h=scatter(rms,norm_std_med,'.k'); hold on;
    h.MarkerFaceAlpha = 0.3;
    h.MarkerEdgeAlpha = 0.3;
end
xlabel('RMS Vibration (g)'); % 1s bin (caption)
ylabel('Median Ionic Current Noise (norm_{std})');

% Save figure
fn = [fig_prefix '_fig_rms_vibration_vs_median_ionic_current_noise'];
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

% Save figure (zoom)
set(gca,'ylim',[0.15 0.4]);
fig.PaperUnits = 'inches'; w = 5; h = 4;
fig.PaperPosition = [(8.5-w)/2 (11-h)/2 w h];
fn = [fig_prefix '_fig_rms_vibration_vs_median_ionic_current_noise_zoom'];
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

% Figure useful for visualizing 3-dimensional relationship
figure('color',[1 1 1]);  set(gca,'fontsize',fontsize,'fontname',fontname);
%t_max = numel(rms);
t_min = min(t);
plot3(norm_std_med,rms,t,'.k')
xlabel('Ionic Current Noise (norm_{std})');
ylabel('RMS vibration (g)');
zlabel('Time (s)');
grid on;

% Save figure
fn = [fig_prefix '_fig_rms_vibration_median_ionic_current_noise_and_time'];
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

% Any relationship between median ionic current noise and  number of bases?
figure('color',[1 1 1]); set(gca,'fontsize',fontsize,'fontname',fontname);
scatter(n_bases,norm_std_med,'ok'); hold on;
xlabel('Number of aligned bases (1s bin)');
ylabel('Median Ionic Current Noise (norm_{std})');

% Save figure
fn = [fig_prefix '_fig_median_ionic_current_noise_vs_number_of_aligned_bases'];
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

% Load g-level data if available and perform appropriate stepwise linear
% regression
g_1s = NaN(1,numel(t));
if ~isempty(g_data)
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
    y = norm_std_med(1:find(t==t_max))';
    % Usage: [b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(X,y)
    mdl1 = stepwiselm(X,y,'VarNames',{'Time' 'Vibration' 'g-level' 'Ionic Current Noise'},'Verbose',2)
else
    % Stepwise Linear Regression using Time, RMS vibration
    disp('Stepwise Linear Regression using Time, RMS vibration');
    X = [(t_min:t_max)' rms(t_min:t_max)'];
    y = norm_std_med(1:find(t==t_max))';
    % Usage: [b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(X,y)
    mdl1 = stepwiselm(X,y,'VarNames',{'Time' 'Vibration' 'Ionic Current Noise'},'Verbose',2)
end

diary off;

% Identify outliers in ionic current noise data (norm_std_med)
% Considered: [norm_std_med,bRemoved]=rmoutliers(norm_std_med,'method')
% Used: median aboslute deviation (standard outlier definition is MAD>3)
norm_std_med_mad = mad(norm_std_med,1);
score = (norm_std_med-nanmedian(norm_std_med))/norm_std_med_mad;
bOut = score>15;
sum(bOut)
% Based on this score, there is only one outlier in flight data

if strcmpi(dataset,'Flight')

%% Does g-level/phase of flight have an impact on read quality?

% Use periods information to compare across groups of Parabola, Transition,
% Hypergravity, and Other

% Do this on the basis of reads, because quality is defined at the read
% level. Compare quality for all reads wholly within a given type of phase
% of flight.
load(read_data)
% Categorize reads as 0=parabola, 1=other(1g), 2=hypergravity,
% 3=transition, NaN=not completely within a phase of flight
bOther = bitand(reads.within_period,~bitor(bitor(reads.transition,reads.parabola),reads.hypergravity));
group = bOther + 0*reads.parabola + 2*reads.hypergravity + 3*reads.transition;
group(~reads.within_period)=NaN;
x = reads.basecalled_quality_phred;
% Exclude transition due to low number (7) and short length, which will
% bias results towards poor quality reads
group(reads.transition) = NaN;
[p,tbl,stats]=anova1(x(~isnan(group)),group(~isnan(group)))
set(gca,'xticklabels',{'Parabola' 'Other (1g)' 'Hypergravity'});
set(gcf,'color',[1 1 1]);
ylabel('Read Quality (Phred)');

% Save figure
fn = [fig_prefix '_fig_read_quality_Flight_ANOVA'];
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

figure; [c,m,h,nms] = multcompare(stats)
set(gca,'yticklabels',{'Hypergravity' 'Other (1g)' 'Parabola'});
set(gcf,'color',[1 1 1]);
xlabel('Read Quality (Phred)');

% Save figure
fn = [fig_prefix '_fig_read_quality_Flight_HSD'];
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});


%% Does g-level/phase of flight have an impact on ionic current noise?

% Use periods information to compare across groups of Parabola, Transition,
% Hypergravity, and Other

% Do this on the basis of bases, because ionic current noise is defined at
% the base level. Compare ionic current quality for all bases wholly within
% a given phase of flight.

% Categorize reads as 0=parabola, 1=other(1g), 2=hypergravity,
% 3=transition, NaN=not completely within a phase of flight
bOther = bitand(bases.within_period,~bitor(bitor(bases.transition,bases.parabola),bases.hypergravity));
group = bOther + 0*bases.parabola + 2*bases.hypergravity + 3*bases.transition;
group(~bases.within_period)=NaN;
x = bases.norm_std;
disp('ANOVA on ionic current, groups = phases of flight');
[p,tbl,stats]=anova1(x(~isnan(group)),group(~isnan(group)))
set(gca,'xticklabels',{'Parabola' 'Other (1g)' 'Hypergravity' 'Transition'});
set(gcf,'color',[1 1 1]);
ylabel('Ionic Current Noise (1s bin)');

% Save figure
fn = [fig_prefix '_fig_ionic_current_noise_Flight_ANOVA'];
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

figure; [c,m,h,nms] = multcompare(stats)
set(gca,'yticklabels',{'Transition' 'Hypergravity' 'Other (1g)' 'Parabola'});
set(gcf,'color',[1 1 1]);
xlabel('Ionic Current Noise (1s bin)');
% pause;

% Save figure
fn = [fig_prefix '_fig_ionic_current_noise_Flight_HSD'];
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

end % Flight dataset

