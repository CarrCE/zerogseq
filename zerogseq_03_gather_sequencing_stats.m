% zerogseq_03_gather_sequencing_stats.m
%
% Collect sequencing statistics gathered from read tables and base tables,
% which were produced by zerogseq_02_make_basetables.m

% start fresh
clear all; close all; clc;

%% Parameters for coverage analysis
genome_length = 48502;
base_data = './analysis/MinION/Flight/flight_run/Bases.mat';

%% Plot options

% Graphics format for saving figures {file_extension, MATLAB print option}
% Type 'help print' for examples
figformat = {'eps' '-depsc' '-painters'};
%figformat = {'pdf' '-dpdf'};

% LineWidth for figure plotting
lw = 0.5;
fontsize = 7;
fontname = 'Helvetica';

%% Read in read data
fn = {'./analysis/MinION/Ground/ground_mux/Reads.mat', ...
      './analysis/MinION/Ground/ground_run/Reads.mat', ...
      './analysis/MinION/Flight/flight_mux/Reads.mat', ...
      './analysis/MinION/Flight/flight_run/Reads.mat'};

for k=1:4  
    load(fn{k});
    summarystats.ReadsFile{k,1} = fn{k};
    summarystats.N_reads(k,1) = numel(reads.filename);
    summarystats.N_basecalled_reads(k,1) = sum(~isnan(reads.basecalled_bases));
    summarystats.N_basecalled_reads_Q65(k,1) = sum([reads.basecalled_quality_phred>6.5]);
    summarystats.N_basecalled_bases(k,1) = nansum(reads.basecalled_bases);
    summarystats.N_tombo_reads(k,1) = sum(reads.tombo_has_data);
    summarystats.N_tombo_bases(k,1) = sum(reads.tombo_bases);
end

save('./analysis/MinION/SummaryStats.mat','-struct','summarystats');
writetable(struct2table(summarystats),'./analysis/MinION/SummaryStats.csv');

%% Coverage analysis of flight periods

% Start logging coverage results
if ~exist('./analysis/Combined/Flight'), mkdir('./analysis/Combined/Flight'); end
diary('./analysis/Combined/Flight/zerogseq_03_coverage.txt');

periods_file = './analysis/acceleration/periods.txt';
periods = readtable(periods_file);
period_ids = periods.period;

%% Coverage estimates based on reads within periods
%period_ids = [5 9 13 17 21];
fprintf('Coverage based on reads wholly within periods\n');
for k=1:numel(period_ids)
    period_id = period_ids(k);
    bPeriod = bitand(reads.start_period==period_id,reads.stop_period==period_id);
    bBasecalledPeriod = bitand(~isnan(reads.basecalled_bases),bPeriod);
    fns = reads.filename(bBasecalledPeriod);
    [S,Q]=extract_fast5_helper(fns);
    % estimate coverage
    basecount = sum(cellfun(@numel,S));
    cov_est = basecount/genome_length;
    % store results
    coverage.period(k,1) = period_id;
    coverage.reads_basecount(k,1) = basecount;
    coverage.reads_coverage(k,1) = cov_est;
    % Display for user
    fprintf('Period:\t%d\tBases:\t%d\tCoverage:\t%0.2f\n',period_id,basecount,cov_est);
end
fprintf('\n\n');

% Coverage estimates based on Tombo bases
% Load base info
bases = load(base_data,'within_period','start_period','stop_period');
fprintf('Coverage based on Tombo-aligned bases:\n');
for k=1:numel(period_ids)
    period_id = period_ids(k);
    bPeriod = bitand(bases.start_period==period_id,bases.stop_period==period_id);
    % Compute coverage
    basecount = sum(bPeriod);
    cov_est = basecount/genome_length;
    % store results
    coverage.period(k,1) = period_id;
    coverage.tombo_basecount(k,1) = basecount;
    coverage.tombo_coverage(k,1) = cov_est;
    % Display for user
    fprintf('Period:\t%d\tBases:\t%d\tCoverage:\t%0.2f\n',period_id,basecount,cov_est);
end
fprintf('\n\n');

save('./analysis/Combined/Flight/Coverage.mat','-struct','coverage');
writetable(struct2table(coverage),'./analysis/Combined/Flight/Coverage.csv');

%% Figure for coverage

% Periods for first five parabolas are 5, 9, 13, 17, 21
p = [5, 9, 13, 17, 21];

bases_reads = coverage.reads_basecount(p);
bases_tombo = coverage.tombo_basecount(p);
cov_reads = coverage.reads_coverage(p);
cov_tombo = coverage.tombo_coverage(p);

fig = figure('color',[1 1 1]);
set(gca,'fontsize',fontsize,'fontname',fontname);

X = categorical({'Mars1','Mars2','Lunar/Europa','Zero1','Zero2'});
X = reordercats(X,{'Mars1','Mars2','Lunar/Europa','Zero1','Zero2'});
Y = [cov_reads';cov_tombo'];
bar(X,Y);
ylabel('Coverage');

% Save figure
fig.PaperPositionMode = 'manual';
orient(fig,'Landscape');
fig.PaperUnits = 'inches'; w = 11; h = 2;
fig.PaperPosition = [(11-w)/2 (8.5-h)/2 w h];
fn = 'zerogseq_03_Mars_Lunar-Europa_Zero_coverage';
print(fullfile('./analysis/Combined/Flight',[fn '.' figformat{1}]),figformat{2});

%% Periods duration vs. coverage regression

% Extract data
p = periods.period(periods.parabola==1);
dt = periods.duration_s(p);
g = periods.g_bar_norm(p);
cov_tombo = coverage.tombo_coverage(p);

% Stepwise Linear Regression using parabola duration, tombo coverage
disp('Stepwise Linear Regression using Parabola Duration (s), Tombo Coverage');
X = dt;
y = cov_tombo;
% Regression
mdl1 = stepwiselm(X,y,'VarNames',{'Duration' 'Coverage'},'Verbose',2)

% Plot params
pointsize = 16;
% Do plot
fig = figure('color',[1 1 1]);
set(gca,'fontsize',fontsize,'fontname',fontname);
% Add regression line
beta = mdl1.Coefficients.Estimate;
dt_fit = 10:40;
cov_fit = beta(1) + beta(2)*dt_fit;
plot(dt_fit,cov_fit,'-k'); hold on;
% Mars
scatter(dt(1:2),cov_tombo(1:2),pointsize,'r^','filled'); hold on;
% Lunar
scatter(dt(3),cov_tombo(3),pointsize,'bs','filled'); hold on;
% Zero
scatter(dt(4:end),cov_tombo(4:end),pointsize,'ko','filled'); hold on;

% Labels
xlabel('Parabola duration (s)'); ylabel('Coverage');

% Save figure
fig.PaperUnits = 'inches'; w = 2; h = 2.5;
fig.PaperPosition = [(8.5-w)/2 (11-h)/2 w h];
fn = 'zerogseq_03_Mars_Lunar-Europa_Zero_coverage_regression';
print(fullfile('./analysis/Combined/Flight',[fn '.' figformat{1}]),figformat{2});

% Legend
legend({'Fit','Mars','Lunar/Europa','Zero-g'},'Location','NorthWest');
% Save figure legend
fn = 'zerogseq_03_Mars_Lunar-Europa_Zero_coverage_regression_legend';
print(fullfile('./analysis/Combined/Flight',[fn '.' figformat{1}]),figformat{2});

% Stop logging
diary off;

