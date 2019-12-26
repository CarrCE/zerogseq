% zerogseq_07_translocation_times.m

% Uncomment below two lines if using this as a standalone script
% clear all; close all; clc;
% addpath('./code');

outfolder = './analysis/MinION/Both';

% Check if output folder exists
if ~(exist(outfolder,'dir')==7), mkdir(outfolder); end

% Graphics format for saving figures {file_extension, MATLAB print option}
% Type 'help print' for examples
figformat = {'eps' '-depsc' '-painters'};
%figformat = {'pdf' '-dpdf'};

% LineWidth for figure plotting
lw = 0.5;
fontsize = 24;
fontname = 'Helvetica';

% get ground duration in samples
d_g_mux = load('./analysis/MinION/Ground/ground_mux/Bases.mat','duration_samples');
d_g_run = load('./analysis/MinION/Ground/ground_run/Bases.mat','duration_samples');
d_g = [d_g_mux.duration_samples;d_g_run.duration_samples];
d_g = [d_g_run.duration_samples];

% get flight duration in samples
d_f_mux = load('./analysis/MinION/Flight/flight_mux/Bases.mat','duration_samples');
d_f_run = load('./analysis/MinION/Flight/flight_run/Bases.mat','duration_samples');
d_f = [d_f_mux.duration_samples;d_f_run.duration_samples];
d_f = [d_f_run.duration_samples];

% sample frequency
fs = 4000;

% edges for histogram in samples
edges = 0:10000;
% get histogram counts in samples
h1 = histcounts(d_g,edges);
h2 = histcounts(d_f,edges);
% edges for plotting in milliseconds
edges_ms = edges/fs*1000; 

figure('color',[1 1 1]);
ea = 0.05;
%h = bar(edges_ms(1:end-1),h1./sum(h1),'b'); hold on;
%h = bar(edges_ms(1:end-1),h2./sum(h2),'r');
histogram('BinEdges',edges_ms,'BinCounts',h1./sum(h1),'EdgeAlpha',ea); hold on;
histogram('BinEdges',edges_ms,'BinCounts',h2./sum(h2),'EdgeAlpha',ea)
set(gca,'yscale','log')
legend('ground','flight');
set(gca,'xlim',[0 2500]);
xlabel('Translocation time (ms)');
ylabel('Frequency');

% Save figure
fn = 'zerogseq_07_fig_translocation_times';
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

figure('color',[1 1 1]);
ea = 1;
%h = bar(edges_ms(1:end-1),h1./sum(h1),'b'); hold on;
%h = bar(edges_ms(1:end-1),h2./sum(h2),'r');
histogram('BinEdges',edges_ms,'BinCounts',h1./sum(h1),'EdgeAlpha',ea); hold on;
histogram('BinEdges',edges_ms,'BinCounts',h2./sum(h2),'EdgeAlpha',ea)
legend('ground','flight');
set(gca,'xlim',[0 10]);
xlabel('Translocation time (ms)');
ylabel('Frequency');
set(gca,'fontsize',fontsize);

% Save figure
fn = 'zerogseq_07_fig_translocation_times_zoom';
print(fullfile(outfolder,[fn '.' figformat{1}]),figformat{2});

% Statistical test:
% Kolmogorov-Smirnov tests whether two sets of samples are from the same
% arbitrary distribution (H0) or different (H1) distributions.
% https://www.mathworks.com/help/stats/kstest2.html
[h,p,ks2stat] = kstest2(d_g,d_f,'tail','unequal')

% results with run (no mux data), 'tail','unequal'
% h = 1; (different)
% p = 0; (highly significant)
% ks2stat = 0.0306;

% results with run (no mux data), 'tail','smaller'
% h = 0; (cdf of ground not smaller than flight)
% p = 0; (highly significant)
% ks2stat = 0;

% results with run (no mux data), 'tail','larger'
% h = 1; (cdf of ground larger than flight, e.g. left shifted translocation
% times thus shorter translocation times)
% p = 0; (highly significant)
% ks2stat = 0.0306;




