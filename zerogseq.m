% zerogseq.m
addpath('./code');

zerogseq_00_download;
zerogseq_01_unzip;
zerogseq_02_make_basetables;
zerogseq_03_gather_sequencing_stats;

clear all; close all; clc; dataset = 'Ground';
zerogseq_04_vibration_analysis;

clear all; close all; clc; dataset = 'Flight';
zerogseq_04_gather_sequencing_stats

clear all; close all; clc; dataset = 'Ground';
zerogseq_05_vibration_and_sequence_quality

clear all; close all; clc; dataset = 'Flight';
zerogseq_05_vibration_and_sequence_quality

clear all; close all; clc; dataset = 'Ground';
zerogseq_06_vibration_and_sequence_quality

clear all; close all; clc; dataset = 'Flight';
zerogseq_06_vibration_and_sequence_quality

clear all; close all; clc;
zerogseq_07_translocation_times