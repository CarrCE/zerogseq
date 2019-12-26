% zerogseq_02_make_basetables.m

% Start fresh
clear all; close all; clc;

% Add code to our path
addpath('./code');

% Periods table to use in binning reads and bases for flight
periods_file = './analysis/acceleration/periods.txt';    

% offset: add to sequencing time to get accel time
flight_mux_offset = 329; 
flight_run_offset = 1061;

% Make the read tables and base tables

% ground_mux
basetable('./data/MinION/Ground/ground_mux','./analysis/MinION/Ground/ground_mux');

% ground_run
basetable('./data/MinION/Ground/ground_run','./analysis/MinION/Ground/ground_run');

% flight_mux
basetable('./data/MinION/Flight/flight_mux','./analysis/MinION/Flight/flight_mux',...
    'periods_file',periods_file,'periods_offset',flight_mux_offset);

% flight_run
basetable('./data/MinION/Flight/flight_run','./analysis/MinION/Flight/flight_run',...
    'periods_file',periods_file,'periods_offset',flight_run_offset);