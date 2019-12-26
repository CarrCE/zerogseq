% BASETABLE.m
%
% Makes a base table for a set of FAST5 files that have been processed
% using the TOMBO RESQUIGGLE command to produce internal alignments to a
% genomic reference. In our case, these internal alignments are used to
% derive a table that contains timing and base identity data for each base
% in the set of reads.
%
% Authors: Christopher E. Carr, Noelle Bryan
% 2019/01/15
%

function [reads,bases,options] = basetable(in,out,varargin)

    %% OPTIONS
    % Set any user specified options
    useroptions = args2options(varargin);
    % Set all options to defaults or user specified options
    options = []; % initial empty options
    options = fieldcheck(options,'base_table', fullfile(out,'Bases.csv'), useroptions);
    options = fieldcheck(options,'base_table_mat', fullfile(out,'Bases.mat'), useroptions);
    options = fieldcheck(options,'read_table', fullfile(out,'Reads.csv'), useroptions);
    options = fieldcheck(options,'read_table_mat', fullfile(out,'Reads.mat'), useroptions);
    options = fieldcheck(options,'periods_file','',useroptions);
    options = fieldcheck(options,'periods_offset',NaN,useroptions); % added to seq time to get periods time
    
    % Define default outputs
    reads = struct([]);
    bases = struct([]);
    
    % see if we have periods into which we should bin the reads and bases
    bPeriods = false;
    if and(~isempty(options.periods_file),~isnan(options.periods_offset))
        % Check that file exists
        if exist(options.periods_file,'file')
           % Read in the file
           periods = readtable(options.periods_file);
           bPeriods = true;
        else
            % File doesn't exist
            warning('basetable: period_file %s not found',options.period_file);
        end        
    end
    
    % Verify output folder    
    if ~exist(out,'dir') 
        mkdir(out);
        
        % Require input folder to exist
        assert(exist(in,'dir')>0,'Input folder %s not found',in);
        
        % Get set of fast5 files to process from the input directory
        filespec = fullfile(in,'**/*.fast5');
        
        % Get list of fast5 files to process
        d = dir(filespec);
        
        % number of files
        N = numel(d);

        %% Build Read Table
        reads(1).filename = repmat({''},[N 1]);
        reads.experiment_start_time = repmat({''},[N 1]);
        reads.sampling_rate_Hz = NaN([N 1]);
        reads.start_time_samples = NaN([N 1]);
        reads.duration_samples = NaN([N 1]);
        reads.read_number = NaN([N 1]);
        reads.read_id = repmat({''},[N 1]);
        reads.basecalled_bases = NaN([N 1]);
        reads.basecalled_quality_phred = NaN([N 1]);   
        reads.tombo_has_data = NaN([N 1]);
        reads.tombo_alignment_offset_samples = NaN([N 1]);
        reads.tombo_duration_samples = NaN([N 1]);
        reads.tombo_bases = NaN([N 1]);
        
        reads.start_period = zeros(N,1,'uint8');
        reads.stop_period = zeros(N,1,'uint8');
        reads.within_period = false(N,1);
        reads.transition = false(N,1);
        reads.parabola = false(N,1);
        reads.hypergravity = false(N,1);
        
        % Loop through all files, store read info, deterine which have
        % valid tombo alignments and total number of bases.
        for k=1:N
            % Process kth file
            % Get filename
            fn = fullfile(d(k).folder,d(k).name);
            
            % Get Experiment Start Time
            experiment_start_time=h5readatt(fn,'/UniqueGlobalKey/tracking_id','exp_start_time');
            
            % Get Sampling Rate
            sampling_rate_Hz = h5readatt(fn,'/UniqueGlobalKey/channel_id','sampling_rate');

            % Read start time
            info = h5info(fn,'/Raw/Reads');
            key = info.Groups.Name;
            read_start_time = h5readatt(fn,key,'start_time');
            read_duration_samples = h5readatt(fn,key,'duration');
            read_number = h5readatt(fn,key,'read_number');
            read_id = h5readatt(fn,key,'read_id');
            
            % convert filename to local filename
            lp = cd;
            fl = ['.' fn((numel(lp)+1):end)];
            
            % Save results expected to be in all reads
            reads.filename(k) = {fl};
            reads.experiment_start_time(k) = {experiment_start_time};
            reads.sampling_rate_Hz(k) = sampling_rate_Hz;
            reads.start_time_samples(k) = read_start_time;
            reads.duration_samples(k) = read_duration_samples;
            reads.read_number(k) = read_number;
            reads.read_id(k) = {deblank(read_id)};
            
            % Get basecalled read attributes
            try
                data = hdf5read(fn,'/Analyses/Basecall_1D_000/BaseCalled_template/Fastq');
                fq = data.Data;
                rows = regexp(fq,'\n','split');
                seq = rows{2};
                qual = rows{4};
                % Process quality data
                % Phred + 33 encoding as documented here: https://goo.gl/rj2LwB
                q = qual-33;
                % Transform to per base error probability
                p = 10.^(-q/10);
                % Estimate mean per base error probability
                p_bar = mean(p);
                % Transform back to phred quality score
                q_bar = -10*log10(p_bar);
                % Add statistics
                reads.basecalled_bases(k) = numel(seq);
                reads.basecalled_quality_phred(k) = q_bar;
            catch
                % Could not load fastq; basecalling failed on this read
                reads.basecalled_bases(k) = NaN;
                reads.basecalled_quality_phred(k) = NaN;
            end
            
            % Tombo data may or may not exist in this read
            try
                % Read in tombo data example
                data = h5read(fn,'/Analyses/RawGenomeCorrected_000/BaseCalled_template/Events');
                
                % Get zero-based offset indicating the beginning of the read genomic
                % sequence within the raw signal, read_start_rel_to_raw
                info = h5info(fn,'/Analyses/RawGenomeCorrected_000/BaseCalled_template/Events');
                tombo_alignment_offset_samples = info.Attributes.Value;
                
                % Save tombo data
                reads.tombo_has_data(k) = true;
                reads.tombo_alignment_offset_samples(k) = tombo_alignment_offset_samples; 
                reads.tombo_duration_samples(k) = sum(data.length);
                reads.tombo_bases(k) = numel(data.base);
            catch
                % No tombo data in read, or failed above for another reason 
                reads.tombo_has_data(k) = false;
                reads.tombo_alignment_offset_samples(k) = NaN; 
                reads.tombo_duration_samples(k) = NaN;
                reads.tombo_bases(k) = 0;
            end
            
            % If appropriate, check whether this read falls within a
            % specific period
            if bPeriods
                % Put start and stop time of read on periods timeline
                t_read_start = (reads.start_time_samples(k) ./ reads.sampling_rate_Hz(k)) + options.periods_offset;
                t_read_stop = t_read_start + (reads.duration_samples(k) ./ reads.sampling_rate_Hz(k));
                % Check the read against the periods
                % Does the read fall entirely within a period or not?
                start_period = find(periods.time_start<=t_read_start,1,'last');
                stop_period = find(periods.time_stop>=t_read_stop,1,'first');
                % Handle case where start_period or stop_period is not found
                if isempty(start_period), start_period = NaN; end
                if isempty(stop_period), stop_period = NaN; end
                % Store results
                reads.start_period(k) = uint8(start_period);
                reads.stop_period(k) = uint8(stop_period);    
                % Check if start_period = stop_period (in the same period)
                if start_period==stop_period
                    % Read is completely within the start period = stop
                    % period. Mark period as the specific period it is
                    % within.
                    reads.within_period(k) = true;
                    bRow = (periods.period==start_period);
                    reads.transition(k) = periods.transition(bRow);
                    reads.parabola(k) = periods.parabola(bRow);
                    reads.hypergravity(k) = periods.hypergravity(bRow);
                else
                    reads.within_period(k) = false;
                    reads.transition(k) = false;
                    reads.parabola(k) = false;
                    reads.hypergravity(k) = false;
                end
            end
        end
        
        % Save read table
        save(options.read_table_mat,'reads','-v7.3');
        writetable(struct2table(reads),options.read_table);
        
        %% Build Base Table
        
        % Determine number of bases in base table
        N_b = nansum(reads.tombo_bases);
        % Initialize base table
        bases(1).read_id = NaN(N_b,1);
        bases.basecalled_length = NaN(N_b,1);
        bases.basecalled_quality = NaN(N_b,1);
        bases.base_id = NaN(N_b,1);
        bases.norm_mean = NaN(N_b,1);
        bases.norm_std = NaN(N_b,1);
        bases.base = repmat(' ',[N_b,1]);
        bases.duration_samples = NaN(N_b,1);
        bases.t_start_s = NaN(N_b,1);
        bases.t_stop_s = NaN(N_b,1);
        
        bases.start_period = zeros(N_b,1,'uint8');
        bases.stop_period = zeros(N_b,1,'uint8');
        bases.within_period = false(N_b,1);
        bases.transition = false(N_b,1);
        bases.parabola = false(N_b,1);
        bases.hypergravity = false(N_b,1);

        % Initialize base counter
        ctr = 0;
        
        % Process kth file
        for k=1:N
            % Check if has tombo data
            if reads.tombo_has_data(k)
                % Has tombo data; process this read
                % Get filename
                fn = reads.filename{k};
                % Get Sampling Rate
                fs = reads.sampling_rate_Hz(k);
                % Get Read Start Time in Samples
                t_read_start = reads.start_time_samples(k);
                % Calculate Read Start Time in Seconds
                t_read_start_s = t_read_start / reads.sampling_rate_Hz(k);
                % Get Tombo Alignment Offset Time in Samples
                t_tombo_offset = reads.tombo_alignment_offset_samples(k);
                % Calculate Alignment Offset Time in Seconds
                t_tombo_offset_s = t_tombo_offset / reads.sampling_rate_Hz(k);
                % Get number of bases
                N_bases_k = reads.tombo_bases(k);
                % Read in data for all bases
                data = h5read(fn,'/Analyses/RawGenomeCorrected_000/BaseCalled_template/Events');
                % data is a struct with fields norm_mean, norm_stdev, start, length, and base
                for j = 1:N_bases_k
                    % Augment base counter
                    ctr = ctr + 1;
                    % Create base entry for each base in the tombo alignment
                    bases.read_id(ctr) = k;
                    bases.basecalled_length(ctr) = reads.basecalled_bases(k);
                    bases.basecalled_quality(ctr) = reads.basecalled_quality_phred(k);
                    bases.base_id(ctr) = j;
                    bases.norm_mean(ctr) = data.norm_mean(j);
                    bases.norm_std(ctr) = data.norm_stdev(j);
                    bases.base(ctr) = char(data.base(j));
                    bases.duration_samples(ctr) = data.length(j);
                    % Calculate base start time in seconds since start of experiment
                    t_base_offset_s = double(data.start(j)) / reads.sampling_rate_Hz(k);
                    bases.t_start_s(ctr) = t_read_start_s + t_tombo_offset_s + t_base_offset_s;
                    % Calculate base stop time in seconds since start of experiment
                    t_base_duration_s = double(data.length(j)) / reads.sampling_rate_Hz(k);
                    bases.t_stop_s(ctr) = bases.t_start_s(ctr) + t_base_duration_s;                  
                    
                    % If appropriate, assign this base to a specific period
                    if bPeriods
                        % Put start and stop time of read on periods timeline
                        t_base_start = bases.t_start_s(ctr) + options.periods_offset;
                        t_base_stop = t_base_start + t_base_duration_s;
                        % Check the base against the periods
                        % Does the base fall entirely within a period or not?
                        start_period = find(periods.time_start<=t_base_start,1,'last');
                        stop_period = find(periods.time_stop>=t_base_stop,1,'first');
                        % Handle case where start_period or stop_period is not found
                        if isempty(start_period), start_period = NaN; end
                        if isempty(stop_period), stop_period = NaN; end
                        % Store results
                        bases.start_period(ctr) = uint8(start_period);
                        bases.stop_period(ctr) = uint8(stop_period);    
                        % Check if start_period = stop_period (in the same period)
                        if start_period==stop_period
                            % Base is completely within the start period = stop
                            % period. Mark period as the specific period it is
                            % within.
                            bases.within_period(ctr) = true;
                            bRow = (periods.period==start_period);
                            bases.transition(ctr) = periods.transition(bRow);
                            bases.parabola(ctr) = periods.parabola(bRow);
                            bases.hypergravity(ctr) = periods.hypergravity(bRow);
                        else
                            bases.within_period(k) = false;
                            bases.transition(k) = false;
                            bases.parabola(k) = false;
                            bases.hypergravity(k) = false;
                        end
                    end
                    
                end
            else
                % Has no tombo data; ignore this read
            end
        end
        
        % Save base table
        %save(options.base_table_mat,'bases','-v7.3');
        save(options.base_table_mat,'-struct','bases');
        writetable(struct2table(bases),options.base_table);
        
    else
        % Output folder exists
        fprintf('Output folder %s exists.\n',out);
        fprintf('Please remove existing folder or modify output folder name and rerun.\n\n',out);
    end
end
