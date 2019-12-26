% extract_fast5_helper extracts sequence data from a FAST5 file. By
% default, data is extracted from 
% /Analyses/Basecall_1D_000/BaseCalled_template/Fastq
% and parsed into sequencing data and quality data.
function [S,Q]= extract_fast5_helper(fns,entrypoint)
    if nargin<2
        entrypoint = '/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'; 
    end

    % handle argument as single file specified as string
    if ischar(fns), fns = {fns}; end
    
    % preallocate memory
    N = numel(fns);
    S = repmat({''},N,1);
    Q = repmat({''},N,1);

    % grab data for each file in file list
    for k=1:numel(fns)
        fq = h5read(fns{k},entrypoint);
        rows = regexp(fq,'\n','split');
        S{k} = rows{2};
        Q{k} = rows{4};
    end
end


