# ZeroGSeq
Code for working with nanopore sequencing data collected in altered gravitointertial environments, such as during parabolic flight. For details on processing acceleration data from parabolic flights, see the related repository: <https://github.com/CarrCE/zerog>

# Compatibility
Tested with MATLAB 2019b on OS X 10.13. Expected to work with MATLAB 2019a. Requires the Statistics and Bioinformatics toolboxes.

# Citation
Carr CE, Bryan NC, Saboda K, Bhattaru SA, Ruvkun G, Zuber MT. Nanopore Sequencing at Mars, Europa, and Microgravity Conditions. In preparation.

See also: Carr CE, Bryan NC, Saboda KN, Bhattaru SA, Ruvkun G, Zuber MT. Acceleration Profiles and Processing Methods for Parabolic Flight. npj Microgravity 4, Article number: 14 (2018) <https://doi.org/10.1038/s41526-018-0050-3>. Preprint: arXiv:1712.05737 <https://arxiv.org/abs/1712.05737>

# Installation
## Get Scripts
Download: <https://github.com/CarrCE/zerogseq/archive/master.zip> or use command line ```git clone git@github.com:CarrCE/zerogseq.git```.

Unzip to preferred location, here denoted ```/zerogseq-master```.

## Get Data
Data is downloaded automatically by a matlab script. This is the recommended approach to ensure data is in the expected locations.

Direct access to accelerometer data:
Download: <https://osf.io/5rqu9/download>. This 1.0 GB (compressed ZIP) dataset has a CC BY 4.0 US license. More details at: <https://osf.io/nk2w4/>

Direct access to accelerometer analysis results:
The analysis results are also available (188 MB ZIP) at: <https://osf.io/pmhj4/download>

## Run Analysis
In MATLAB, go to your ```/zerogseq-master``` path, and run the main script: ```zerogseq_00_download```. To perform the same analysis as in the publication (see citation, above), run each script in turn. Some require changing the dataset between "Flight" or "Ground" options. See each script for details and instructions.

The results of running this analysis in MATLAB include a series of EPS and/or PDF figures, replicating those in the paper, and various tab-delimited files. All times are elapsed time, and for reference, the start time is: 2017-11-17 18:28:51 UTC.

# License
Distributed under an MIT license. See LICENSE for details.

