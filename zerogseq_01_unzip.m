% zerogseq_01_unzip.m

%% Unzip data from OSF

% Unzip Accelerometer Data Analysis File, generated from analysis of
% zerog parabolic data using the https://github.com/CarrCE/zerog package. 
% These results are available at: https://osf.io/pmhj4/download
disp('Unzipping accelerometer data analysis file...');
system('unzip ./OSF/ZeroG_Analysis.zip -d ./OSF');
system('rm ./OSF/ZeroG_Analysis.zip');
system('mkdir ./analysis');
system('mv ./OSF/analysis ./analysis/acceleration');

% Unzip flight vibration and sequencing reads files, available from
% https://osf.io/n6krq/

% Unzip Flight Sequencing Reads
disp('Unzipping Flight MinION data...');
mkdir('./data/MinION/Flight');
system('cat ./OSF/MinION/Flight/Flight.zip.* > ./data/MinION/Flight/Flight.zip');
system('unzip ./data/MinION/Flight/Flight.zip -d ./data/MinION/Flight');
system('rm ./data/MinION/Flight/Flight.zip');

% Unzip Ground Sequencing Reads
disp('Unzipping Ground MinION data...');
mkdir('./data/MinION/Ground');
system('cat ./OSF/MinION/Ground/Ground.zip.* > ./data/MinION/Ground/Ground.zip');
system('unzip ./data/MinION/Ground/Ground.zip -d ./data/MinION/Ground');
system('rm ./data/MinION/Ground/Ground.zip');

% Unzip Flight Slamstick Data
disp('Unzipping Flight Slamstick data...');
mkdir('./data/SlamStick/Flight');
system('unzip ./OSF/SlamStick/SlamStick_Flight.zip -d ./data/SlamStick');
system('rm -rf ./data/SlamStick/__MACOSX');

% Unzip Ground Slamstick Data
disp('Unzipping Ground Slamstick data...');
mkdir('./data/SlamStick/Ground');
system('unzip ./OSF/SlamStick/SlamStick_Ground.zip -d ./data/SlamStick');
system('rm -rf ./data/SlamStick/__MACOSX -rf');
