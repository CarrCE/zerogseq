% zerogseq_00_download.m

% Downloading all data... this may take 20-25 minutes on a fast connection.

% Download all data required for analysis
mkdir('./OSF');

disp('Downloading accelerometer g-level results');
websave('./OSF/ZeroG_Analysis.zip','https://osf.io/pmhj4/download');

disp('Downloading vibration data (please be patient)');
mkdir('./OSF/Slamstick');
websave('./OSF/Slamstick/Slamstick_Flight.zip','https://osf.io/cxwj7/?action=download');
websave('./OSF/Slamstick/Slamstick_Ground.zip','https://osf.io/wykdb/?action=download');

disp('Downloading sequencing data (please be *very* patient)');
mkdir('./OSF/MinION/Flight');
websave('./OSF/MinION/Flight/Flight.txt','https://osf.io/erp5h/?action=download');
websave('./OSF/MinION/Flight/Flight.zip.001','https://osf.io/8swcz/?action=download');
websave('./OSF/MinION/Flight/Flight.zip.002','https://osf.io/sn8zc/?action=download');
websave('./OSF/MinION/Flight/Flight.zip.003','https://osf.io/yk7qv/?action=download');
websave('./OSF/MinION/Flight/Flight.zip.004','https://osf.io/a59ub/?action=download');
websave('./OSF/MinION/Flight/Flight.zip.005','https://osf.io/db5yx/?action=download');
websave('./OSF/MinION/Flight/Flight.zip.006','https://osf.io/zjhsb/?action=download');
websave('./OSF/MinION/Flight/Flight.zip.007','https://osf.io/bxa8u/?action=download');
websave('./OSF/MinION/Flight/Flight.zip.008','https://osf.io/b6fvc/?action=download');
websave('./OSF/MinION/Flight/Flight.zip.009','https://osf.io/6vjre/?action=download');
websave('./OSF/MinION/Flight/Flight.zip.010','https://osf.io/hrv96/?action=download');

mkdir('./OSF/MinION/Ground');
websave('./OSF/MinION/Ground/Ground.txt','https://osf.io/cdeys/?action=download');
websave('./OSF/MinION/Ground/Ground.zip.001','https://osf.io/57qwf/?action=download');
websave('./OSF/MinION/Ground/Ground.zip.002','https://osf.io/pq82k/?action=download');
websave('./OSF/MinION/Ground/Ground.zip.003','https://osf.io/45dtk/?action=download');
