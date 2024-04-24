% For specific file
FilePath_Uncompressed = 'F:\Rhys\MATLAB Workspace\231104_Inc_5_5mm3\231104_Inc_5_5mm3_0p_L\data\Intensity.data';
FilePath_Compressed = 'F:\Rhys\MATLAB Workspace\231104_Inc_5_5mm3\231104_Inc_5_5mm3_10p_L\data\Intensity.data';

stresslayer_curve_filename = '20201214_SILPURAN2400_1to1to2_proc_stress_strain_0_100.mat';

ImSize = [1240,1240,1024];

LayerThickness = 5.8;

config.UseSigmoidFit = 1; % set to 1 to use sigmoid fitting for edge detection, 
%only applicable when sample has large OCT contrast conpared to layer
config.HistFilt = 3; % imgaussfilt(N,config.HistFilt)
config.NoHist = 0; % use near bottom of B-Scan
config.NoHistW = 50;
config.InterfacePeakM = 0.1;
config.InterfacePeakW = 20;
config.InterFlatM = 0.0; % use near bottom of B-Scan
config.InterFlatW = 20;
config.MinHistPeakSep = 3.0; % Choose half the average intensity difference between the layer and the top of the sample

E_lin = 16210.6; %16745.59; % 9.3% strain, approximate linear elasticity of the layer

config.RemoveLine = 0; % 1 if want to remove line that is the same in every scan

config.ImSize = ImSize;

config.FOV_mm = [15.5,15.5];

config.SPACING_mm_pix = [15.5/1240,15.5/1240,0.0025];

config.EXTENT_Z = [350, 750];  % All
config.EXTENT_Z_u = [500, 750];  % 0%
config.EXTENT_Z_c = [260, 510];  % 10%
% config.EXTENT_Z_c = [50, 400];  % 20%

config.OCT_NOMINAL_DYNAMIC_RANGE = 60; 

config.OCT_MIN_DB_RANGE = [0, 15]; 

config.OCT_SIGMA = 3; 

config.EDGE_SIGMA = 3; 

config.DETECT_NEG_GRAD = -1; % -1 for oil-sample, 1 for layer(non-scattering)-oil

config.PERCENT_PIX_NOT_EDGES = 0.95; 

config.LAYER_OFFSET_Z = 0;