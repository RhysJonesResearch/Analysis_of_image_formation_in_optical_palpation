% For specific file
FilePath_Uncompressed = 'F:\Rhys\MATLAB Workspace\231104_Inc_5_0p5mm3\231104_Inc_5_0p5mm3_0p_O_R\data\Intensity.data';
FilePath_Compressed = 'F:\Rhys\MATLAB Workspace\231104_Inc_5_0p5mm3\231104_Inc_5_0p5mm3_10p_L_R_O\data\Intensity.data';

stresslayer_curve_filename = '20201214_SILPURAN2400_1to1to2_proc_stress_strain_0_100.mat';

ImSize = [1240,1240,1024];

LayerThickness = 1.098; %1.21; %1.063;

config.UseSigmoidFit = 1; % set to 1 to use sigmoid fitting for edge detection (1 or 0)
%only applicable when sample has large OCT contrast conpared to layer
config.HistFilt = 4; % imgaussfilt(N,config.HistFilt)
config.MinHistPeakSep = 4.0; % Choose half the average intensity difference between the layer and the top of the sample
config.NoHist = 0; % use near bottom of B-Scan (1 or 0)
config.NoHistW = 50;
config.InterfacePeakM = 0.1;
config.InterfacePeakW = 10;
config.InterFlatM = 0.0; % use near bottom of B-Scan
config.InterFlatW = 20;

E_lin = 16410.77; % 10.81% N_strain, approximate linear elasticity of the layer

config.RemoveLine = 1; % 1 if want to remove line that is the same in every scan
%     img(img<36) = 0;
%     img(:,1:300) = 0;
%     img(:,400:end) = 0;

config.ImSize = ImSize;

config.FOV_mm = [15.5,15.5];

config.SPACING_mm_pix = [15.5/1240,15.5/1240,2.56/1024];

config.EXTENT_Z = [350, 750];  % All
config.EXTENT_Z_u = [380, 500];  % 0%
config.EXTENT_Z_c = [250, 500];  % 10%
% config.EXTENT_Z_c = [50, 400];  % 20%

config.OCT_NOMINAL_DYNAMIC_RANGE = 60; 

config.OCT_MIN_DB_RANGE = [0, 15]; 

config.OCT_SIGMA = 3; 

config.EDGE_SIGMA = 3; 

config.DETECT_NEG_GRAD = -1; % -1 for oil-sample, 1 for layer(non-scattering)-oil

config.PERCENT_PIX_NOT_EDGES = 0.95; 

config.LAYER_OFFSET_Z = 0;