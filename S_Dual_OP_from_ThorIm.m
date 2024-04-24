%% Script, Optical Palpation via ThorImage, needs compressed and uncompressed scan
%% Run Init_Config after making relevant changes
%run(Init_Config)
%% Get OP images from a pair of ThorImages OCT scans
Im1 = loadSSOCT(FilePath_Uncompressed,ImSize);
Im2 = loadSSOCT(FilePath_Compressed,ImSize);

%% Initialise variables and stress function
x_min = 0;
x_max = ImSize(1);
y_min = 0;
y_max = ImSize(2);
%%
% initialise 2-d .mat files
layeridx_xy_1              = zeros([x_max-x_min, y_max-y_min]); % Uncompressed
layeridx_xy_2              = zeros([x_max-x_min, y_max-y_min]); % Compressed
layer_delta_xy             = zeros([x_max-x_min, y_max-y_min]);
layer_strain_bulk_xy       = zeros([x_max-x_min, y_max-y_min]);
Stress                     = zeros([x_max-x_min, y_max-y_min]);

% define strain to stress transform
% stresslayer_curve_filename = '20180508_P7676_1-1_stress_strain_0_100.mat'; % set in config
load(stresslayer_curve_filename);
Stress_Func = @(bulk_strain) interp1(layer_stress_strain.strain, ...
    layer_stress_strain.stress_kpa, ...
    bulk_strain, 'pchip');

px = ImSize(1);
py = ImSize(2);

%% Edge detection Debug
for n = 100:100 % DEBUG
    OCT_Bscan_zx = squeeze(Im1(:,n,:))';
    if config.RemoveLine == 1
        OCT_Bscan_zx(key1) = OCT_Bscan_zx(key2);
    end
    if config.UseSigmoidFit == 1
        layeridx_xy_1(:,n) = sigmoid_fit2(OCT_Bscan_zx, config, 0);
    else
        layeridx_xy_1(:,n) = oct_find_surface(OCT_Bscan_zx, config);
    end
    OCT_Bscan_zx = squeeze(Im2(:,n,:))';
    if config.RemoveLine == 1
        OCT_Bscan_zx(key1) = OCT_Bscan_zx(key2);
    end
    if config.UseSigmoidFit == 1
        layeridx_xy_2(:,n) = sigmoid_fit2(OCT_Bscan_zx, config, 1);
    else
        layeridx_xy_2(:,n) = oct_find_surface(OCT_Bscan_zx, config);
    end
    layer_delta_xy(:,n) = layeridx_xy_2(:,n) - layeridx_xy_1(:,n);
end
%% Edge detection and OP
tic 
parfor n = 1:py 
    OCT_Bscan_zx = squeeze(Im1(:,n,:))';
    if config.RemoveLine == 1
        OCT_Bscan_zx(key1) = OCT_Bscan_zx(key2);
    end
    if config.UseSigmoidFit == 1
        layeridx_xy_1(:,n) = sigmoid_fit2(OCT_Bscan_zx, config, 0);
    else
        layeridx_xy_1(:,n) = oct_find_surface(OCT_Bscan_zx, config);
    end
    OCT_Bscan_zx = squeeze(Im2(:,n,:))';
    if config.RemoveLine == 1
        OCT_Bscan_zx(key1) = OCT_Bscan_zx(key2);
    end
    if config.UseSigmoidFit == 1
        layeridx_xy_2(:,n) = sigmoid_fit2(OCT_Bscan_zx, config, 1);
    else
        layeridx_xy_2(:,n) = oct_find_surface(OCT_Bscan_zx, config);
    end
    layer_delta_xy(:,n) = layeridx_xy_2(:,n) - layeridx_xy_1(:,n);
    layer_strain_bulk_xy(:,n) = layer_delta_xy(:,n) ./ (LayerThickness/config.SPACING_mm_pix(3));
    Stress(:,n) = Stress_Func(-layer_strain_bulk_xy(:,n));
end
toc
%% Debugged individual lines?
n = 960;
layer_delta_xy(:,n) = layeridx_xy_2(:,n) - layeridx_xy_1(:,n);
layer_strain_bulk_xy(:,n) = layer_delta_xy(:,n) ./ (LayerThickness/config.SPACING_mm_pix(3));
Stress(:,n) = Stress_Func(-layer_strain_bulk_xy(:,n));

%% True Linear Layer Material Stress Map
Stress_lin = -log(1+layer_strain_bulk_xy)*E_lin/1000.0;

%% Save Other
save([FilePath_Compressed(1:(end-14)),'OtherSig2.mat'],'layer_delta_xy','layer_strain_bulk_xy','layeridx_xy_2','layeridx_xy_1');

%%
Stress = Stress_Func(-layer_strain_bulk_xy);
%% Make Stress Figure
% createfigure(Stress')
S1 = Stress;
%T1 = (Stress.*(1+layer_strain_bulk_xy));
T1 = Stress_lin;
y_int = 1:y_max;
for n2 = 1:x_max
    for n1 = 1:y_max
        if isnan(S1(n1,n2))
            temp = S1(:,n2);
            temp(isnan(temp)) = interp1(y_int(~isnan(temp)),temp(~isnan(temp)),y_int(isnan(temp)),'linear','extrap'); %, 'spline');
            S1(:,n2) = temp;
        end
        if isnan(T1(n1,n2))
            temp = T1(:,n2);
            temp(isnan(temp)) = interp1(y_int(~isnan(temp)),temp(~isnan(temp)),y_int(isnan(temp)),'linear','extrap'); %, 'spline');
            T1(:,n2) = temp;
        end
    end
end
% createfigure(S1')
createfigure(T1')

% T1 = medfilt2(T1,[1 5],'symmetric');
% 
% createfigure(T1')

disp(mean2(T1))
disp(std2(T1))
StdNum = 4.5; % adjust this to only remove artefacts
RepIm = imgaussfilt(T1,[100 100]);
RepImS = imgaussfilt(S1,[100 100]);
logical = T1 > mean2(T1)+StdNum*std2(T1) |  T1 < mean2(T1)-StdNum*std2(T1);
logicalS = S1 > mean2(S1)+StdNum*std2(S1) |  S1 < mean2(S1)-StdNum*std2(S1);
T1(logical) = RepIm(logical);
S1(logicalS) = RepImS(logicalS);

createfigure(T1')

%S2 = medfilt2(S1,[1 3],'symmetric');
%T2 = medfilt2(T1,[1 3],'symmetric');
S2 = medfilt2(S1,[3 3],'symmetric');
T2 = medfilt2(T1,[3 3],'symmetric');
%createfigure(S2')
% createfigure(T2')
S3 = imgaussfilt(S2,[2 2]);
T3 = imgaussfilt(T2,[2 2]);
%createfigure(S3')
createfigure(T3')
T3 = imgaussfilt(T3,[4 4]);
createfigure(T3')
%%
save([FilePath_Compressed(1:(end-14)),'StressImageSig2.mat'],'T3','T2','T1','S3','S2','S1','Stress');
%% Smoothing Displacements for COP
Lp = LayerThickness+layer_delta_xy*config.SPACING_mm_pix(3);
%createfigure(Lp')
Lp1 = Lp;
y_int = 1:y_max;
for n2 = 1:x_max
    for n1 = 1:y_max
        if isnan(Lp1(n1,n2))
            temp = Lp1(:,n2);
            temp(isnan(temp)) = interp1(y_int(~isnan(temp)),temp(~isnan(temp)),y_int(isnan(temp)),'linear','extrap'); %, 'spline');
            Lp1(:,n2) = temp;
        end
    end
end
%createfigure(Lp1')

createfigure(Lp1')
StdNum = 4.0;
RepIm = imgaussfilt(Lp1,[100 100]);
logical = Lp1 > mean2(Lp1)+StdNum*std2(Lp1) |  Lp1 < mean2(Lp1)-StdNum*std2(Lp1);
Lp1(logical) = RepIm(logical);
createfigure(Lp1')

Lp2 = medfilt2(Lp1,[3 3],'symmetric');
%Lp2 = medfilt2(Lp1,[1 3],'symmetric');
%createfigure(Lp2')
Lp3 = imgaussfilt(Lp2,[3 3]);
createfigure(Lp3')
%%
save([FilePath_Compressed(1:(end-14)),'LTImageSig2.mat'],'Lp3','Lp2','Lp1','Lp');

%% Make Strain Figure
% createfigure(layer_strain_bulk_xy')
NE = layer_strain_bulk_xy;
y_int = 1:y_max;
for n2 = 1:x_max
    for n1 = 1:y_max
        if isnan(NE(n1,n2))
            temp = NE(:,n2);
            temp(isnan(temp)) = interp1(y_int(~isnan(temp)),temp(~isnan(temp)),y_int(isnan(temp)),'linear','extrap'); %, 'spline');
            NE(:,n2) = temp;
        end
    end
end

%createfigure(NE')
disp(mean2(NE))
disp(std2(NE))
StdNum = 2.0;
RepIm = medfilt2(NE,[20 20],'symmetric');
logical = NE > mean2(NE)+StdNum*std2(NE) |  NE < mean2(NE)-StdNum*std2(NE);
NE(logical) = RepIm(logical);

NE2 = medfilt2(NE,[3 3],'symmetric');
NE3 = imgaussfilt(NE2,[2 2]);

createfigure(NE3')

%% Make Uncompressed Layer Figure
% createfigure(layer_strain_bulk_xy')
L = layeridx_xy_1*config.SPACING_mm_pix(3);
y_int = 1:y_max;
for n2 = 1:x_max
    for n1 = 1:y_max
        if isnan(L(n1,n2))
            temp = L(:,n2);
            temp(isnan(temp)) = interp1(y_int(~isnan(temp)),temp(~isnan(temp)),y_int(isnan(temp)),'linear','extrap'); %, 'spline');
            L(:,n2) = temp;
        end
    end
end

%createfigure(NE')
disp(mean2(L))
disp(std2(L))
StdNum = 3.0;
createfigure(L')
%RepIm = medfilt2(L,[20 20],'symmetric');
RepIm = imgaussfilt(L,[100 100]);
logical = L > mean2(L)+StdNum*std2(L) |  L < mean2(L)-StdNum*std2(L);
L(logical) = RepIm(logical);

createfigure(L')

L2 = medfilt2(L,[3 3],'symmetric');
L3 = imgaussfilt(L2,[2 2]);

createfigure(L3')

%% Extending edges for COP (to reduce corner effects)
Lpe = LayerThickness+layer_delta_xy*config.SPACING_mm_pix(3);
%createfigure(Lpe')

Lpe1 = Lpe;
y_int = 1:y_max;
for n2 = 1:x_max
    for n1 = 1:y_max
        if isnan(Lpe1(n1,n2))
            temp = Lpe1(:,n2);
            temp(isnan(temp)) = interp1(y_int(~isnan(temp)),temp(~isnan(temp)),y_int(isnan(temp)),'linear','extrap'); %, 'spline');
            Lpe1(:,n2) = temp;
        end
    end
end
temp1 = Lpe1;
padpixels = 100;
padval = (mean(mean(Lpe1(1:padpixels,1:end)))+mean(mean(Lpe1(x_max-padpixels:end,1:end)))+mean(mean(Lpe1(1:end,1:padpixels)))+mean(mean(Lpe1(1:end,y_max-padpixels))))/4;
Lpe1 = padarray(Lpe1,[padpixels padpixels],padval,'both');
%createfigure(Lpe1')
temp = imgaussfilt(Lpe1,[50 50]);
temp(padpixels+1:x_max+padpixels, padpixels+1:y_max+padpixels) = temp1;
Lpe2 = temp;
%createfigure(Lpe2')

%createfigure(Lpe2')
StdNum = 4.7;                  % MIGHT NEED TO INCREASE DEPENDING ON DETECTABILITY, LOOK FOR TOO MUCH SMOOTHING OVER INCLUSION 
RepIm = imgaussfilt(Lpe2,[100 100]);
logical = Lpe2 > mean2(Lpe2)+StdNum*std2(Lpe2) |  Lpe2 < mean2(Lpe2)-StdNum*std2(Lpe2);
Lpe2(logical) = RepIm(logical);
%createfigure(Lpe2')

% all this is specifically for 5mm COP, COP functions best if non-physical imaging artefacts are removed (check OCT scans to see if edge detection error or physical)
logical = false([n1+2*padpixels n2+2*padpixels]);
%logical(670:720,1030:1070) = true;
%logical(650:675,987:1005) = true;
%logical(140:200,1240:1270) = true;
%logical(1000:1030,1244:1270) = true;
%logical(1120:1180,280:360) = true;
%logical(945:970,674:694) = true;
%logical(832:843,778:781) = true;
%logical(882:892,699:700) = true;
%logical(130:180,270:325) = true;
%logical(140:180,100:130) = true;
%logical(156:172,376:390) = true;
%logical(905:926,290:340) = true;
%logical(420:470,280:335) = true;
%logical(450:490,155:200) = true;
%logical(945:970,670:696) = true;
%logical(950:970,800:820) = true;
%logical(1020:1045,420:455) = true;
%logical(765:790,1002:1024) = true;
%logical(950:970,1000:1017) = true;
%logical(488:502,826:825) = true;
%logical(500:512,887:895) = true;
%logical(506:522,904:916) = true;
%logical(910:925,525:545) = true;
%logical(1230:1270,914:928) = true;
Lpe2(logical) = RepIm(logical);
createfigure(logical')
Lpe2(logical) = nan;

RepImEdge = imgaussfilt(Lpe2,[10 10]);
logical = false([n1+2*padpixels n2+2*padpixels]);
%logical(120:180,275:325) = true;
edge = 50;
logical(padpixels:n1+padpixels,padpixels-edge:edge+padpixels) = true;
logical(padpixels:n1+padpixels,n2+padpixels-edge:edge+n2+padpixels) = true;
logical(padpixels-edge:edge+padpixels,padpixels:n2+padpixels) = true;
logical(n1+padpixels-edge:edge+n1+padpixels,padpixels:n2+padpixels) = true;
%logical(830:850,772:784) = true;
Lpe2(logical) = RepImEdge(logical);
%createfigure(Lpe2')
%Lpe2(logical) = nan;

y_int = 1:y_max+2*padpixels;
for n2 = 1:x_max+2*padpixels
    for n1 = 1:y_max+2*padpixels
        if isnan(Lpe2(n1,n2))
            temp = Lpe2(:,n2);
            temp(isnan(temp)) = interp1(y_int(~isnan(temp)),temp(~isnan(temp)),y_int(isnan(temp)),'linear','extrap'); %, 'spline');
            Lpe2(:,n2) = temp;
        end
    end
end

Lpe3 = medfilt2(Lpe2,[3 3],'symmetric');
%Lpe3 = medfilt2(Lpe2,[1 3],'symmetric');
%createfigure(Lpe3')
%Lpe4 = imgaussfilt(Lpe3,[3 3]);
%createfigure(Lpe4')
%Lpe4 = medfilt2(Lpe4,[5 5],'symmetric');
Lpe4 = imgaussfilt(Lpe3,[3 3]);
createfigure(Lpe4')
%%
save([FilePath_Compressed(1:(end-14)),'LTeImageSig7.mat'],'Lpe4','Lpe3','Lpe2','Lpe1','Lpe','padpixels');


%% Generate inputs to COP
% Meta data
meta.job = '231104_Inc_5_5mm3_10p_highres_sig7';
meta.model = 'Model-1';
meta.gen = 'Rhys';
meta.output_path =  'F:\Rhys\COP_Outputs';
meta.abaqus_work_dir = 'F:\Rhys\AbaqusModels\COP_Workspace';
% Thickness data
L0 = LayerThickness;
load([FilePath_Compressed(1:(end-14)),'LTeImageSig7.mat']);
% Material
% UNCOMMENT BELOW FOR HYPERELASTIC, COMMENT BELOW THAT
% clear material
% material.C10 = 0.002402;  % 0.0022;
% material.C01 = 4.1371e-04; % 7.03E-4;
% material.Den = 0.01; 
% IF WANT LINEAR ELASTIC COMMENT OUT ABOVE, UNCOMMENT BELOW
clear material
material.E = E_lin; %16320.0; 
material.v = 0.48; 
material.Den = 0.005; % need to fix rf to layer density ratio to reflect mesh size ratio
% Meshing THIS TEXT IS HERE TO MAKE SEEING THINGS EASIER
opt.size = [config.FOV_mm(1)+padpixels*config.SPACING_mm_pix(1),config.FOV_mm(2)+padpixels*config.SPACING_mm_pix(2)]; % FOV (mm)
opt.Padpixels_mm = (padpixels*config.SPACING_mm_pix(1));
opt.finalPix = [config.ImSize(1),config.ImSize(2)];
opt.mesh.mesh_type = 'hex'; %{‘hex’,’tet’,'tet10'}
opt.mesh.elem_type = 'C3D8';
opt.mesh.mesh_size = 0.1; % 0.5,0.25,0.1,   0.05 % specifies mesh size, isotropic if length=1, in [x,y,z] if length=3;
% Friction
opt.frict.plate = 0.01; %friction coefficient of layer with plate, set to 0.01 if not written on job name
opt.frict.sample = 0.000001; %friction coefficient of layer with sample (should be small as sample is supposed to expand with layer)
% Interaction Properties, 1.0 makes the surface in question the master surface with respect to the layer
opt.int.weight.plate = 1.0; %0.9, 0.5
opt.int.weight.sample = 1.0; % 0.75, 0.5

%% Run COP
save([meta.output_path,'\','inputs_',meta.job,'.mat'],'meta','L0','Lpe4','material','opt');
[T] = main_cop_explicit( meta, L0, Lpe4, material, opt ); % Only have one in path

%% Read COP
meta.output_path =  'F:\Rhys\COP_Outputs';
meta.job = '231104_Inc_5_5mm3_10p_highres_sig7';
load([meta.output_path,'\','inputs_',meta.job,'.mat'],'meta','L0','Lpe4','material','opt');
[T] = read_cop_explicit( meta, L0, Lpe4, material, opt );
%%
T = T/1000;
%% Visualise Results
%figure;imagesc(T);colormap(gray(256));
createfigure(-T/1000)
loc = [FilePath_Compressed(1:(end-14)),'\',meta.job];
mkdir(loc);
save([loc,'\','COP_T.mat'],'T');
save([loc,'\','inputs_',meta.job,'.mat'],'meta','L0','Lpe4','material','opt');
