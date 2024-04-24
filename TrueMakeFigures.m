% Make Nice Figures
% FilePath_Uncompressed = 'F:\Rhys\MATLAB Workspace\230907_Inc_4_2_5mm2\230907_Inc_4_2_5mm2_0p\Intensity.data';
% FilePath_Compressed = 'F:\Rhys\MATLAB Workspace\230907_Inc_4_2_5mm2\230907_Inc_4_2_5mm2_10p\Intensity.data';

%% Stress Map
% load([FilePath_Compressed(1:(end-14)),'StressImage.mat'], 'S3');
cdata1 = S3';
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
image(cdata1,'Parent',axes1,'CDataMapping','scaled');
S = size(cdata1);
xlim(axes1,[0.5 S(1)+0.5]);
ylim(axes1,[0.5 S(2)+0.5]);
title('Nominal Stress (kPa) - L5_2,S4_2,P10,ROI15p5');
box(axes1,'on');
axis(axes1,'ij');
set(axes1,'Layer','top');
colorbar('peer',axes1);

%% True Stress Map with circle for real inclusion size
cT = rot90(0.825*(T3),2);
%cT = rot90(-T/1000,0);
centre = [632,600];  %[624,632];                    % NORMALLY DON'T USE THIS
inner_det_rad = 1240/15.5*1.9; % pixels radius for centre point
output = rotateAround(cT, centre(1), centre(2), 0);
S = size(output);
cdata1 = output; % '
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
image(cdata1,'Parent',axes1,'CDataMapping','scaled');
xlim(axes1,[0.5 S(1)+0.5]);
ylim(axes1,[0.5 S(2)+0.5]);
title('L: 5mm#2, S: Inc 4, P: 10%, 15.5mmx15.5mm, kPa');
box(axes1,'on');
axis(axes1,'ij');
set(axes1,'Layer','top');
colorbar('peer',axes1);
caxis([0,3.55]);
x = centre(1) - inner_det_rad; y = centre(2) - inner_det_rad;
d = inner_det_rad * 2;
pos = [x y d d];
%rectangle('Position',pos,'Curvature',[1 1])

%% True Stress Map
%load([FilePath_Compressed(1:(end-14)),'StressImage.mat'], 'T3');
cdata1 = rot90(0.825*(T3'),4);
%cdata1 = rot90((T3'),2);
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
image(cdata1,'Parent',axes1,'CDataMapping','scaled');
S = size(cdata1);
xlim(axes1,[0.5 S(1)+0.5]);
ylim(axes1,[0.5 S(2)+0.5]);
caxis([0,3.6]);
title('5');
%title('True Stress (kPa) - L5_2,S4_2,P10,ROI15p5');
box(axes1,'on');
axis(axes1,'ij');
set(axes1,'Layer','top');
colorbar('peer',axes1);

%% Layer Thickness
% load([FilePath_Compressed(1:(end-14)),'LTImage.mat'], 'Lp3');
cdata1 = Lpe4';
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
image(cdata1,'Parent',axes1,'CDataMapping','scaled');
S = size(cdata1);
xlim(axes1,[0.5 S(1)+0.5]);
ylim(axes1,[0.5 S(2)+0.5]);
title('L: 2mm#2, S: Inc 2, P: 10%, 15mmx15mm, mm');
box(axes1,'on');
axis(axes1,'ij');
set(axes1,'Layer','top');
colorbar('peer',axes1);
%% Resolution find centre OP
% load([FilePath_Compressed(1:(end-14)),'StressImage.mat'], 'T3');
T3 = T3;

Square_TL = 450; % Want to choose square that is slightly inside the inclusion circle, for gaussian fitting to work well
Square_BR = 750;
MPW = 0.8; % min peak width (mm)
MPP = 10; %0.0015; % min peak prominence
centre = [0,0];
widths_row = zeros(1,Square_BR-Square_TL);
widths_col = zeros(1,Square_BR-Square_TL);
h = [1/2 1/2];
binomialCoeff = conv(h,h);
n_max = 2000;
for n = 1:n_max
    binomialCoeff = conv(binomialCoeff,h);
end
fDelay = (length(binomialCoeff)-1)/2;
for row = Square_TL:Square_BR
    S0 = size(T3(:,row));
    cdata1 = [T3(:,row); zeros([1000,1])];
    binomialMA = filter(binomialCoeff, 1, cdata1);
    cdata1 = diff(binomialMA);
    S1 = size(cdata1);
    x = (((fDelay+ceil(0.1*S0):S1-ceil(0.1*S0)))-fDelay)*config.SPACING_mm_pix(1);
    y = abs(cdata1(fDelay+ceil(0.1*S0):S1-ceil(0.1*S0)));
    [~,locs,~,~] = findpeaks(y,x,'SortStr','descend','MinPeakWidth',MPW,'MinPeakProminence',MPP);
    if size(locs,2)>1
            width = abs(locs(1)-locs(2));
            widths_row(row-Square_TL+1) = width;
            centre(1) = row;
    elseif row ~= Square_TL
            widths_row(row-Square_TL+1) = 0; % widths_row(row-Square_TL);
    end
    if row == 630
        figure
        findpeaks(y,x,'Annotate','extents')
    end
end
S4 = T3';
for col = Square_TL:Square_BR
    S0 = size(S4(:,col));
    cdata1 = [S4(:,col); zeros([1000,1])];
    binomialMA = filter(binomialCoeff, 1, cdata1);
    cdata1 = diff(binomialMA);
    S1 = size(cdata1);
    x = (((fDelay+ceil(0.1*S0):S1-ceil(0.1*S0)))-fDelay)*config.SPACING_mm_pix(1);
    y = abs(cdata1(fDelay+ceil(0.1*S0):S1-ceil(0.1*S0)));
    [~,locs,~,~] = findpeaks(y,x,'SortStr','descend','MinPeakWidth',MPW,'MinPeakProminence',MPP);
    if size(locs,2)>1
            width = abs(locs(1)-locs(2));
            widths_col(col-Square_TL+1) = width;
            centre(2) = col;
    elseif col ~= Square_TL
            widths_col(col-Square_TL+1) = 0; % widths_col(row-Square_TL);
    end
%     if col == 471
%         figure
%         findpeaks(y,x,'Annotate','extents')
%     end
end
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
startPoints = [mean(widths_row) (Square_BR-Square_TL)/2 (Square_BR-Square_TL)/2 0];
LowBounds = [0, Square_TL, 0, 0];
UpBounds = [max(widths_row)*2, Square_BR, Square_BR-Square_TL, 0];
f_row = fit((Square_TL:Square_BR)',(widths_row)',gaussEqn,'Start', startPoints,'Lower', LowBounds,'Upper', UpBounds); %, 'Exclude', x < 800);
f_col = fit((Square_TL:Square_BR)',(widths_col)',gaussEqn,'Start', startPoints,'Lower', LowBounds,'Upper', UpBounds); %, 'Exclude', x < 800);
row_cs = coeffvalues(f_row);
col_cs = coeffvalues(f_col);
centre = [ceil(row_cs(2)) ceil(col_cs(2))];
figure
plot(f_row,(Square_TL:Square_BR)',(widths_row)')
figure
plot(f_col,(Square_TL:Square_BR)',(widths_col)')
disp(centre)
disp('reverse order for next section')

%% Stress Map - Choose ROI's
T3 = T;
centre = [625,624];                    % NORMALLY DON'T USE THIS
Square_TL = 450;
Square_BR = 810;
h = [1/2 1/2];
binomialCoeff = conv(h,h);
n_max = 2000;% 2000
for n = 1:n_max
    binomialCoeff = conv(binomialCoeff,h);
end
fDelay = (length(binomialCoeff)-1)/2;
inner_det_rad = 30; % pixels radius for centre point
outer_det_min_r = 450; % radial range for far point
outer_det_max_r = 550; % radial range for far point
output = rotateAround(T3, centre(1), centre(2), 0);
S = size(output);
cdata1 = output; % '
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
image(cdata1,'Parent',axes1,'CDataMapping','scaled');
xlim(axes1,[0.5 S(1)+0.5]);
ylim(axes1,[0.5 S(2)+0.5]);
title('L: 5mm#2, S: Inc 4, P: 10%, 15.5mmx15.5mm, kPa');
box(axes1,'on');
axis(axes1,'ij');
set(axes1,'Layer','top');
colorbar('peer',axes1);
x = centre(1) - inner_det_rad; y = centre(2) - inner_det_rad;
d = inner_det_rad * 2;
pos = [x y d d];
rectangle('Position',pos,'Curvature',[1 1])
x = centre(1) - outer_det_min_r; y = centre(2) - outer_det_min_r;
d = outer_det_min_r * 2;
pos = [x y d d];
rectangle('Position',pos,'Curvature',[1 1])
x = centre(1) - outer_det_max_r; y = centre(2) - outer_det_max_r;
d = outer_det_max_r * 2;
pos = [x y d d];
rectangle('Position',pos,'Curvature',[1 1])

%% Resolution - full circle OP - Just centre line
% load([FilePath_Compressed(1:(end-14)),'LTImage.mat'], 'Lp3');
% load([FilePath_Compressed(1:(end-14)),'Other.mat'], 'layer_strain_bulk_xy');
LSB = layer_strain_bulk_xy;
MPW = 0.7; % min peak width (mm)
MPP = 8; %8;%0.0012; % min peak prominence
nmax1 = 500; %500;
res_1 = zeros(1,nmax1);
res_2 = zeros(1,nmax1);
S_centre_1 = zeros(1,nmax1);
S_centre_2 = zeros(1,nmax1);
S_far_1 = zeros(1,nmax1);
S_far_2 = zeros(1,nmax1);
L_centre_1 = zeros(1,nmax1);
L_centre_2 = zeros(1,nmax1);
L_far_1 = zeros(1,nmax1);
L_far_2 = zeros(1,nmax1);
LSB_far_1 = zeros(1,nmax1);
LSB_far_2 = zeros(1,nmax1);
smoothing = 0;
if smoothing == 1
    %T4 = medfilt2(T3,[10 10]); % T3
    T4 = imgaussfilt(T3,[5 5]);
    %Lp4 = medfilt2(Lp3,[10 10]);
    %Lp4 = imgaussfilt(Lp4,[10 10]);
    %LSB = medfilt2(LSB,[10 10]);
    %LSB = imgaussfilt(LSB,[10 10]);
else
    T4 = T3; %T3
    Lp4 = Lp3;
    LSB = layer_strain_bulk_xy;
end

h = [1/2 1/2];
binomialCoeff = conv(h,h);
n_max = 1750;% 2000
for n = 1:n_max
    binomialCoeff = conv(binomialCoeff,h);
end
fDelay = (length(binomialCoeff)-1)/2;

parfor n1 = 1:nmax1
    output = rotateAround(T4, centre(1), centre(2), n1*180/nmax1);
    outputL = rotateAround(Lp4, centre(1), centre(2), n1*180/nmax1);
    outputLSB = rotateAround(LSB, centre(1), centre(2), n1*180/nmax1);
    
    S_centre_1(n1) = mean(output(centre(1)-inner_det_rad:centre(1),centre(2)));
    S_centre_2(n1) = mean(output(centre(1):centre(1)+inner_det_rad,centre(2)));
    S_far_1(n1) = mean(output(centre(1)-outer_det_max_r:centre(1)-outer_det_min_r,centre(2)));
    S_far_2(n1) = mean(output(centre(1)+outer_det_min_r:centre(1)+outer_det_max_r,centre(2)));
    L_centre_1(n1) = mean(outputL(centre(1)-inner_det_rad:centre(1),centre(2)));
    L_centre_2(n1) = mean(outputL(centre(1):centre(1)+inner_det_rad,centre(2)));
    L_far_1(n1) = mean(outputL(centre(1)-outer_det_max_r:centre(1)-outer_det_min_r,centre(2)));
    L_far_2(n1) = mean(outputL(centre(1)+outer_det_min_r:centre(1)+outer_det_max_r,centre(2)));
    LSB_far_1(n1) = mean(outputLSB(centre(1)-outer_det_max_r:centre(1)-outer_det_min_r,centre(2)));
    LSB_far_2(n1) = mean(outputLSB(centre(1)+outer_det_min_r:centre(1)+outer_det_max_r,centre(2)));
    
    flag = 0;
    n2 = 0;
    while flag == 0 && n2 < 20; % 20
        row = ceil(centre(1)) + n2*(-1)^n2;
        S0 = size(output(row,:),2);
        cdata1 = [output(row,:), zeros([1,n_max/2])];
        binomialMA = filter(binomialCoeff, 1, cdata1);
        cdata1 = diff(binomialMA);
        S1 = size(cdata1,2);
        x = (((fDelay+ceil(0.2*S0):S1-ceil(0.2*S0)))-fDelay)*config.SPACING_mm_pix(1);
        y = abs(cdata1(fDelay+ceil(0.2*S0):S1-ceil(0.2*S0)));
%         figure
%         plot(output(row,:))
%         hold on
%         plot(binomialMA(n_max/2:end))
        [~,m] = min(y(floor(size(x,2)*2/6):floor(size(x,2)*4/6)));
        m = m+floor(size(x,2)*2/6);
        [pks_l,locs_l,widths_l,proms_l] = findpeaks(y(1:m),'SortStr','descend','MinPeakWidth',MPW/config.SPACING_mm_pix(1),'MinPeakProminence',MPP);
        [pks_r,locs_r,widths_r,proms_r] = findpeaks(y(m:end),'SortStr','descend','MinPeakWidth',MPW/config.SPACING_mm_pix(1),'MinPeakProminence',MPP);
%         figure
%         findpeaks(y,'Annotate','extents');
        if size(widths_l,2)>0 && size(widths_r,2)>0
            z = y(1:floor(locs_l(1)));
            [~,~,widths_l1,~] = findpeaks([min(z),z,min(z)],'SortStr','descend','MinPeakWidth',MPW/config.SPACING_mm_pix(1)/3,'MinPeakProminence',MPP);
            z = y(ceil(m+locs_r(1)):end);
            [~,~,widths_r1,~] = findpeaks([min(z),z,min(z)],'SortStr','descend','MinPeakWidth',MPW/config.SPACING_mm_pix(1)/3,'MinPeakProminence',MPP);
            z = y(ceil(locs_l(1)):m);
            [~,~,widths_l2,~] = findpeaks([min(z),z,min(z)],'SortStr','descend','MinPeakWidth',MPW/config.SPACING_mm_pix(1)/3,'MinPeakProminence',MPP);
            z = y(m:floor(m+locs_r(1)));
            [~,~,widths_r2,~] = findpeaks([min(z),z,min(z)],'SortStr','descend','MinPeakWidth',MPW/config.SPACING_mm_pix(1)/3,'MinPeakProminence',MPP);
%            figure
%            findpeaks(y,'Annotate','extents');
%             figure
%             findpeaks([0;y(1:floor(locs_l(1)));0],'Annotate','extents');
%             figure
%             findpeaks([0;y(ceil(m+locs_r(1)):end);0],'Annotate','extents');
%             figure
%             findpeaks([0;y(ceil(locs_l(1)):m);0],'Annotate','extents');
%             figure
%             findpeaks([0;y(m:floor(m+locs_r(1)));0],'Annotate','extents');
            if size(widths_l1,2)>0 && size(widths_r1,2)>0 && size(widths_l2,2)>0 && size(widths_r2,2)>0
                widths_l = [widths_l1(1), widths_l2(1)];
                widths_r = [widths_r1(1), widths_r2(1)];
                res_1(n1) = min(widths_l)*2; %min
                res_2(n1) = min(widths_r)*2; %min
                flag = 1;
            else
                n2 = n2 + 1;
            end
        else
            n2 = n2 + 1;
            
        end
    end
end
%% Stress Map - Choose ROI's
T3 = -T/1000;
%T3(600:601,:) = 00;
centre = [628,622];                    % NORMALLY DON'T USE THIS
Square_TL = 450;
Square_BR = 810;
h = [1/2 1/2];
binomialCoeff = conv(h,h);
n_max = 2000;% 2000
for n = 1:n_max
    binomialCoeff = conv(binomialCoeff,h);
end
fDelay = (length(binomialCoeff)-1)/2;
inner_det_rad = 30; % pixels radius for centre point
outer_det_min_r = 450; % radial range for far point
outer_det_max_r = 550; % radial range for far point
output = rotateAround(T3, centre(1), centre(2), 0*180/500);
S = size(output);
cdata1 = output; % '
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
image(cdata1,'Parent',axes1,'CDataMapping','scaled');
%xlim(axes1,[0.5 S(1)+0.5]);
%ylim(axes1,[0.5 S(2)+0.5]);
title('L: 5mm#2, S: Inc 4, P: 10%, 15.5mmx15.5mm, kPa');
box(axes1,'on');
axis(axes1,'ij');
set(axes1,'Layer','top');
colorbar('peer',axes1);
x = centre(2) - inner_det_rad; y = centre(1) - inner_det_rad;
d = inner_det_rad * 2;
pos = [x y d d];
rectangle('Position',pos,'Curvature',[1 1])
x = centre(2) - outer_det_min_r; y = centre(1) - outer_det_min_r;
d = outer_det_min_r * 2;
pos = [x y d d];
rectangle('Position',pos,'Curvature',[1 1])
x = centre(2) - outer_det_max_r; y = centre(1) - outer_det_max_r;
d = outer_det_max_r * 2;
pos = [x y d d];
rectangle('Position',pos,'Curvature',[1 1])
%% Resolution - full circle OP - Just centre line - NEWRES
% load([FilePath_Compressed(1:(end-14)),'LTImage.mat'], 'Lp3');
% load([FilePath_Compressed(1:(end-14)),'Other.mat'], 'layer_strain_bulk_xy');
LSB = layer_strain_bulk_xy;
MPW = 0.7; % min peak width (mm)
MPP = 5 ;% *(0.001); %8;%0.0012; % min peak prominence
nmax1 = 500; %500;
res_1 = zeros(1,nmax1);
res_2 = zeros(1,nmax1);
S_centre_1 = zeros(1,nmax1);
S_centre_2 = zeros(1,nmax1);
S_far_1 = zeros(1,nmax1);
S_far_2 = zeros(1,nmax1);
L_centre_1 = zeros(1,nmax1);
L_centre_2 = zeros(1,nmax1);
L_far_1 = zeros(1,nmax1);
L_far_2 = zeros(1,nmax1);
LSB_far_1 = zeros(1,nmax1);
LSB_far_2 = zeros(1,nmax1);
smoothing = 0;
if smoothing == 1
    %T4 = medfilt2(T3,[10 10]); % T3
    T4 = imgaussfilt(T3,[5 5]);
    %Lp4 = medfilt2(Lp3,[10 10]);
    %Lp4 = imgaussfilt(Lp4,[10 10]);
    %LSB = medfilt2(LSB,[10 10]);
    %LSB = imgaussfilt(LSB,[10 10]);
else
    T4 = T3; %T3
    %T4S = imgaussfilt(T3,[20 20]);
    Lp4 = Lp3;
    %Lp4S = imgaussfilt(Lp3,[20 20]);
    LSB = layer_strain_bulk_xy;
end

h = [1/2 1/2];
binomialCoeff = conv(h,h);
n_max = 100; %1750;% 2000
for n = 1:n_max
    binomialCoeff = conv(binomialCoeff,h);
end
fDelay = (length(binomialCoeff)-1)/2;

y_save = {};
%parfor n1 = 1:nmax1
nmax1 = 4;
for n1 = 1:nmax1
    output = rotateAround(T4, centre(1), centre(2), n1*180/nmax1);
    outputL = rotateAround(Lp4, centre(1), centre(2), n1*180/nmax1);
    %outputS = rotateAround(T4S, centre(1), centre(2), n1*180/nmax1);
    %outputLS = rotateAround(Lp4S, centre(1), centre(2), n1*180/nmax1);
    outputLSB = rotateAround(LSB, centre(1), centre(2), n1*180/nmax1);
 
    row = output(centre(1),:);
    rowL = outputL(centre(1),:);
    rowLSB = outputLSB(centre(1),:);
    
    S_centre_1(n1) = mean(row(centre(2)-inner_det_rad:centre(2)));
    S_centre_2(n1) = mean(row(centre(2):centre(2)+inner_det_rad));
    S_far_1(n1) = mean(row(centre(2)-outer_det_max_r:centre(2)-outer_det_min_r));
    S_far_2(n1) = mean(row(centre(2)+outer_det_min_r:centre(2)+outer_det_max_r));
    L_centre_1(n1) = mean(rowL(centre(2)-inner_det_rad:centre(2)));
    L_centre_2(n1) = mean(rowL(centre(2):centre(2)+inner_det_rad));
    L_far_1(n1) = mean(rowL(centre(2)-outer_det_max_r:centre(2)-outer_det_min_r));
    L_far_2(n1) = mean(rowL(centre(2)+outer_det_min_r:centre(2)+outer_det_max_r));
    LSB_far_1(n1) = mean(rowLSB(centre(2)-outer_det_max_r:centre(2)-outer_det_min_r));
    LSB_far_2(n1) = mean(rowLSB(centre(2)+outer_det_min_r:centre(2)+outer_det_max_r));
    
    flag = 0;
    n2 = 0;
    while flag == 0 && n2 < 20; % 20
        row = ceil(centre(1)) + n2*(-1)^n2;
        S0 = size(output(row,:),2);
        cdata1 = [output(row,:), zeros([1,n_max/2])];
        binomialMA = filter(binomialCoeff, 1, cdata1);
        cdata1 = diff(binomialMA);
        S1 = size(cdata1,2);
        x = (((fDelay+ceil(0.2*S0):S1-ceil(0.2*S0)))-fDelay)*config.SPACING_mm_pix(1);
        y = abs(cdata1(fDelay+ceil(0.2*S0):S1-ceil(0.2*S0)));
        figure
        plot(output(row,:))
        hold on
        plot(binomialMA(n_max/2:end))
        [~,m] = min(y(floor(size(x,2)*2/6):floor(size(x,2)*4/6)));
        m = m+floor(size(x,2)*2/6);
        [pks_l,locs_l,widths_l,proms_l] = findpeaks(y(1:m),'SortStr','descend','MinPeakWidth',MPW/config.SPACING_mm_pix(1),'MinPeakProminence',MPP);
        [pks_r,locs_r,widths_r,proms_r] = findpeaks(y(m:end),'SortStr','descend','MinPeakWidth',MPW/config.SPACING_mm_pix(1),'MinPeakProminence',MPP);
        if n1 == 1 && n2 == 0
            y_save{n1} = y;
        end
        if size(widths_l,2)>0 && size(widths_r,2)>0
            z = y(1:floor(locs_l(1)));
            [~,~,widths_l1,~] = findpeaks([min(z),z,min(z)],'SortStr','descend','MinPeakWidth',MPW/config.SPACING_mm_pix(1)/3,'MinPeakProminence',MPP);
            z = y(ceil(m+locs_r(1)):end);
            [~,~,widths_r1,~] = findpeaks([min(z),z,min(z)],'SortStr','descend','MinPeakWidth',MPW/config.SPACING_mm_pix(1)/3,'MinPeakProminence',MPP);
            z = y(ceil(locs_l(1)):m);
            [~,~,widths_l2,~] = findpeaks([min(z),z,min(z)],'SortStr','descend','MinPeakWidth',MPW/config.SPACING_mm_pix(1)/3,'MinPeakProminence',MPP);
            z = y(m:floor(m+locs_r(1)));
            [~,~,widths_r2,~] = findpeaks([min(z),z,min(z)],'SortStr','descend','MinPeakWidth',MPW/config.SPACING_mm_pix(1)/3,'MinPeakProminence',MPP);
%            figure
%            findpeaks(y,'Annotate','extents');
%             figure
%             findpeaks([0;y(1:floor(locs_l(1)));0],'Annotate','extents');
%             figure
%             findpeaks([0;y(ceil(m+locs_r(1)):end);0],'Annotate','extents');
%             figure
%             findpeaks([0;y(ceil(locs_l(1)):m);0],'Annotate','extents');
%             figure
%             findpeaks([0;y(m:floor(m+locs_r(1)));0],'Annotate','extents');
            if size(widths_l1,2)>0 && size(widths_r1,2)>0 && size(widths_l2,2)>0 && size(widths_r2,2)>0
                widths_l = [widths_l1(1), widths_l2(1)];
                widths_r = [widths_r1(1), widths_r2(1)];
                res_1(n1) = mean(widths_l)*2; %min
                res_2(n1) = mean(widths_r)*2; %min
                flag = 1;
            else
                n2 = n2 + 1;
            end
        else
            n2 = n2 + 1;
            
        end
    end
    figure
    findpeaks(y,'Annotate','extents');
end
figure
findpeaks(y_save{1},'Annotate','extents');
%%
res_s = [res_1, res_2]*config.SPACING_mm_pix(1);
res = mean(res_s);
StdNum = 7; %2;
logical = res_s > mean(res_s)+StdNum*std(res_s) |  res_s < mean(res_s)-StdNum*std(res_s);
res_s(logical) = res;
res_s(res_s == res) = mean(res_s);

a = 160;
b = 161;
%res_s = double(medfilt1(res_s,7,'truncate'));
%res_s = double(imgaussfilt(res_s,5));
%res_s = filter(binomialCoeff, 1, res_s);
%angle = [(180/nmax1):(180/nmax1):(180/nmax1*a),(180/nmax1*b):(180/nmax1):360];
res_s = [res_s(1:a),res_s(b:end)];
figure
plot(res_s)
res = mean(res_s); % resolution
res_std = std(res_s);
disp(res)
disp(res_std)

%%
L_centre = [L_centre_1, L_centre_2];
L_far = [L_far_1, L_far_2];
figure
plot(L_centre)
hold on
plot(L_far)
det = mean(L_far)-mean(L_centre); % detectability
det_std = std(L_far-L_centre);
%disp(det)
%disp(det_std)

S_centre = [S_centre_1, S_centre_2];
S_far = [S_far_1, S_far_2];
figure
plot(S_centre)
hold on
plot(S_far)
con = mean(S_centre)/mean(S_far); % stress contrast ratio
con_std = mean([std(S_centre./mean(S_far)),std(mean(S_centre)./S_far)]);
%disp(con)
%disp(con_std)

LSB_far = [LSB_far_1, LSB_far_2];
BS = mean(mean(LSB_far)); % bulk strain
BS_std = std2(LSB_far);
%disp(BS)
%disp(BS_std)
%%
disp([det,det_std,con,con_std,res,res_std])
disp([LayerThickness,BS,det,con,res])

%%
save([FilePath_Compressed(1:(end-14)),'TrueMetricsOP_Lin_NewRes.mat'],'BS','res','con','det','BS_std','res_std','con_std','det_std','centre','nmax1');

%%
pathpart = '231104_Inc_5_5mm3\231104_Inc_5_5mm3_10p_L';
load(['F:\Rhys\MATLAB Workspace\',pathpart,'\data\TrueMetricsOP_Lin_NewNewRes.mat'])
opcon = con;
opconstd = con_std;
opres = res;
opresstd = res_std;
load(['F:\Rhys\MATLAB Workspace\',pathpart,'\data\TrueMetricsCOP_Lin_NewNewRes.mat'])
disp([5.8,-BS,det,det_std,opcon,opconstd,opres,opresstd,con,con_std,res,res_std])

%% Resolution - full circle OP
% load([FilePath_Compressed(1:(end-14)),'LTImage.mat'], 'Lp3');
% load([FilePath_Compressed(1:(end-14)),'Other.mat'], 'layer_strain_bulk_xy');
LSB = layer_strain_bulk_xy;
nmax1 = 500; %500;
res_1 = zeros(1,nmax1);
res_2 = zeros(1,nmax1);
S_centre = zeros(1,nmax1);
S_far = zeros(1,nmax1);
L_centre = zeros(1,nmax1);
L_far = zeros(1,nmax1);
LSB_far = zeros(1,nmax1);

T4 = medfilt2(T3,[10 10]);
T4 = imgaussfilt(T4,[10 10]);
Lp4 = medfilt2(Lp3,[10 10]);
Lp4 = imgaussfilt(Lp4,[10 10]);
LSB = medfilt2(LSB,[10 10]);
LSB = imgaussfilt(LSB,[10 10]);

%T4 = T; % for COP

parfor n1 = 1:nmax1
    output = rotateAround(T4, centre(1), centre(2), n1*180/nmax1);
    outputL = rotateAround(Lp3, centre(1), centre(2), n1*180/nmax1);
    outputLSB = rotateAround(LSB, centre(1), centre(2), n1*180/nmax1);
    widths_row = zeros(1,Square_BR-Square_TL);
    
    S_centre(n1) = mean(output(centre(1)-inner_det_rad:centre(1)+inner_det_rad,centre(2)));
    S_far(n1) = mean(mean([output(centre(1)+outer_det_min_r:centre(1)+outer_det_max_r,centre(2)),output(centre(1)-outer_det_max_r:centre(1)-outer_det_min_r,centre(2))]));
    L_centre(n1) = mean(outputL(centre(1)-inner_det_rad:centre(1)+inner_det_rad,centre(2)));
    L_far(n1) = mean(mean([outputL(centre(1)+outer_det_min_r:centre(1)+outer_det_max_r,centre(2)),outputL(centre(1)-outer_det_max_r:centre(1)-outer_det_min_r,centre(2))]));
    LSB_far(n1) = mean(mean([outputLSB(centre(1)+outer_det_min_r:centre(1)+outer_det_max_r,centre(2)),outputLSB(centre(1)-outer_det_max_r:centre(1)-outer_det_min_r,centre(2))]));
    
    for row = Square_TL:Square_BR
        S0 = size(output(:,row));
        cdata1 = [output(:,row); zeros([1000,1])];
        binomialMA = filter(binomialCoeff, 1, cdata1);
        cdata1 = diff(binomialMA);
        S1 = size(cdata1);
        x = (((fDelay+ceil(0.1*S0):S1-ceil(0.1*S0)))-fDelay)*config.SPACING_mm_pix(1);
        y = abs(cdata1(fDelay+ceil(0.1*S0):S1-ceil(0.1*S0)));
        [~,locs,~,~] = findpeaks(y,x,'SortStr','descend','MinPeakWidth',MPW,'MinPeakProminence',MPP);
        
        if size(locs,2)>1
                width = abs(locs(1)-locs(2));
                widths_row(row-Square_TL+1) = width;
                %xp = x;
                %yp = y;
        elseif row ~= Square_TL
                widths_row(row-Square_TL+1) = 0;
        end
    end
    %figure
    %findpeaks(yp,xp,'Annotate','extents')
    gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
    startPoints = [mean(widths_row) (Square_BR-Square_TL)/2 (Square_BR-Square_TL)/2 0];
    LowBounds = [0, Square_TL, 0, 0];
    UpBounds = [max(widths_row)*2, Square_BR, Square_BR-Square_TL, 0];
    f_row = fit((Square_TL:Square_BR)',(widths_row)',gaussEqn,'Start', startPoints,'Lower', LowBounds,'Upper', UpBounds);
    row_cs = coeffvalues(f_row);
    row = ceil(row_cs(2));
    %figure
    %plot(f_row,(Square_TL:Square_BR)',(widths_row)')
    
    %row_cs(2) = centre(1);
    
    flag = 0;
    n2 = 0;
    while flag == 0 && n2 < 50;
        row = ceil(row_cs(2)) + n2*(-1)^n2;
        S0 = size(output(:,row));
        cdata1 = [output(:,row); zeros([1000,1])];
        binomialMA = filter(binomialCoeff, 1, cdata1);
        cdata1 = diff(binomialMA);
        S1 = size(cdata1);
        x = (((fDelay+ceil(0.1*S0):S1-ceil(0.1*S0)))-fDelay)*config.SPACING_mm_pix(1);
        y = abs(cdata1(fDelay+ceil(0.1*S0):S1-ceil(0.1*S0)));
        [pks,locs,widths,proms] = findpeaks(y,x,'SortStr','descend','MinPeakWidth',MPW,'MinPeakProminence',MPP);
        if size(locs,2)>1
            if locs(1) < locs(2)
                res_1(n1) = widths(1);
                res_2(n1) = widths(2);
            else
                res_1(n1) = widths(2);
                res_2(n1) = widths(1);
            end
            flag = 1;
        else
            n2 = n2 + 1;
        end
    end
end

%figure
%findpeaks(y,x,'Annotate','extents')

%%
res_s = [res_1, res_2];
figure
plot(res_s)
res = mean(res_s); % resolution
res_std = std(res_s);
%disp(res)
%disp(res_std)

figure
plot(L_centre)
hold on
plot(L_far)
det = mean(L_far)-mean(L_centre); % detectability
det_std = std(L_far-L_centre);
%disp(det)
%disp(det_std)

figure
plot(S_centre)
hold on
plot(S_far)
con = mean(S_centre)/mean(S_far); % stress contrast ratio
con_std = mean([std(S_centre./mean(S_far)),std(mean(S_centre)./S_far)]);
%disp(con)
%disp(con_std)

BS = mean(mean(LSB_far)); % bulk strain
BS_std = std2(LSB_far);
%disp(BS)
%disp(BS_std)
%%
disp([det,det_std,con,con_std,res,res_std])
disp([LayerThickness,BS,det,con,res])
%%
save([FilePath_Compressed(1:(end-14)),'TrueMetricsOPSig_smoothed10x10.mat'],'BS','res','con','det','BS_std','res_std','con_std','det_std','centre','nmax1');

%% Stress Map
%load([FilePath_Compressed(1:(end-14)),'230907_Inc_4_2_5mm2_10p_lowrestest\COP_T.mat'], 'T');
cdata1 = -T*1000';
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
image(cdata1,'Parent',axes1,'CDataMapping','scaled');
S = size(cdata1);
xlim(axes1,[0.5 S(1)+0.5]);
ylim(axes1,[0.5 S(2)+0.5]);
title('L: 5mm#2, S: Inc 4, P: 10%, 15.5mmx15.5mm, kPa');
box(axes1,'on');
axis(axes1,'ij');
set(axes1,'Layer','top');
colorbar('peer',axes1);

%% Resolution find centre COP
load([FilePath_Compressed(1:(end-14)),'231025_Inc_5_5mm2_10p_lowrestest\COP_T.mat'], 'T');
S3 = -T*1000;
Square_TL = 450;
Square_BR = 800;
MPW = 0.6; % min peak width (mm)
MPP = 0.005; %0.01; % min peak prominence
centre = [0,0];
widths_row = zeros(1,Square_BR-Square_TL);
widths_col = zeros(1,Square_BR-Square_TL);
h = [1/2 1/2];
binomialCoeff = conv(h,h);
n_max = 2000;
for n = 1:n_max
    binomialCoeff = conv(binomialCoeff,h);
end
fDelay = (length(binomialCoeff)-1)/2;
for row = Square_TL:Square_BR
    S0 = size(S3(:,row));
    cdata1 = [S3(:,row); zeros([1000,1])];
    binomialMA = filter(binomialCoeff, 1, cdata1);
    cdata1 = diff(binomialMA);
    S1 = size(cdata1);
    x = (((fDelay+ceil(0.1*S0):S1-ceil(0.1*S0)))-fDelay)*config.SPACING_mm_pix(1);
    y = abs(cdata1(fDelay+ceil(0.1*S0):S1-ceil(0.1*S0)));
    [~,locs,~,~] = findpeaks(y,x,'SortStr','descend','MinPeakWidth',MPW,'MinPeakProminence',MPP);
    if size(locs,2)>1
            width = abs(locs(1)-locs(2));
            widths_row(row-Square_TL+1) = width;
            centre(1) = row;
    elseif row ~= Square_TL
            widths_row(row-Square_TL+1) = 0; % widths_row(row-Square_TL);
    end
    if row == 500
        figure
        findpeaks(y,x,'Annotate','extents')
    end
end
S4 = S3';
for col = Square_TL:Square_BR
    S0 = size(S4(:,col));
    cdata1 = [S4(:,col); zeros([1000,1])];
    binomialMA = filter(binomialCoeff, 1, cdata1);
    cdata1 = diff(binomialMA);
    S1 = size(cdata1);
    x = (((fDelay+ceil(0.1*S0):S1-ceil(0.1*S0)))-fDelay)*config.SPACING_mm_pix(1);
    y = abs(cdata1(fDelay+ceil(0.1*S0):S1-ceil(0.1*S0)));
    [~,locs,~,~] = findpeaks(y,x,'SortStr','descend','MinPeakWidth',MPW,'MinPeakProminence',MPP);
    if size(locs,2)>1
            width = abs(locs(1)-locs(2));
            widths_col(col-Square_TL+1) = width;
            centre(2) = col;
    elseif col ~= Square_TL
            widths_col(col-Square_TL+1) = 0; % widths_col(row-Square_TL);
    end
%     if col == 471
%         figure
%         findpeaks(y,x,'Annotate','extents')
%     end
end
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
startPoints = [mean(widths_row) (Square_BR-Square_TL)/2 (Square_BR-Square_TL)/2 0];
LowBounds = [0, Square_TL, 0, 0];
UpBounds = [max(widths_row)*2, Square_BR, Square_BR-Square_TL, 0];
f_row = fit((Square_TL:Square_BR)',(widths_row)',gaussEqn,'Start', startPoints,'Lower', LowBounds,'Upper', UpBounds); %, 'Exclude', x < 800);
f_col = fit((Square_TL:Square_BR)',(widths_col)',gaussEqn,'Start', startPoints,'Lower', LowBounds,'Upper', UpBounds); %, 'Exclude', x < 800);
row_cs = coeffvalues(f_row);
col_cs = coeffvalues(f_col);
centre_C = [ceil(row_cs(2)) ceil(col_cs(2))];
figure
plot(f_row,(Square_TL:Square_BR)',(widths_row)')
figure
plot(f_col,(Square_TL:Square_BR)',(widths_col)')
disp(centre_C)

%% Stress Map - Choose ROI's - COP (should choose same, if I figure out how to not take the extended COP image, and we use same centre?)
% load([FilePath_Compressed(1:(end-14)),'COP_T.mat'], 'T');
% S3 = -T*1000;
%centre = [608,633];
inner_det_rad = 20; % pixels radius for centre point
outer_det_min_r = 500; % radial range for far point
outer_det_max_r = 550; % radial range for far point
output = rotateAround(S3, centre_C(1), centre_C(2), 180);
cdata1 = output';
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
image(cdata1,'Parent',axes1,'CDataMapping','scaled');
S = size(cdata1);
xlim(axes1,[0.5 S(1)+0.5]);
ylim(axes1,[0.5 S(2)+0.5]);
title('L: 5mm#2, S: Inc 4, P: 10%, 15.5mmx15.5mm, kPa');
box(axes1,'on');
axis(axes1,'ij');
set(axes1,'Layer','top');
colorbar('peer',axes1);
x = centre_C(1) - inner_det_rad; y = centre_C(2) - inner_det_rad;
d = inner_det_rad * 2;
pos = [x y d d];
rectangle('Position',pos,'Curvature',[1 1])
x = centre_C(1) - outer_det_min_r; y = centre_C(2) - outer_det_min_r;
d = outer_det_min_r * 2;
pos = [x y d d];
rectangle('Position',pos,'Curvature',[1 1])
x = centre_C(1) - outer_det_max_r; y = centre_C(2) - outer_det_max_r;
d = outer_det_max_r * 2;
pos = [x y d d];
rectangle('Position',pos,'Curvature',[1 1])

%% Resolution - full circle COP
% load([FilePath_Compressed(1:(end-14)),'COP_T.mat'], 'T');
% S3 = -T*1000;
nmax1 = 20;
res_1 = zeros(1,nmax1);
res_2 = zeros(1,nmax1);
S_centre = zeros(1,nmax1);
S_far = zeros(1,nmax1);

for n1 = 1:nmax1
    output = rotateAround(S3, centre_C(1), centre_C(2), n1*180/nmax1);
    widths_row = zeros(1,Square_BR-Square_TL);
    
    S_centre(n1) = mean(output(centre_C(1)-inner_det_rad:centre_C(1)+inner_det_rad,centre_C(2)));
    S_far(n1) = mean(mean([output(centre_C(1)+outer_det_min_r:centre_C(1)+outer_det_max_r,centre_C(2)),output(centre_C(1)-outer_det_max_r:centre_C(1)-outer_det_min_r,centre_C(2))]));
    
    for row = Square_TL:Square_BR
        S0 = size(output(:,row));
        cdata1 = [output(:,row); zeros([1000,1])];
        binomialMA = filter(binomialCoeff, 1, cdata1);
        cdata1 = diff(binomialMA);
        S1 = size(cdata1);
        x = (((fDelay+ceil(0.1*S0):S1-ceil(0.1*S0)))-fDelay)*config.SPACING_mm_pix(1);
        y = abs(cdata1(fDelay+ceil(0.1*S0):S1-ceil(0.1*S0)));
        [~,locs,~,~] = findpeaks(y,x,'SortStr','descend','MinPeakWidth',MPW,'MinPeakProminence',MPP);
        if size(locs,2)>1
                width = abs(locs(1)-locs(2));
                widths_row(row-Square_TL+1) = width;
        elseif row ~= Square_TL
                widths_row(row-Square_TL+1) = 0;
        end
    end
    gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
    startPoints = [mean(widths_row) (Square_BR-Square_TL)/2 (Square_BR-Square_TL)/2 0];
    LowBounds = [0, Square_TL, 0, 0];
    UpBounds = [max(widths_row)*2, Square_BR, Square_BR-Square_TL, 0];
    f_row = fit((Square_TL:Square_BR)',(widths_row)',gaussEqn,'Start', startPoints,'Lower', LowBounds,'Upper', UpBounds);
    row_cs = coeffvalues(f_row);
    row = ceil(row_cs(2));
    
    flag = 0;
    n2 = 0;
    while flag == 0 && n2 < 50;
        row = ceil(row_cs(2)) + n2*(-1)^n2;
        S0 = size(output(:,row));
        cdata1 = [output(:,row); zeros([1000,1])];
        binomialMA = filter(binomialCoeff, 1, cdata1);
        cdata1 = diff(binomialMA);
        S1 = size(cdata1);
        x = (((fDelay+ceil(0.1*S0):S1-ceil(0.1*S0)))-fDelay)*config.SPACING_mm_pix(1);
        y = abs(cdata1(fDelay+ceil(0.1*S0):S1-ceil(0.1*S0)));
        [pks,locs,widths,proms] = findpeaks(y,x,'SortStr','descend','MinPeakWidth',MPW,'MinPeakProminence',MPP);
        if size(locs,2)>1
            res_1(n1) = widths(1);
            res_2(n1) = widths(2);
            flag = 1;
        else
            n2 = n2 + 1;
        end
    end
end

res_C = [res_1, res_2];
figure
plot(res_C)
res_C = mean(res_C);
disp(res_C)

figure
plot(S_centre)
hold on
plot(S_far)
con_C = mean(S_centre)/mean(S_far);
disp(con_C)

save([FilePath_Compressed(1:(end-14)),meta.job,'\','MetricsCOP.mat'],'res_C','con_C','centre_C');
