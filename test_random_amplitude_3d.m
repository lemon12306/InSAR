clear; close all; clc;

% 多组随机幅度无噪声实验：几何保持一致，仅改变散射点幅度。
geometrySeed = 20261012;
amplitudeSeeds = [1101 2202 3303];
Nrun = numel(amplitudeSeeds);

rng(geometrySeed,'twister');
SAR_Parameters;

ThetaTa = deg2rad(55);
Tp = 30e-9;
Nfft_az = 16384;
dx_insar = 0.05;
filterSize = 5;

Rr = Rr(:);
x_rda = V*((-Nfft_az/2:Nfft_az/2-1)/PRF);
ROW_idx = find(Rr>=40 & Rr<=60);
COL_idx = find(x_rda>=-10 & x_rda<=10);
dx_rda = abs(x_rda(2)-x_rda(1));
col_step = max(1,round(dx_insar/dx_rda));
COL_idx = COL_idx(1:col_step:end);

Ax_r = x_rda(COL_idx);
Ay_r = Rr(ROW_idx);

zA = Antenna10(3);
zB = Antenna20(3);
yA = Antenna10(2);
yB = Antenna20(2);

Y_ground = sqrt(max(Ay_r.^2-(Z0-zA)^2,0));
[Xgrid,Ygrid] = meshgrid(Ax_r,Y_ground);

RrxA_Z0 = sqrt((Ygrid-yA).^2+(Z0-zA).^2);
RrxB_Z0 = sqrt((Ygrid-yB).^2+(Z0-zB).^2);
phi_flat = (2*pi/lambda).*(RrxB_Z0-RrxA_Z0);
sens = (2*pi/lambda).*((Z0-zB)./RrxB_Z0- ...
    (Z0-zA)./RrxA_Z0);

h_true = BaseHeight+HeightHill.*exp( ...
    -0.5*((Xgrid-xc)/s_x).^2 ...
    -0.5*((Ygrid-yc)/sigma_y).^2);

ref_mask_x = Ax_r>=-10 & Ax_r<=-8;
ref_mask_y = Ay_r>=40 & Ay_r<=43;
ref_area = ref_mask_y(:) & ref_mask_x(:).';
win = ones(filterSize);
[~,centerCol] = min(abs(Ax_r-xc));

pointGeometry = point;
hFull = cell(Nrun,1);
hBit = cell(Nrun,1);
profileError = zeros(numel(Y_ground),Nrun);

AmpMean = zeros(Nrun,1);
AmpStd = zeros(Nrun,1);
HeightRMSE_Full = zeros(Nrun,1);
HeightRMSE_Bit = zeros(Nrun,1);
RippleStd = zeros(Nrun,1);
DominantWavelength = zeros(Nrun,1);

for irun = 1:Nrun
    fprintf('\n========== 随机幅度实验 %d/%d，seed=%d ==========\n', ...
        irun,Nrun,amplitudeSeeds(irun));

    rng(amplitudeSeeds(irun),'twister');
    randomAmp = 0.5 + rand(size(pointGeometry,1),1);
    pointRun = pointGeometry;
    pointRun(:,6) = randomAmp;

    AmpMean(irun) = mean(randomAmp);
    AmpStd(irun) = std(randomAmp);

    [Echo,~] = RawDataGen_rand( ...
        pointRun,Rr,nan,B,fc,PRF,Ym,H,V,ThetaTa,Bz,0,Tp,'raw');

    Echo_A = Echo(:,:,1);
    Echo_B = Echo(:,:,2);
    clear Echo pointRun

    Echo_A_bit = QuantizeZero(Echo_A);
    Echo_B_bit = QuantizeZero(Echo_B);

    A_full = RunRDAAndCrop( ...
        Echo_A,Rr,B,Tp,C,PRF,V,lambda,Nfft_az,ROW_idx,COL_idx);
    B_full = RunRDAAndCrop( ...
        Echo_B,Rr,B,Tp,C,PRF,V,lambda,Nfft_az,ROW_idx,COL_idx);
    A_bit = RunRDAAndCrop( ...
        Echo_A_bit,Rr,B,Tp,C,PRF,V,lambda,Nfft_az,ROW_idx,COL_idx);
    B_bit = RunRDAAndCrop( ...
        Echo_B_bit,Rr,B,Tp,C,PRF,V,lambda,Nfft_az,ROW_idx,COL_idx);

    clear Echo_A Echo_B Echo_A_bit Echo_B_bit

    [~,phiFilt_full,coh_full] = ProcessInSAR( ...
        A_full,B_full,phi_flat,win);
    [~,phiFilt_bit,coh_bit] = ProcessInSAR( ...
        A_bit,B_bit,phi_flat,win);

    amp_ref = sqrt(abs(A_full).*abs(B_full));
    validMask = amp_ref>=0.05*max(amp_ref(:)) & ...
        coh_full>=0.5 & coh_bit>=0.5 & ...
        abs(sens)>1e-8 & isfinite(h_true);

    phi_ref_full = median(phiFilt_full(ref_area),'omitnan');
    phi_ref_bit = median(phiFilt_bit(ref_area),'omitnan');

    hFull{irun} = Z0+(phiFilt_full-phi_ref_full)./sens;
    hBit{irun} = Z0+(phiFilt_bit-phi_ref_bit)./sens;

    errFull = hFull{irun}-h_true;
    errBit = hBit{irun}-h_true;
    HeightRMSE_Full(irun) = sqrt(mean(errFull(validMask).^2));
    HeightRMSE_Bit(irun) = sqrt(mean(errBit(validMask).^2));

    profileError(:,irun) = errBit(:,centerCol);
    centralMask = Y_ground>=42 & Y_ground<=58;
    RippleStd(irun) = std(profileError(centralMask,irun),'omitnan');
    DominantWavelength(irun) = GetDominantWavelength( ...
        Y_ground(centralMask),profileError(centralMask,irun));

    fprintf('幅度均值/标准差 = %.4f / %.4f\n', ...
        AmpMean(irun),AmpStd(irun));
    fprintf('无噪声全精度/单比特高程RMSE = %.4f / %.4f m\n', ...
        HeightRMSE_Full(irun),HeightRMSE_Bit(irun));
    fprintf('中心剖面误差标准差 = %.4f m，主导波长 = %.4f m\n', ...
        RippleStd(irun),DominantWavelength(irun));

    clear A_full B_full A_bit B_bit phiFilt_full phiFilt_bit
    clear coh_full coh_bit amp_ref validMask errFull errBit
end

ProfileErrorCorr = corrcoef(profileError,'Rows','pairwise');
disp('中心剖面误差的组间相关系数：');
disp(ProfileErrorCorr);

Results = table(amplitudeSeeds(:),AmpMean,AmpStd, ...
    HeightRMSE_Full,HeightRMSE_Bit,RippleStd,DominantWavelength, ...
    'VariableNames',{'AmplitudeSeed','AmplitudeMean','AmplitudeStd', ...
    'HeightRMSE_Full','HeightRMSE_Bit','ProfileRippleStd', ...
    'DominantWavelength'});

writetable(Results,'随机幅度_无噪声_统计.csv');
writematrix(ProfileErrorCorr,'随机幅度_剖面误差相关系数.csv');
save('随机幅度_无噪声_高程结果.mat', ...
    'amplitudeSeeds','Ax_r','Y_ground','h_true','hFull','hBit', ...
    'profileError','ProfileErrorCorr','Results');

figure('Color','w','Position',[120 60 1200 920]);
for irun = 1:Nrun
    subplot(Nrun,1,irun);
    plot(Y_ground,h_true(:,centerCol),'k-','LineWidth',1.8);
    hold on;
    plot(Y_ground,hFull{irun}(:,centerCol), ...
        'Color',[0.8500 0.3250 0.0980],'LineWidth',1.3);
    plot(Y_ground,hBit{irun}(:,centerCol),'--', ...
        'Color',[0.9290 0.6940 0.1250],'LineWidth',1.2);
    title(sprintf('随机幅度 seed=%d，主导波长=%.2f m', ...
        amplitudeSeeds(irun),DominantWavelength(irun)));
    xlabel('距离向/m');
    ylabel('高程/m');
    legend('真实高程','全精度','零阈值单比特','Location','best');
    grid on;
end
sgtitle('不同随机幅度实现下的无噪声山丘中心高程剖面');
exportgraphics(gcf,'随机幅度_无噪声_中心高程剖面对比.png','Resolution',300);

figure('Color','w','Position',[60 120 1500 520]);
for irun = 1:Nrun
    subplot(1,Nrun,irun);
    PlotHeightSurface(Ax_r,Y_ground,hBit{irun});
    title(sprintf('单比特三维高程，seed=%d',amplitudeSeeds(irun)));
end
sgtitle('不同随机幅度实现下的无噪声零阈值单比特三维高程');
exportgraphics(gcf,'随机幅度_无噪声_单比特三维高程.png','Resolution',300);

figure('Color','w','Position',[60 120 1500 520]);
subplot(1,3,1);
PlotHeightSurface(Ax_r,Y_ground,h_true);
title('真实三维高程');
subplot(1,3,2);
PlotHeightSurface(Ax_r,Y_ground,hFull{1});
title('无噪声全精度三维高程');
subplot(1,3,3);
PlotHeightSurface(Ax_r,Y_ground,hBit{1});
title(sprintf('无噪声单比特三维高程，seed=%d',amplitudeSeeds(1)));
sgtitle('真实、全精度与零阈值单比特三维高程对比');
exportgraphics(gcf,'随机幅度_真实全精度单比特_三维高程对比.png', ...
    'Resolution',300);

fprintf('\n随机幅度多次无噪声实验及三维绘图完成。\n');

function Data = QuantizeZero(Data)
Data = complex(2*double(real(Data)>=0)-1, ...
    2*double(imag(Data)>=0)-1);
end

function cropped = RunRDAAndCrop( ...
    Echo,Rr,B,Tp,C,PRF,V,lambda,Nfft_az,ROW_idx,COL_idx)
[img,~] = RunRDA(Echo,Rr,B,Tp,C,PRF,V,lambda,Nfft_az);
cropped = img(ROW_idx,COL_idx);
end

function HRRP = RangeCompress(Echo,Rr,B,Tp,C)
nrn = size(Echo,1);
dR = mean(diff(Rr));
dt = 2*dR/C;
Kr = B/Tp;
halfN = ceil(Tp/(2*dt));
t_ref = (-halfN:halfN).'*dt;
s_tx = complex(zeros(length(t_ref),1));
mask = abs(t_ref)<=Tp/2;
s_tx(mask) = exp(1j*pi*Kr*t_ref(mask).^2);
h_mf = conj(flipud(s_tx));
M = length(h_mf);
Nfft = 2^nextpow2(nrn+M-1);
Y_full = ifft(fft(Echo,Nfft,1).*fft(h_mf,Nfft,1),Nfft,1);
startIdx = floor(M/2)+1;
HRRP = Y_full(startIdx:startIdx+nrn-1,:);
end

function [img,x_rda] = RunRDA( ...
    Echo,Rr,B,Tp,C,PRF,V,lambda,Nfft_az)
HRRP = RangeCompress(Echo,Rr,B,Tp,C);
[nrn,naz] = size(HRRP);
az_win = hann(naz,'periodic').';
HRRP = HRRP.*az_win;
A_pad = complex(zeros(nrn,Nfft_az));
i0 = floor((Nfft_az-naz)/2);
A_pad(:,i0+1:i0+naz) = HRRP;
fa = (-Nfft_az/2:Nfft_az/2-1)*PRF/Nfft_az;
S = fftshift(fft(ifftshift(A_pad,2),Nfft_az,2),2);
clear A_pad HRRP
S_rcmc = complex(zeros(size(S)));
for k = 1:Nfft_az
    D2 = 1-(lambda*fa(k)/(2*V))^2;
    if D2<=0
        continue;
    end
    D = sqrt(D2);
    R_query = Rr/D;
    S_rcmc(:,k) = interp1(Rr,S(:,k),R_query,'linear',0);
end
clear S
for r = 1:nrn
    Ka = 2*V^2/(lambda*Rr(r));
    H_az = exp(-1j*pi*fa.^2/Ka);
    S_rcmc(r,:) = S_rcmc(r,:).*H_az;
end
img = fftshift(ifft(ifftshift(S_rcmc,2),Nfft_az,2),2);
x_rda = V*((-Nfft_az/2:Nfft_az/2-1)/PRF);
end

function [phi_deflat,phi_filt,coh] = ProcessInSAR(A,B,phi_flat,win)
I_reg = A.*conj(B);
I_deflat = I_reg.*exp(-1j*phi_flat);
phi_deflat = angle(I_deflat);
I_unit = I_deflat./(abs(I_deflat)+eps);
I_filt = conv2(I_unit,win,'same');
phi_filt = angle(I_filt);
num = abs(conv2(A.*conj(B),win,'same'));
denA = conv2(abs(A).^2,win,'same');
denB = conv2(abs(B).^2,win,'same');
coh = num./sqrt(denA.*denB+eps);
end

function wavelength = GetDominantWavelength(y,e)
valid = isfinite(y) & isfinite(e);
y = y(valid);
e = e(valid);
e = e-mean(e);
N = numel(e);
dy = mean(diff(y));
w = hann(N,'periodic');
powerSpectrum = abs(fft(e.*w)).^2;
f = (0:N-1)'/(N*dy);
band = f>=0.2 & f<=2.5;
idx = find(band);
[~,imax] = max(powerSpectrum(idx));
wavelength = 1/f(idx(imax));
end

function PlotHeightSurface(Ax,Yg,Hmap)
surf(Ax,Yg,Hmap,'EdgeColor','none');
shading interp;
axis tight;
view(45,35);
colormap(parula);
colorbar;
xlabel('方位向/m');
ylabel('距离向/m');
zlabel('高程/m');
grid on;
end
