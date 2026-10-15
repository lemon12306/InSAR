
clear; close all; clc;
rng(20261012);
SAR_Parameters;

%%%%%%%%%%%%%%%%%%%%%%%%%% 1. 实验参数
ThetaTa = deg2rad(55);
Tp = 30e-9;

SNR_list = [10 15 20 25 30 Inf];
Nsnr = length(SNR_list);

SNR_plot = 1:Nsnr;
SNR_labels = {'10','15','20','25','30','无噪声'};

dx_insar = 0.05;
filterSize = 5;
Nfft_az = 16384;

showSNR = [10 20 30 Inf];
Nshow = length(showSNR);

%%%%%%%%%%%%%%%%%%%%%%%%%% 2. 生成原始回波
fprintf('\n========== 生成原始回波 ==========\n');

[Echo,~] = RawDataGen_rand( ...
    point,Rr,nan,B,fc,PRF,Ym,H,V,ThetaTa,Bz,0,Tp,'raw');

Echo_A = Echo(:,:,1);
Echo_B = Echo(:,:,2);

Rr = Rr(:);

%%%%%%%%%%%%%%%%%%%%%%%%%% 3. 接收噪声设置
maskA = Echo_A~=0;
maskB = Echo_B~=0;

PowerA = mean(abs(Echo_A(maskA)).^2);
PowerB = mean(abs(Echo_B(maskB)).^2);

P_ref = (PowerA+PowerB)/2;

Noise_A0 = (randn(size(Echo_A))+ ...
    1j*randn(size(Echo_A)))/sqrt(2);

Noise_B0 = (randn(size(Echo_B))+ ...
    1j*randn(size(Echo_B)))/sqrt(2);

Noise_A0 = Noise_A0/ ...
    sqrt(mean(abs(Noise_A0(:)).^2));

Noise_B0 = Noise_B0/ ...
    sqrt(mean(abs(Noise_B0(:)).^2));

fprintf('参考信号功率 = %.6e\n',P_ref);
fprintf('SNR设置：10、15、20、25、30dB及无噪声\n');

%%%%%%%%%%%%%%%%%%%%%%%%%% 4. 固定评价坐标
x_rda = V*((-Nfft_az/2:Nfft_az/2-1)/PRF);

ROW_idx = find(Rr>=40 & Rr<=60);
COL_idx = find(x_rda>=-10 & x_rda<=10);

dx_rda = abs(x_rda(2)-x_rda(1));
col_step = max(1,round(dx_insar/dx_rda));

COL_idx = COL_idx(1:col_step:end);

Ax_r = x_rda(COL_idx);
Ay_r = Rr(ROW_idx);

%%%%%%%%%%%%%%%%%%%%%%%%%% 5. 地形几何及真实高程
zA = Antenna10(3);
zB = Antenna20(3);

yA = Antenna10(2);
yB = Antenna20(2);

Y_ground = sqrt(max( ...
    Ay_r.^2-(Z0-zA)^2,0));

[Xgrid,Ygrid] = meshgrid(Ax_r,Y_ground);

RrxA_Z0 = sqrt((Ygrid-yA).^2+(Z0-zA).^2);
RrxB_Z0 = sqrt((Ygrid-yB).^2+(Z0-zB).^2);

phi_flat = (2*pi/lambda).* ...
    (RrxB_Z0-RrxA_Z0);

sens = (2*pi/lambda).* ...
    ((Z0-zB)./RrxB_Z0- ...
     (Z0-zA)./RrxA_Z0);

h_true = BaseHeight+HeightHill.*exp( ...
    -0.5*((Xgrid-xc)/s_x).^2 ...
    -0.5*((Ygrid-yc)/sigma_y).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%% 6. 参考区域
ref_mask_x = Ax_r>=-10 & Ax_r<=-8;
ref_mask_y = Ay_r>=40 & Ay_r<=43;

ref_area = ref_mask_y(:) & ref_mask_x(:).';

%%%%%%%%%%%%%%%%%%%%%%%%%% 7. 固定评价掩码
fprintf('\n========== 构建固定评价掩码 ==========\n');

[imgA_ref,~] = RunRDA( ...
    Echo_A,Rr,B,Tp,C,PRF,V,lambda,Nfft_az);

[imgB_ref,~] = RunRDA( ...
    Echo_B,Rr,B,Tp,C,PRF,V,lambda,Nfft_az);

A_ref = imgA_ref(ROW_idx,COL_idx);
B_ref = imgB_ref(ROW_idx,COL_idx);

clear imgA_ref imgB_ref

win = ones(filterSize);

[~,~,coh_ref] = ProcessInSAR( ...
    A_ref,B_ref,phi_flat,win);

amp_ref = sqrt(abs(A_ref).*abs(B_ref));

validMask = amp_ref>=0.05*max(amp_ref(:)) & ...
            coh_ref>=0.5 & ...
            abs(sens)>1e-8 & ...
            isfinite(h_true);

fprintf('固定有效像素数 = %d\n',nnz(validMask));
fprintf('有效像素比例 = %.4f%%\n', ...
    100*mean(validMask(:)));

clear A_ref B_ref

%%%%%%%%%%%%%%%%%%%%%%%%%% 8. 统计变量
AmpCorr = zeros(Nsnr,2);
ComplexCorr = zeros(Nsnr,2);

BackgroundFull = zeros(Nsnr,1);
BackgroundBit = zeros(Nsnr,1);

CohFull = zeros(Nsnr,1);
CohBit = zeros(Nsnr,1);

PhaseRMSE_before = zeros(Nsnr,1);
PhaseRMSE_after = zeros(Nsnr,1);

HeightRMSE_full = zeros(Nsnr,1);
HeightRMSE_bit = zeros(Nsnr,1);

HeightMAE_full = zeros(Nsnr,1);
HeightMAE_bit = zeros(Nsnr,1);

HeightP95_full = zeros(Nsnr,1);
HeightP95_bit = zeros(Nsnr,1);

Outlier_full = zeros(Nsnr,1);
Outlier_bit = zeros(Nsnr,1);

HillRMSE_full = zeros(Nsnr,1);
HillRMSE_bit = zeros(Nsnr,1);

FlatRMSE_full = zeros(Nsnr,1);
FlatRMSE_bit = zeros(Nsnr,1);

h_full_show = cell(Nshow,1);
h_bit_show = cell(Nshow,1);

[Xall,Rall] = meshgrid(x_rda,Rr);
targetMask = abs(Xall)<=10 & ...
             Rall>=40 & Rall<=60;

%%%%%%%%%%%%%%%%%%%%%%%%%% 9. SNR循环
for k = 1:Nsnr

    snr_dB = SNR_list(k);

    if isinf(snr_dB)

        fprintf('\n========== 无噪声 ==========\n');

        Input_A = Echo_A;
        Input_B = Echo_B;

    else

        fprintf('\n========== SNR = %.1f dB ==========\n',snr_dB);

        NoisePower = P_ref/10^(snr_dB/10);

        Noise_A = sqrt(NoisePower)*Noise_A0;
        Noise_B = sqrt(NoisePower)*Noise_B0;

        Input_A = Echo_A+Noise_A;
        Input_B = Echo_B+Noise_B;

        clear Noise_A Noise_B

    end

    %%%%%%%%%%%%%%%%%%%%% 9.1 零阈值单比特量化
    Input_A_bit = QuantizeZero(Input_A);
    Input_B_bit = QuantizeZero(Input_B);

    %%%%%%%%%%%%%%%%%%%%% 9.2 含噪全精度RDA
    [imgA_full,~] = RunRDA( ...
        Input_A,Rr,B,Tp,C,PRF,V,lambda,Nfft_az);

    [imgB_full,~] = RunRDA( ...
        Input_B,Rr,B,Tp,C,PRF,V,lambda,Nfft_az);

    %%%%%%%%%%%%%%%%%%%%% 9.3 零阈值单比特RDA
    [imgA_bit,~] = RunRDA( ...
        Input_A_bit,Rr,B,Tp,C,PRF,V,lambda,Nfft_az);

    [imgB_bit,~] = RunRDA( ...
        Input_B_bit,Rr,B,Tp,C,PRF,V,lambda,Nfft_az);

    clear Input_A Input_B Input_A_bit Input_B_bit

    %%%%%%%%%%%%%%%%%%%%% 9.4 RDA相关系数
    AmpCorr(k,1) = GetAmplitudeCorr( ...
        imgA_full,imgA_bit);

    AmpCorr(k,2) = GetAmplitudeCorr( ...
        imgB_full,imgB_bit);

    ComplexCorr(k,1) = GetComplexCorr( ...
        imgA_full,imgA_bit);

    ComplexCorr(k,2) = GetComplexCorr( ...
        imgB_full,imgB_bit);

    %%%%%%%%%%%%%%%%%%%%% 9.5 背景能量
    pFullA = abs(imgA_full).^2;
    pFullB = abs(imgB_full).^2;

    pBitA = abs(imgA_bit).^2;
    pBitB = abs(imgB_bit).^2;

    BackgroundFull(k) = 50*( ...
        sum(pFullA(~targetMask))/sum(pFullA(:)) + ...
        sum(pFullB(~targetMask))/sum(pFullB(:)));

    BackgroundBit(k) = 50*( ...
        sum(pBitA(~targetMask))/sum(pBitA(:)) + ...
        sum(pBitB(~targetMask))/sum(pBitB(:)));

    clear pFullA pFullB pBitA pBitB

    %%%%%%%%%%%%%%%%%%%%% 9.6 裁剪RDA
    A_full = imgA_full(ROW_idx,COL_idx);
    B_full = imgB_full(ROW_idx,COL_idx);

    A_bit = imgA_bit(ROW_idx,COL_idx);
    B_bit = imgB_bit(ROW_idx,COL_idx);

    clear imgA_full imgB_full imgA_bit imgB_bit

    %%%%%%%%%%%%%%%%%%%%% 9.7 InSAR处理
    [phiDef_full,phiFilt_full,coh_full] = ...
        ProcessInSAR(A_full,B_full,phi_flat,win);

    [phiDef_bit,phiFilt_bit,coh_bit] = ...
        ProcessInSAR(A_bit,B_bit,phi_flat,win);

    CohFull(k) = mean(coh_full(validMask));
    CohBit(k) = mean(coh_bit(validMask));

    %%%%%%%%%%%%%%%%%%%%% 9.8 相位误差
    e_before = angle(exp(1j*( ...
        phiDef_bit-phiDef_full)));

    e_after = angle(exp(1j*( ...
        phiFilt_bit-phiFilt_full)));

    PhaseRMSE_before(k) = CircularRMSE( ...
        e_before(validMask));

    PhaseRMSE_after(k) = CircularRMSE( ...
        e_after(validMask));

    %%%%%%%%%%%%%%%%%%%%% 9.9 参考相位
    phi_ref_full = median( ...
        phiFilt_full(ref_area),'omitnan');

    phi_ref_bit = median( ...
        phiFilt_bit(ref_area),'omitnan');

    %%%%%%%%%%%%%%%%%%%%% 9.10 高程反演
    h_full = Z0+(phiFilt_full-phi_ref_full)./sens;
    h_bit = Z0+(phiFilt_bit-phi_ref_bit)./sens;

    %%%%%%%%%%%%%%%%%%%%% 9.11 高程误差
    err_full = h_full-h_true;
    err_bit = h_bit-h_true;

    metrics_full = HeightMetrics(err_full,validMask);
    metrics_bit = HeightMetrics(err_bit,validMask);

    HeightRMSE_full(k) = metrics_full.RMSE;
    HeightRMSE_bit(k) = metrics_bit.RMSE;

    HeightMAE_full(k) = metrics_full.MAE;
    HeightMAE_bit(k) = metrics_bit.MAE;

    HeightP95_full(k) = metrics_full.P95;
    HeightP95_bit(k) = metrics_bit.P95;

    Outlier_full(k) = metrics_full.Outlier2m;
    Outlier_bit(k) = metrics_bit.Outlier2m;

    %%%%%%%%%%%%%%%%%%%%% 9.12 分区域高程误差
    hillMask = validMask & ...
        h_true>=BaseHeight+0.5*HeightHill;

    flatMask = validMask & ...
        h_true<=BaseHeight+0.1*HeightHill;

    HillRMSE_full(k) = sqrt( ...
        mean(err_full(hillMask).^2));

    HillRMSE_bit(k) = sqrt( ...
        mean(err_bit(hillMask).^2));

    FlatRMSE_full(k) = sqrt( ...
        mean(err_full(flatMask).^2));

    FlatRMSE_bit(k) = sqrt( ...
        mean(err_bit(flatMask).^2));

    %%%%%%%%%%%%%%%%%%%%% 9.13 保存典型结果
    showIdx = find(showSNR==snr_dB,1);

    if ~isempty(showIdx)
        h_full_show{showIdx} = h_full;
        h_bit_show{showIdx} = h_bit;
    end

    %%%%%%%%%%%%%%%%%%%%% 9.14 输出统计
    fprintf('RDA幅度相关A = %.4f\n',AmpCorr(k,1));
    fprintf('RDA幅度相关B = %.4f\n',AmpCorr(k,2));

    fprintf('全精度相干系数 = %.4f\n',CohFull(k));
    fprintf('单比特相干系数 = %.4f\n',CohBit(k));

    fprintf('滤波后相位RMSE = %.5f rad\n', ...
        PhaseRMSE_after(k));

    fprintf('全精度高程RMSE = %.5f m\n', ...
        HeightRMSE_full(k));

    fprintf('单比特高程RMSE = %.5f m\n', ...
        HeightRMSE_bit(k));

end

%%%%%%%%%%%%%%%%%%%%%%%%%% 10. RDA相关性
figure('Color','w','Position',[100 100 950 480]);

plot(SNR_plot,mean(AmpCorr,2),'o-', ...
    'LineWidth',1.6,'MarkerSize',7);
hold on;

plot(SNR_plot,mean(ComplexCorr,2),'s--', ...
    'LineWidth',1.6,'MarkerSize',7);

SetSNRAxis(SNR_plot,SNR_labels);

ylabel('相关系数');
title('不同SNR下的RDA成像相关性');

legend('幅度相关系数','复数相关系数', ...
    'Location','best');

ylim([0 1.05]);
grid on;

exportgraphics(gcf,'图1_SNR_RDA相关性.png', ...
    'Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%% 11. 背景能量占比
figure('Color','w','Position',[110 110 950 480]);

plot(SNR_plot,BackgroundFull,'o-', ...
    'LineWidth',1.6,'MarkerSize',7);
hold on;

plot(SNR_plot,BackgroundBit,'s--', ...
    'LineWidth',1.6,'MarkerSize',7);

SetSNRAxis(SNR_plot,SNR_labels);

ylabel('背景能量占比/%');
title('不同SNR下的RDA背景能量');

legend('含噪全精度','零阈值单比特', ...
    'Location','best');

grid on;

exportgraphics(gcf,'图2_SNR_背景能量.png', ...
    'Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%% 12. 干涉相干系数
figure('Color','w','Position',[120 120 950 480]);

plot(SNR_plot,CohFull,'o-', ...
    'LineWidth',1.6,'MarkerSize',7);
hold on;

plot(SNR_plot,CohBit,'s--', ...
    'LineWidth',1.6,'MarkerSize',7);

SetSNRAxis(SNR_plot,SNR_labels);

ylabel('平均局部干涉相干系数');
title('不同SNR下的双通道干涉相干性');

legend('含噪全精度','零阈值单比特', ...
    'Location','best');

ylim([0 1.05]);
grid on;

exportgraphics(gcf,'图3_SNR_干涉相干性.png', ...
    'Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%% 13. 相位RMSE
figure('Color','w','Position',[130 130 950 480]);

plot(SNR_plot,PhaseRMSE_before,'o-', ...
    'LineWidth',1.6,'MarkerSize',7);
hold on;

plot(SNR_plot,PhaseRMSE_after,'s--', ...
    'LineWidth',1.6,'MarkerSize',7);

SetSNRAxis(SNR_plot,SNR_labels);

ylabel('相位RMSE/rad');
title('不同SNR下的单比特干涉相位误差');

legend('滤波前','滤波后', ...
    'Location','best');

grid on;

exportgraphics(gcf,'图4_SNR_相位RMSE.png', ...
    'Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%% 14. 高程RMSE和MAE
figure('Color','w','Position',[140 100 1200 480]);

subplot(1,2,1);

plot(SNR_plot,HeightRMSE_full,'o-', ...
    'LineWidth',1.6,'MarkerSize',7);
hold on;

plot(SNR_plot,HeightRMSE_bit,'s--', ...
    'LineWidth',1.6,'MarkerSize',7);

SetSNRAxis(SNR_plot,SNR_labels);

ylabel('高程RMSE/m');
title('不同SNR下的高程RMSE');

legend('含噪全精度','零阈值单比特', ...
    'Location','best');

grid on;

subplot(1,2,2);

plot(SNR_plot,HeightMAE_full,'o-', ...
    'LineWidth',1.6,'MarkerSize',7);
hold on;

plot(SNR_plot,HeightMAE_bit,'s--', ...
    'LineWidth',1.6,'MarkerSize',7);

SetSNRAxis(SNR_plot,SNR_labels);

ylabel('高程MAE/m');
title('不同SNR下的高程MAE');

legend('含噪全精度','零阈值单比特', ...
    'Location','best');

grid on;

exportgraphics(gcf,'图5_SNR_高程精度.png', ...
    'Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%% 15. 典型SNR高程剖面
[~,centerCol] = min(abs(Ax_r-xc));

figure('Color','w','Position',[150 50 1150 950]);

for k = 1:Nshow

    subplot(Nshow,1,k);

    plot(Y_ground,h_true(:,centerCol), ...
        'k-','LineWidth',1.8);
    hold on;

    plot(Y_ground,h_full_show{k}(:,centerCol), ...
        'LineWidth',1.3);

    plot(Y_ground,h_bit_show{k}(:,centerCol),'--', ...
        'LineWidth',1.2);

    if isinf(showSNR(k))
        title('无噪声');
    else
        title(sprintf('SNR=%.0fdB',showSNR(k)));
    end

    xlabel('距离向/m');
    ylabel('高程/m');

    legend('真实高程','含噪全精度', ...
        '零阈值单比特','Location','best');

    grid on;

end

sgtitle('不同SNR下的山丘中心高程剖面');

exportgraphics(gcf,'图6_SNR_高程剖面.png', ...
    'Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%% 16. 典型高程二维对比
figure('Color','w','Position',[80 50 1300 1050]);

h_limits = [BaseHeight-2, ...
            BaseHeight+HeightHill+2];

for k = 1:Nshow

    subplot(Nshow,3,3*k-2);
    PlotHeight(Ax_r,Y_ground,h_true,h_limits);

    if isinf(showSNR(k))
        title('无噪声：真实高程');
    else
        title(sprintf('SNR=%.0fdB：真实高程',showSNR(k)));
    end

    subplot(Nshow,3,3*k-1);
    PlotHeight(Ax_r,Y_ground,h_full_show{k},h_limits);
    title('全精度高程');

    subplot(Nshow,3,3*k);
    PlotHeight(Ax_r,Y_ground,h_bit_show{k},h_limits);
    title('零阈值单比特高程');

end

sgtitle('不同SNR下的高程反演对比');

exportgraphics(gcf,'图7_SNR_高程二维对比.png', ...
    'Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%% 17. 分区域高程RMSE
figure('Color','w','Position',[160 120 1050 470]);

subplot(1,2,1);

plot(SNR_plot,HillRMSE_full,'o-', ...
    'LineWidth',1.5);
hold on;

plot(SNR_plot,HillRMSE_bit,'s--', ...
    'LineWidth',1.5);

SetSNRAxis(SNR_plot,SNR_labels);

ylabel('高程RMSE/m');
title('山丘区域RMSE');

legend('含噪全精度','零阈值单比特');
grid on;

subplot(1,2,2);

plot(SNR_plot,FlatRMSE_full,'o-', ...
    'LineWidth',1.5);
hold on;

plot(SNR_plot,FlatRMSE_bit,'s--', ...
    'LineWidth',1.5);

SetSNRAxis(SNR_plot,SNR_labels);

ylabel('高程RMSE/m');
title('平坦区域RMSE');

legend('含噪全精度','零阈值单比特');
grid on;

exportgraphics(gcf,'图8_SNR_分区域高程RMSE.png', ...
    'Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%% 18. 汇总数据
Results = table( ...
    SNR_list(:), ...
    mean(AmpCorr,2), ...
    mean(ComplexCorr,2), ...
    BackgroundFull, ...
    BackgroundBit, ...
    CohFull, ...
    CohBit, ...
    PhaseRMSE_before, ...
    PhaseRMSE_after, ...
    HeightRMSE_full, ...
    HeightRMSE_bit, ...
    HeightMAE_full, ...
    HeightMAE_bit, ...
    HeightP95_full, ...
    HeightP95_bit, ...
    Outlier_full, ...
    Outlier_bit, ...
    HillRMSE_full, ...
    HillRMSE_bit, ...
    FlatRMSE_full, ...
    FlatRMSE_bit, ...
    'VariableNames',{ ...
    'SNR_dB', ...
    'AmplitudeCorr', ...
    'ComplexCorr', ...
    'BackgroundFull_pct', ...
    'BackgroundBit_pct', ...
    'CohFull', ...
    'CohBit', ...
    'PhaseRMSE_Before', ...
    'PhaseRMSE_After', ...
    'HeightRMSE_Full', ...
    'HeightRMSE_Bit', ...
    'HeightMAE_Full', ...
    'HeightMAE_Bit', ...
    'HeightP95_Full', ...
    'HeightP95_Bit', ...
    'OutlierFull_pct', ...
    'OutlierBit_pct', ...
    'HillRMSE_Full', ...
    'HillRMSE_Bit', ...
    'FlatRMSE_Full', ...
    'FlatRMSE_Bit'});

disp('========== SNR实验统计 ==========');
disp(Results);

writetable(Results,'OneBit_InSAR_SNR_Results.csv');

fprintf('\n全部SNR及无噪声实验完成。\n');

%%%%%%%%%%%%%%%%%%%%%%%%%% 局部函数

function Data = QuantizeZero(Data)

Data = complex( ...
    2*double(real(Data)>=0)-1, ...
    2*double(imag(Data)>=0)-1);

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

Y_full = ifft( ...
    fft(Echo,Nfft,1).*fft(h_mf,Nfft,1),Nfft,1);

startIdx = floor(M/2)+1;

HRRP = Y_full(startIdx:startIdx+nrn-1,:);

end

function [img,x_rda] = RunRDA( ...
    Echo,Rr,B,Tp,C,PRF,V,lambda,Nfft_az)

%%%%%%%%%%%%%%%%% 1. 距离脉压
HRRP = RangeCompress(Echo,Rr,B,Tp,C);

[nrn,naz] = size(HRRP);

%%%%%%%%%%%%%%%%% 2. 方位加窗
az_win = hann(naz,'periodic').';

HRRP = HRRP.*az_win;

%%%%%%%%%%%%%%%%% 3. 方位补零
A_pad = complex(zeros(nrn,Nfft_az));

i0 = floor((Nfft_az-naz)/2);

A_pad(:,i0+1:i0+naz) = HRRP;

%%%%%%%%%%%%%%%%% 4. 方位FFT
fa = (-Nfft_az/2:Nfft_az/2-1)*PRF/Nfft_az;

S = fftshift( ...
    fft(ifftshift(A_pad,2),Nfft_az,2),2);

clear A_pad HRRP

%%%%%%%%%%%%%%%%% 5. RCMC
S_rcmc = complex(zeros(size(S)));

for k = 1:Nfft_az

    D2 = 1-(lambda*fa(k)/(2*V))^2;

    if D2<=0
        continue;
    end

    D = sqrt(D2);

    R_query = Rr/D;

    S_rcmc(:,k) = interp1( ...
        Rr,S(:,k),R_query,'linear',0);

end

clear S

%%%%%%%%%%%%%%%%% 6. 方位匹配滤波
for r = 1:nrn

    Ka = 2*V^2/(lambda*Rr(r));

    H_az = exp(-1j*pi*fa.^2/Ka);

    S_rcmc(r,:) = S_rcmc(r,:).*H_az;

end

%%%%%%%%%%%%%%%%% 7. 方位IFFT
img = fftshift( ...
    ifft(ifftshift(S_rcmc,2),Nfft_az,2),2);

x_rda = V*((-Nfft_az/2:Nfft_az/2-1)/PRF);

end

function [phi_deflat,phi_filt,coh] = ...
    ProcessInSAR(A,B,phi_flat,win)

%%%%%%%%%%%%%%%%% 1. 干涉及去平地
I_reg = A.*conj(B);

I_deflat = I_reg.*exp(-1j*phi_flat);

phi_deflat = angle(I_deflat);

%%%%%%%%%%%%%%%%% 2. 复数回转均值滤波
I_unit = I_deflat./(abs(I_deflat)+eps);

I_filt = conv2(I_unit,win,'same');

phi_filt = angle(I_filt);

%%%%%%%%%%%%%%%%% 3. 局部干涉相干系数
num = abs(conv2(A.*conj(B),win,'same'));

denA = conv2(abs(A).^2,win,'same');
denB = conv2(abs(B).^2,win,'same');

coh = num./sqrt(denA.*denB+eps);

end

function rho = GetAmplitudeCorr(A,B)

a = abs(A(:));
b = abs(B(:));

R = corrcoef(a,b);

rho = R(1,2);

end

function gamma = GetComplexCorr(A,B)

a = A(:);
b = B(:);

gamma = abs(sum(a.*conj(b)))/ ...
    sqrt(sum(abs(a).^2)*sum(abs(b).^2));

end

function rmse = CircularRMSE(e)

e = e(isfinite(e));

bias = angle(mean(exp(1j*e)));

e = angle(exp(1j*(e-bias)));

rmse = sqrt(mean(e.^2));

end

function s = HeightMetrics(err,mask)

e = err(mask);
e = e(isfinite(e));

s.RMSE = sqrt(mean(e.^2));
s.MAE = mean(abs(e));
s.P95 = prctile(abs(e),95);
s.Outlier2m = 100*mean(abs(e)>2);

end

function SetSNRAxis(x,labels)

set(gca,'XTick',x,'XTickLabel',labels);
xlim([min(x)-0.3 max(x)+0.3]);
xlabel('输入参考SNR/dB');

end

function PlotHeight(Ax,Yg,Hmap,limits)

imagesc(Ax,Yg,Hmap);

axis xy;
clim(limits);
colormap(parula);
colorbar;

xlabel('方位向/m');
ylabel('距离向/m');

end
