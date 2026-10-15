
clear; close all; clc;
rng(20261011);
SAR_Parameters;

ThetaTa = deg2rad(55);
Tp = 30e-9;
SNR_dB = 30;

outDir = 'RDA_Results';
if ~exist(outDir,'dir')
    mkdir(outDir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%% 1. 生成原始回波
[Echo,~] = RawDataGen_rand( ...
    point,Rr,nan,B,fc,PRF,Ym,H,V,ThetaTa,Bz,0,Tp,'raw');

Echo_A = Echo(:,:,1);
Echo_B = Echo(:,:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%% 2. 添加接收噪声
maskA = Echo_A~=0;
maskB = Echo_B~=0;

PowerA = mean(abs(Echo_A(maskA)).^2);
PowerB = mean(abs(Echo_B(maskB)).^2);
P_ref = (PowerA+PowerB)/2;

NoisePower = P_ref/10^(SNR_dB/10);

Noise_A = sqrt(NoisePower/2)*( ...
    randn(size(Echo_A))+1j*randn(size(Echo_A)));

Noise_B = sqrt(NoisePower/2)*( ...
    randn(size(Echo_B))+1j*randn(size(Echo_B)));

Input_A = Echo_A+Noise_A;
Input_B = Echo_B+Noise_B;

%%%%%%%%%%%%%%%%%%%%%%%%%% 3. 零阈值单比特量化
Input_A_1bit = complex( ...
    2*double(real(Input_A)>=0)-1, ...
    2*double(imag(Input_A)>=0)-1);

Input_B_1bit = complex( ...
    2*double(real(Input_B)>=0)-1, ...
    2*double(imag(Input_B)>=0)-1);

%%%%%%%%%%%%%%%%%%%%%%%%%% 4. 双模式RDA成像
[imgA_full,imgB_full,x_rda] = RunRDA( ...
    Input_A,Input_B,Rr,nan,B,Tp,C,PRF,V,lambda);

[imgA_1bit,imgB_1bit,~] = RunRDA( ...
    Input_A_1bit,Input_B_1bit,Rr,nan,B,Tp,C,PRF,V,lambda);

%%%%%%%%%%%%%%%%%%%%%%%%%% 5. 成像结果对比
peakFull = max([abs(imgA_full(:));abs(imgB_full(:))]);
peakBit = max([abs(imgA_1bit(:));abs(imgB_1bit(:))]);

A_full_dB = 20*log10(abs(imgA_full)/peakFull+eps);
B_full_dB = 20*log10(abs(imgB_full)/peakFull+eps);
A_bit_dB = 20*log10(abs(imgA_1bit)/peakBit+eps);
B_bit_dB = 20*log10(abs(imgB_1bit)/peakBit+eps);

figure('Color','w','Position',[100 100 1200 800]);

subplot(2,2,1);
PlotRDA(x_rda,Rr,A_full_dB);
title('含噪全精度RDA - 通道A');

subplot(2,2,2);
PlotRDA(x_rda,Rr,B_full_dB);
title('含噪全精度RDA - 通道B');

subplot(2,2,3);
PlotRDA(x_rda,Rr,A_bit_dB);
title('零阈值单比特RDA - 通道A');

subplot(2,2,4);
PlotRDA(x_rda,Rr,B_bit_dB);
title('零阈值单比特RDA - 通道B');



%%%%%%%%%%%%%%%%%%%%%%%%%% 6. 成像相关性
ampCorrA = GetAmplitudeCorr(imgA_full,imgA_1bit);
ampCorrB = GetAmplitudeCorr(imgB_full,imgB_1bit);

complexCorrA = GetComplexCorr(imgA_full,imgA_1bit);
complexCorrB = GetComplexCorr(imgB_full,imgB_1bit);

fprintf('\n========== RDA成像结果 ==========\n');
fprintf('输入参考SNR = %.2f dB\n',SNR_dB);
fprintf('通道A幅度相关系数 = %.6f\n',ampCorrA);
fprintf('通道B幅度相关系数 = %.6f\n',ampCorrB);
fprintf('通道A复数相关系数 = %.6f\n',complexCorrA);
fprintf('通道B复数相关系数 = %.6f\n',complexCorrB);

figure('Color','w','Position',[200 150 850 500]);

metrics = [ampCorrA complexCorrA; ...
           ampCorrB complexCorrB];

b = bar(metrics,'grouped','BarWidth',0.65);
b(1).FaceColor = [0.10 0.62 0.72];
b(2).FaceColor = [0.94 0.51 0.44];
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';

hold on;
for j = 1:2
    for i = 1:2
        text(b(j).XEndPoints(i),metrics(i,j)+0.02, ...
            sprintf('%.4f',metrics(i,j)), ...
            'HorizontalAlignment','center', ...
            'FontName','Times New Roman','FontSize',11);
    end
end

set(gca,'XTick',1:2,'XTickLabel',{'通道A','通道B'});
ylabel('相关系数');
title('RDA成像相关性对比');
legend('幅度相关系数','复数相关系数','Location','northwest');
ylim([0 1.1]);
grid on;



%%%%%%%%%%%%%%%%%%%%%%%%%% 7. 保存结果

% 全精度RDA
RDAimage = cat(3,imgA_full,imgB_full);
save RDA_Full_Noisy.mat

% 零阈值单比特RDA
RDAimage = cat(3,imgA_1bit,imgB_1bit);
save RDA_Zero_OneBit.mat



%%%%%%%%%%%%%%%%%%%%%%%%%% 局部函数
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

function [imgA,imgB,x_rda] = RunRDA( ...
    Echo_A,Echo_B,Rr,nan,B,Tp,C,PRF,V,lambda)

%%%%%%%%%%%%%%%%% 1. 距离脉压
HRRP_A = RangeCompress(Echo_A,Rr,B,Tp,C);
HRRP_B = RangeCompress(Echo_B,Rr,B,Tp,C);

%%%%%%%%%%%%%%%%% 2. 方位FFT
nrn = size(Echo_A,1);
x_az_max = 12;

Nfft_need = ceil(2*x_az_max*PRF/V);
Nfft_az = 2^nextpow2(max(Nfft_need,nan));

fa = (-Nfft_az/2:Nfft_az/2-1)*PRF/Nfft_az;
az_win = hann(nan,'periodic').';

HRRP_A_w = HRRP_A.*az_win;
HRRP_B_w = HRRP_B.*az_win;

A_pad = complex(zeros(nrn,Nfft_az));
B_pad = complex(zeros(nrn,Nfft_az));

i0 = floor((Nfft_az-nan)/2);
A_pad(:,i0+1:i0+nan) = HRRP_A_w;
B_pad(:,i0+1:i0+nan) = HRRP_B_w;

S_A = fftshift(fft(ifftshift(A_pad,2),Nfft_az,2),2);
S_B = fftshift(fft(ifftshift(B_pad,2),Nfft_az,2),2);

%%%%%%%%%%%%%%%%% 3. RCMC
S_rcmc_A = complex(zeros(size(S_A)));
S_rcmc_B = complex(zeros(size(S_B)));

for k = 1:Nfft_az
    D2 = 1-(lambda*fa(k)/(2*V))^2;
    if D2<=0
        continue;
    end
    D = sqrt(D2);
    R_query = Rr/D;

    S_rcmc_A(:,k) = interp1(Rr,S_A(:,k),R_query,'linear',0);
    S_rcmc_B(:,k) = interp1(Rr,S_B(:,k),R_query,'linear',0);
end

%%%%%%%%%%%%%%%%% 4. 方位匹配滤波
for r = 1:nrn
    Ka = 2*V^2/(lambda*Rr(r));
    H_az = exp(-1j*pi*fa.^2/Ka);

    S_rcmc_A(r,:) = S_rcmc_A(r,:).*H_az;
    S_rcmc_B(r,:) = S_rcmc_B(r,:).*H_az;
end

%%%%%%%%%%%%%%%%% 5. 方位IFFT
imgA = fftshift(ifft(ifftshift(S_rcmc_A,2),Nfft_az,2),2);
imgB = fftshift(ifft(ifftshift(S_rcmc_B,2),Nfft_az,2),2);

x_rda = V*((-Nfft_az/2:Nfft_az/2-1)/PRF);

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
gamma = abs(sum(a.*conj(b))) / ...
    sqrt(sum(abs(a).^2)*sum(abs(b).^2));
end

function PlotRDA(x_rda,Rr,img_dB)
imagesc(x_rda,Rr,img_dB);
axis xy;
xlim([-20 20]);
ylim([min(Rr) max(Rr)]);
clim([-14 0]);
colormap(parula);
colorbar;
xlabel('方位向/m');
ylabel('等效斜距/m');
end
