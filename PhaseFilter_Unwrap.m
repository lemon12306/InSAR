clear; close all; clc;

load('RDA_Zero_OneBit.mat');
SAR_Parameters;

imgA = RDAimage(:,:,1);
imgB = RDAimage(:,:,2);

Ax = x_rda;
Ay = Rr;

[nr, nc] = size(imgA);

I_raw = imgA .* conj(imgB);
phi_raw = angle(I_raw);

figure;
imagesc(Ax, Ay, phi_raw);
title('配准前干涉相位');xlabel('方位向 (m)');ylabel('斜距 (m)');

%%%%%%%%%% 粗配准

A = abs(imgA);
B = abs(imgB);

A0 = A - mean(A(:));
B0 = B - mean(B(:));

FA = fft2(A0);
FB = fft2(B0);

C = ifft2(FB .* conj(FA));

[~, ind] = max(abs(C(:)));
[rpk, cpk] = ind2sub(size(C), ind);

dr = rpk - 1;
dc = cpk - 1;

if dr > nr/2
    dr = dr - nr;
end
if dc > nc/2
    dc = dc - nc;
end

imgB_reg = circshift(imgB, [-dr, -dc]);

%%%%%%%%%% 裁剪有效区域

ROW_min = find(Ay >= 40, 1, 'first');
ROW_max = find(Ay <= 60, 1, 'last');

COL_min = find(Ax >= -10, 1, 'first');
COL_max = find(Ax <= 10, 1, 'last');

dx_insar = 0.05;

dx_rda = abs(Ax(2) - Ax(1));

col_step = max(1, round(dx_insar / dx_rda));

COL_idx = COL_min:col_step:COL_max;
ROW_idx = ROW_min:ROW_max;

imgA_cut = imgA(ROW_idx, COL_idx);
imgB_cut = imgB_reg(ROW_idx, COL_idx);

Ax_r = Ax(COL_idx);
Ay_r = Ay(ROW_idx);

[nr_cut, nc_cut] = size(imgA_cut);

I_reg = imgA_cut .* conj(imgB_cut);

figure;
imagesc(Ax_r, Ay_r, angle(I_reg));
title('配准后干涉相位');xlabel('方位向 (m)');ylabel('斜距 (m)');

%%%%%%%%%% 去平地相位

zA = Antenna10(3);
zB = Antenna20(3);

yA = Antenna10(2);
yB = Antenna20(2);

Y_ground = sqrt(max(Ay_r.^2 - (Z0-zA)^2, 0));

[~, gy] = meshgrid(Ax_r, Y_ground);

RrxA_Z0 = sqrt((gy-yA).^2 + (Z0-zA).^2);
RrxB_Z0 = sqrt((gy-yB).^2 + (Z0-zB).^2);

phi_flat = (2*pi/lambda) .* (RrxB_Z0 - RrxA_Z0);

I_deflat = I_reg .* exp(-1j * phi_flat);
phi_deflat = angle(I_deflat);

figure;
imagesc(Ax_r, Ay_r, phi_deflat);
title('去平地后干涉相位');xlabel('方位向 (m)');ylabel('斜距 (m)');

%%%%%%%%%% 回转均值滤波

% win = ones(5)/25;
% 
% I_filt = conv2(I_deflat, win, 'same');
% 
% phi_filt = angle(I_filt);
% 
% figure;
% imagesc(Ax_r, Ay_r, phi_filt);axis xy; clim([-pi pi]);
% title('回转均值滤波');xlabel('方位向 (m)');ylabel('斜距 (m)');
%%%%%%%%%% 回转均值滤波

win = ones(5)/25;

I_filt = conv2(I_deflat, win, 'same');

phi_filt = angle(I_filt);

figure;
imagesc(Ax_r, Ay_r, phi_filt);axis xy; clim([-pi pi]); colorbar;
title('回转均值滤波');xlabel('方位向 (m)');ylabel('斜距 (m)');

%%%%%%%%%%%%%%%%二维相位解缠绕

ampA = abs(imgA_cut);
ampB = abs(imgB_cut);

coh_win = ones(5);

num = abs(conv2(imgA_cut .* conj(imgB_cut), coh_win, 'same'));

denA = conv2(ampA.^2, coh_win, 'same');
denB = conv2(ampB.^2, coh_win, 'same');

coh = num ./ sqrt(denA .* denB + eps);

amp = sqrt(ampA .* ampB);%几何平均

amp_norm = amp ./ (max(amp(:)) + eps);

quality = coh .* sqrt(amp_norm);%减弱幅度权重

quality = quality ./ (max(quality(:)) + eps);



%phi_unw_deflat = QualityGuidedUnwrap2D(phi_filt, quality);
phi_unw_deflat = LeastSquaresUnwrap2D(phi_filt);


figure;
imagesc(Ax_r, Ay_r, phi_unw_deflat);
title('二维解缠相位');xlabel('方位向 (m)');ylabel('斜距 (m)');


% %%%%%%%%% 高程反演
% 
% sens = (2*pi/lambda) .* ((Z0-zB)./RrxB_Z0 - (Z0-zA)./RrxA_Z0);% 高程灵敏度
% 
% ref_mask_x = Ax_r >= -10 & Ax_r <= -8;
% ref_mask_y = Ay_r >= 40 & Ay_r <= 43;
% ref_area = ref_mask_y(:) & ref_mask_x(:).';
% 
% phi_ref = median(phi_unw_deflat(ref_area), 'omitnan');
% 
% h_est = Z0 + (phi_unw_deflat - phi_ref) ./ sens;
% 
% % phi_ref = median(phi_filt(ref_area), 'omitnan');
% % 
% % h_est = Z0 + (phi_filt - phi_ref) ./ sens;
% 
% 
% figure;
% surf(Ax_r, Y_ground, h_est, 'EdgeColor', 'none');
% title('高程反演3D');xlabel('方位向 (m)');ylabel('距离向 (m)');zlabel('高度 (m)');

%%%%%%%%%% 高程反演

sens = (2*pi/lambda) .* ...
    ((Z0-zB)./RrxB_Z0 - (Z0-zA)./RrxA_Z0);

ref_mask_x = Ax_r >= -10 & Ax_r <= -8;
ref_mask_y = Ay_r >= 40 & Ay_r <= 43;

ref_area = ref_mask_y(:) & ref_mask_x(:).';

phi_ref = median(phi_filt(ref_area), 'omitnan');

h_est = Z0 + (phi_filt - phi_ref) ./ sens;

figure;
surf(Ax_r, Y_ground, h_est, 'EdgeColor', 'none');
title('高程反演3D');xlabel('方位向 (m)');ylabel('距离向 (m)');zlabel('高度 (m)');


function phi_unw = LeastSquaresUnwrap2D(phi_wrapped)

[nr, nc] = size(phi_wrapped);

gx = phi_wrapped(:,2:end) - phi_wrapped(:,1:end-1);% 计算x方向缠绕相位梯度
gx = angle(exp(1j * gx));

gy = phi_wrapped(2:end,:) - phi_wrapped(1:end-1,:);% 计算y方向缠绕相位梯度
gy = angle(exp(1j * gy));

rho = zeros(nr,nc);% 构造最小二乘Poisson方程右端项

rho(:,1) = rho(:,1) - gx(:,1);% x方向散度

if nc > 2

    rho(:,2:nc-1) = rho(:,2:nc-1) + gx(:,1:end-1) - gx(:,2:end);

end

rho(:,nc) = rho(:,nc) + gx(:,end);
rho(1,:) = rho(1,:) - gy(1,:);% y方向散度

if nr > 2

    rho(2:nr-1,:) = rho(2:nr-1,:) + gy(1:end-1,:) - gy(2:end,:);

end

rho(nr,:) = rho(nr,:) + gy(end,:);
rho_dct = dct2(rho);% DCT求解二维Poisson方程
[col_idx, row_idx] = meshgrid(0:nc-1, 0:nr-1);

den = 2 - 2*cos(pi*col_idx/nc) + 2 - 2*cos(pi*row_idx/nr);
den(1,1) = 1;% 去除直流奇异点
phi_dct = rho_dct ./ den;
phi_dct(1,1) = 0;
phi_unw = idct2(phi_dct);% 恢复二维连续相位

end
