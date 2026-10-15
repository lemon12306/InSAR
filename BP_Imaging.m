clear; close all; clc;
SAR_Parameters; 
%BP 成像网格 
ThetaTa = deg2rad(55); % 方位向波束宽度
BW = rad2deg(ThetaTa); % 波束宽度（度）
Rho = 0.25; % BP 成像像素间隔
Ax = -15:Rho:15;% 方位向x成像坐标
Ay = (35:Rho:65).'; % 距离向y成像坐标
NX = length(Ax);
NY = length(Ay);
[AX, AY] = meshgrid(Ax, Ay);% 成像网格，NY × NX
% 波束掩膜，只保留波束覆盖方位角范围内的像素
az_angle = atan2(AX, AY);
Mask = (AY > 0) & (abs(az_angle) <= ThetaTa/2) & ...
       (sqrt(AX.^2 + AY.^2) <= Rmax);
%生成双通道理想距离压缩回波HRRP
[HRRP, simInfo] = RawDataGen_rand(point, Rr, nan, B, fc, PRF, Ym, H, V, ThetaTa, Bz);
figure;imagesc(1:nan, Rr, abs(HRRP(:,:,1)));
xlabel('脉冲序号'); ylabel('等效斜距R(m)'); title('通道A HRRP幅度');
figure;imagesc(1:nan, Rr, abs(HRRP(:,:,2)));
xlabel('脉冲序号'); ylabel('等效斜距R(m)'); title('通道B HRRP幅度');
%双通道BP成像
BPimage = complex(zeros(NY, NX, nchannel));% 对每个通道、脉冲，将回波能量按几何关系反向投影到成像网格
Zref = 0;% BP 参考高度平面
T_y = 0; % 发射天线y坐标
Z = [H, H+Bz];% 两接收天线z坐标
for nc = 1:nchannel
    for np = 1:nan
        xplat = V*ta(np); % 当前脉冲平台x位置
        Rtx = sqrt((AX - xplat).^2 + (AY - T_y).^2 + (Zref - H).^2);% 发射天线到每个成像像素的斜距
        Rrx = sqrt((AX - xplat).^2 + (AY - Ym(nc)).^2 + (Zref - Z(nc)).^2);% 接收天线到每个成像像素的斜距
        Rbp = 0.5*(Rtx + Rrx);% 单发双收等效斜距
        range_line = HRRP(:,np,nc);
        BPimage(:,:,nc) = BPimage(:,:,nc) + ...
            interp1(Rr, range_line, Rbp, 'linear', 0) .* ...
            exp(1j*4*pi/lambda * Rbp) .* Mask;
    end
end
imgA = BPimage(:,:,1); % 通道A复图像
imgB = BPimage(:,:,2); % 通道B复图像
imgB_shift = [2, 3];% 模拟辅图配准误差
imgB = circshift(imgB, imgB_shift);
% 显示 BP 幅度图像
figure; imagesc(abs(imgA)); 
title('通道A BP'); xlabel('X方位向(m)'); ylabel('Y侧视距离向(m)');
figure; imagesc(abs(imgB));
title('通道B BP'); xlabel('X方位向(m)'); ylabel('Y侧视距离向(m)');
save BP_complex_img.mat