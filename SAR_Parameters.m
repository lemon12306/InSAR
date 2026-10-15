% 坐标系定义：
% x：车辆运动方向/方位向
% y：距离向/侧视方向
% z：高度方向
% 侧视几何：车辆沿x方向运动，目标位于车辆侧方y=40~60 m
C = 3e8;% 光速
fc = 77e9;% 载频，77 GHz 毫米波
lambda = C/fc;% 波长
B = 300e6;% 距离向带宽
PRF = 2000;% 脉冲重复频率
V = 5;% 车速
deltaR = C/(2*B);% 距离分辨率
Rmax = 100;% 最大有效侧视距离
dR = 0.1;% HRRP距离采样间隔
Rr = (32:dR:70).'; % 等效斜距轴
nrn = length(Rr);% 距离采样点数
nan = 512; % 方位脉冲数
ta = (-(nan-1)/2:(nan-1)/2)/PRF; % 慢时间轴
H = 0.66; % 天线距地面高度

B_y = 0; % y向基线（两接收天线间距）
Bz = 0.00585;% z向基线（高度差）
x_min = -10; x_max = 10; % 目标区域x方向范围
y_min =  40; y_max = 60; % 目标区域y方向范围
Y0 = 50; % 目标区域中心侧视距离
Z0 = 0.1;% 地表高度
R0 = sqrt(Y0^2 + (H-Z0)^2); % 场景中心参考斜距
Ym = [-B_y/2, B_y/2];% 两个接收天线的y坐标
% Zm = [0, Bz];% 两个接收天线的z向相对高度

% 天线初始位置
AntennaT0 = [0,0,H];% 发射天线
Antenna10 = [0, -B_y/2, H+Bz];% 接收通道A
Antenna20 = [0,  B_y/2, H+Bz+Bz];% 接收通道B
nchannel = 2;% 通道数

%构建高斯丘
xc = 0;  yc = 50;% 中心坐标
dx = 0.10;  dy = 0.10;  % 模型网格采样间隔
x_v = x_min:dx:x_max;
y_v = y_min:dy:y_max;
[X_g, Y_g] = meshgrid(x_v, y_v);  % 41×41 模型网格
% 二维高斯函数
BaseHeight = 0.1;% 基底高度
HeightHill = 5;% 丘峰值高度
s_x = 2.8;  sigma_y = 2.8;% 丘宽度
Hill = exp(-0.5*((X_g-xc)/s_x).^2 - ...
           0.5*((Y_g-yc)/sigma_y).^2);
Z_g = BaseHeight + HeightHill * Hill;  % 真实高程
A_cone = ones(size(Z_g));  % 均匀散射幅度

% 加大抖动破坏规则网格 + 随机相位破坏相干叠加
% jitter = 0;
% Xp = X_g(:) + jitter*randn(numel(X_g),1);
% Yp = Y_g(:) + jitter*randn(numel(Y_g),1);
% Zp = BaseHeight * ones(numel(Z_g),1);
% Zp = Z_g(:);
% %Amp = A_cone(:) .* exp(1j*2*pi*rand(numel(Z_g),1));
% Amp = A_cone(:);
jitter = 0.02;

Xp = X_g(:) + jitter*randn(numel(X_g),1);
Yp = Y_g(:) + jitter*randn(numel(Y_g),1);
Zp = BaseHeight + HeightHill .* exp(-0.5*((Xp-xc)/s_x).^2 - 0.5*((Yp-yc)/sigma_y).^2);

% 随机正幅度，均值约为1；不引入随机相位，便于单独考察幅度影响
Amp = A_cone(:) .* (0.5 + rand(numel(Z_g),1));
Vx = zeros(numel(Zp),1);% x向速度（静止目标）
Vy = zeros(numel(Zp),1);% y向速度（静止目标）
point = [Xp, Yp, Zp, Vx, Vy, Amp];% 散射点列表，1681 点


% figure; 
% mesh(X_g, Y_g, Z_g);
% xlabel('X方位向(m)'); ylabel('Y侧视距离向(m)'); zlabel('Z高度(m)');
% title('车载目标高程模型'); grid on;


% 
% point = [
%      0, 45, 0.1, 0, 0, 1;
%      0, 55, 0.1, 0, 0, 1;
%     -2, 50, 0.1, 0, 0, 1;
%      2, 50, 0.1, 0, 0, 1;
%     -5, 50, 0.1, 0, 0, 1;
%      5, 50, 0.1, 0, 0, 1;
%     -8, 50, 0.1, 0, 0, 1;
%      8, 50, 0.1, 0, 0, 1
% ];
