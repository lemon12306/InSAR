function [Data, simInfo] = RawDataGen_rand(point, Rr, nan, B, fc, PRF, Ym, H, V, ThetaTa, Bz, Bzy, Tp, Mode)
% 输入：
%   point(:,1) — x 坐标，m
%   point(:,2) — y 坐标，m
%   point(:,3) — z 高度，m
%   point(:,4) — x 向速度（静止目标为 0）
%   point(:,5) — y 向速度（静止目标为 0）
%   point(:,6) — 复散射系数
%   Rr         — 等效斜距轴，m
%   nan        — 方位脉冲数
%   B          — 距离向带宽，Hz
%   fc         — 载频，Hz
%   PRF        — 脉冲重复频率，Hz
%   Ym         — 两接收天线 y 坐标 [y1, y2]
%   H          — 接收天线高度，m
%   V          — 平台运动速度，m/s
%   ThetaTa    — 方位向波束宽度，rad
%   Bz         — z 向基线（两接收天线高度差），m（可选，默认 0）
%   Bzy        — 发射天线相对于接收天线 H 的高度偏移量，m（可选，默认 0）
%   Tp         — 发射脉冲宽度，s（Mode='raw' 时使用，可选，默认 20 ns）
%   Mode       — 'rc'  输出距离压缩 HRRP（默认，兼容原 BP 流程）
%                'raw' 输出原始 LFM 回波（脉压前，供 RDA 做距离脉压）
%
% 输出：
%   Data(:,:,1) — 通道 A 数据
%   Data(:,:,2) — 通道 B 数据
%   simInfo     — 仿真几何信息

if nargin < 11 || isempty(Bz)
    Bz = 0;
end
if nargin < 12 || isempty(Bzy)
    Bzy = 0;
end
if nargin < 14 || isempty(Mode)
    Mode = 'rc';
end
rawMode = strcmpi(Mode, 'raw');
if rawMode && (nargin < 13 || isempty(Tp))
    Tp = 20e-9;   % 默认脉宽
end

C = 3e8;
lambda = C/fc;
nrange = length(Rr);
nchannel = 2;
ta = (-(nan-1)/2:(nan-1)/2)/PRF;

if rawMode
    Kr = B / Tp;                 % 距离向调频率
    dR = Rr(2) - Rr(1);
    dTau = 2*dR / C;
    Rref = Rr(round(nrange/2));  % 快时间参考斜距
    tau_r = 2*(Rr - Rref) / C;   % 快时间轴，nrange×1
end

% 天线初始位置
Tx0  = [0, 0, H];
Rx10 = [0, Ym(1), H + Bz];
Rx20 = [0, Ym(2), H + 2*Bz];

% 散射点参数
pointn = size(point, 1);
if size(point, 2) >= 6
    amp = point(:, 6);
else
    amp = ones(pointn, 1);
end

Data = complex(zeros(nrange, nan, nchannel));

xp0 = real(point(:,1));
yp0 = real(point(:,2));
zp0 = real(point(:,3));
vx  = real(point(:,4));
vy  = real(point(:,5));

for np = 1:nan
    t = ta(np);

    % 目标位置（静止目标，vx=vy=0，位置不随脉冲变化）
    xp = xp0 + vx*t;
    yp = yp0 + vy*t;
    zp = zp0;

    % 当前脉冲时刻的天线位置
    tx  = Tx0  + [V*t, 0, 0];
    rx1 = Rx10 + [V*t, 0, 0];
    rx2 = Rx20 + [V*t, 0, 0];

    % 侧视波束，只保留波束方位角范围内的散射点
    look_angle = atan2(xp - tx(1), yp - tx(2));
    beam_win = abs(look_angle) <= ThetaTa/2;

    % 发射和接收斜距
    Rtx  = sqrt((xp - tx(1)).^2  + (yp - tx(2)).^2  + (zp - tx(3)).^2);
    Rrx1 = sqrt((xp - rx1(1)).^2 + (yp - rx1(2)).^2 + (zp - rx1(3)).^2);
    Rrx2 = sqrt((xp - rx2(1)).^2 + (yp - rx2(2)).^2 + (zp - rx2(3)).^2);

    % 单发双收等效斜距
    Req1 = 0.5*(Rtx + Rrx1);
    Req2 = 0.5*(Rtx + Rrx2);

    % 载频相位（视频相位）
    sig1 = amp(:) .* beam_win(:) .* exp(-1j*4*pi/lambda * Req1(:));
    sig2 = amp(:) .* beam_win(:) .* exp(-1j*4*pi/lambda * Req2(:));

    if rawMode
        % ---- 原始 LFM 回波（脉压前）----
        % 相对参考斜距的回波延迟
        del1 = 2*(Req1(:) - Rref) / C;
        del2 = 2*(Req2(:) - Rref) / C;
        % 快时间网格上的发射chirp，按回波延迟平移后叠加
        D1 = tau_r - del1.';
        D2 = tau_r - del2.';
        env1 = abs(D1) <= Tp/2;     % 脉内窗
        env2 = abs(D2) <= Tp/2;
        chirp1 = exp(1j*pi*Kr * D1.^2);
        chirp2 = exp(1j*pi*Kr * D2.^2);
        Data(:,np,1) = (chirp1 .* env1) * sig1;   % 通道 A
        Data(:,np,2) = (chirp2 .* env2) * sig2;   % 通道 B
    else
        % ---- 距离压缩 HRRP（理想 sinc 点扩散函数叠加）----
        S1 = sinc_local((2*B/C) * (Rr(:) - Req1(:).'));
        S2 = sinc_local((2*B/C) * (Rr(:) - Req2(:).'));
        Data(:,np,1) = S1 * sig1;
        Data(:,np,2) = S2 * sig2;
    end

    if mod(np, max(1, round(nan/5))) == 0
        fprintf('回波生成进度：%4.1f %%\n', 100*np/nan);
    end
end

% 保存仿真信息
if rawMode
    simInfo.model = 'raw LFM echo (pre range-compression)';
    simInfo.Kr = Kr;
    simInfo.Tp = Tp;
else
    simInfo.model = 'ideal range-compressed HRRP';
end
simInfo.rangeDefinition = 'R_eq = 0.5*(R_tx + R_rx)';
simInfo.Rr = Rr;
simInfo.ta = ta;
simInfo.Tx0 = Tx0;
simInfo.Rx10 = Rx10;
simInfo.Rx20 = Rx20;
simInfo.lambda = lambda;
simInfo.B = B;
end

function y = sinc_local(x)
% Hann 窗加权 sinc，第一旁瓣 -32 dB（原 sinc 为 -13 dB）
% 主瓣展宽 ~1.44×，但可有效抑制距离旁瓣条纹
y = sinc_base(x);
y = 0.5*y + 0.25*(sinc_base(x+1) + sinc_base(x-1));
end

function y = sinc_base(x)
y = ones(size(x));
idx = abs(x) > 1e-12;
y(idx) = sin(pi*x(idx)) ./ (pi*x(idx));
end
