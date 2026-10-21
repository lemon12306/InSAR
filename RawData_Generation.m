clear; close all; clc;
rng(20261011);
SAR_Parameters;

ThetaTa = deg2rad(55);
Tp = 30e-9;

%%%%%%%%%%%%%%%%%%%%%%%%%% 1. 生成原始回波
[Echo,~] = RawDataGen_rand( ...
    point,Rr,nan,B,fc,PRF,Ym,H,V,ThetaTa,Bz,0,Tp,'raw');

%%%%%%%%%%%%%%%%%%%%%%%%%% 2. 保存原始回波及参数
NoiseState = rng;% 保留添加接收噪声前的随机状态

save Raw_Echo.mat
