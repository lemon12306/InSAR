% 读取已完成的实验，在交互式 MATLAB 中保留图窗，不重新生成回波。
dataDir = fileparts(mfilename('fullpath'));
S = load(fullfile(dataDir,'随机幅度_无噪声_高程结果.mat'));
[~,centerCol] = min(abs(S.Ax_r));
y = S.Y_ground(:);
region = y>=42 & y<=58;
yr = y(region);
t = yr-mean(yr);
nRuns = numel(S.amplitudeSeeds);
analysisRows = zeros(nRuns,8);
fitCurves = cell(nRuns,1);
localCycles = cell(nRuns,1);

for k = 1:nRuns
    e = S.profileError(region,k);
    trendDesign = [ones(size(t)),t];
    trend = trendDesign*(trendDesign\e);
    baselineSSE = sum((e-trend).^2);
    bestSSE = Inf;
    % 在与前次频谱分析相同的频率范围内搜索整段最优正弦。
    for f = 0.2:0.001:2.5
        design = [trendDesign,sin(2*pi*f*t),cos(2*pi*f*t)];
        fitted = design*(design\e);
        sse = sum((e-fitted).^2);
        if sse<bestSSE
            bestSSE = sse;
            bestF = f;
            fitCurves{k} = fitted;
        end
    end
    sineR2 = 1-bestSSE/baselineSSE;

    % 用去线性趋势后的相邻向上过零点定义候选局部周期。
    % 只保留至少8个原始采样点的区间，避免极短片段过拟合。
    residual = e-trend;
    crossings = find(residual(1:end-1)<=0 & residual(2:end)>0);
    crossingY = yr(crossings)-residual(crossings).* ...
        (yr(crossings+1)-yr(crossings))./ ...
        (residual(crossings+1)-residual(crossings));
    cycles = zeros(0,4);
    for j = 1:numel(crossingY)-1
        a = crossingY(j);
        b = crossingY(j+1);
        use = yr>=a & yr<=b;
        if nnz(use)<8
            continue;
        end
        q = 2*pi*(yr(use)-a)/(b-a);
        design = [ones(nnz(use),1),sin(q),cos(q)];
        segment = residual(use);
        fitted = design*(design\segment);
        sst = sum((segment-mean(segment)).^2);
        r2 = 1-sum((segment-fitted).^2)/max(sst,eps);
        cycles(end+1,:) = [a,b,b-a,r2]; %#ok<SAGROW>
    end
    localCycles{k} = cycles;
    if isempty(cycles)
        cycleStats = [NaN NaN NaN NaN];
    else
        cycleStats = [min(cycles(:,3)),max(cycles(:,3)), ...
            median(cycles(:,4)),mean(cycles(:,4)>=0.9)];
    end
    analysisRows(k,:) = [S.amplitudeSeeds(k),1/bestF,sineR2, ...
        size(cycles,1),cycleStats];
end
SineAnalysis = array2table(analysisRows,'VariableNames', ...
    {'Seed','BestGlobalWavelength_m','GlobalSineR2','LocalCycleCount', ...
    'MinLocalPeriod_m','MaxLocalPeriod_m','MedianLocalSineR2', ...
    'FractionLocalR2Above09'});
disp(SineAnalysis);
save(fullfile(dataDir,'onebit_sine_analysis.mat'), ...
    'SineAnalysis','localCycles','yr','fitCurves');

% 三维图直接使用反演矩阵，中心黑线与剖面图取自同一列。
f3d = figure('Name','单比特三维高程：可拖动旋转，黑线为中心剖面', ...
    'NumberTitle','off','Color','w','Position',[80 100 1400 650]);
tl = tiledlayout(f3d,1,nRuns,'TileSpacing','compact');
allHeights = cell2mat(cellfun(@(h) h(:),S.hBit,'UniformOutput',false));
heightLimits = [min(allHeights(:)),max(allHeights(:))];
for k = 1:nRuns
    ax = nexttile(tl);
    showSurface(ax,S.Ax_r,y,S.hBit{k},heightLimits);
    hold(ax,'on');
    plot3(ax,repmat(S.Ax_r(centerCol),size(y)),y, ...
        S.hBit{k}(:,centerCol),'k-','LineWidth',1.7);
    title(ax,sprintf('seed=%d',S.amplitudeSeeds(k)));
end
title(tl,'无噪声零阈值单比特高程（统一坐标和色标，含全部像素）');
rotate3d(f3d,'on');
savefig(f3d,fullfile(dataDir,'onebit_3d_interactive.fig'));

fcmp = figure('Name','真实、全精度和单比特：统一尺度三维比较', ...
    'NumberTitle','off','Color','w','Position',[100 130 1400 650]);
tl = tiledlayout(fcmp,1,3,'TileSpacing','compact');
maps = {S.h_true,S.hFull{1},S.hBit{1}};
names = {'真实高程','全精度反演','零阈值单比特反演'};
for k = 1:3
    ax = nexttile(tl);
    showSurface(ax,S.Ax_r,y,maps{k},heightLimits);
    title(ax,names{k});
end
title(tl,'相同视角、纵轴和色标；seed=1101');
rotate3d(fcmp,'on');

fdiag = figure('Name','中心剖面与正弦拟合：检查每段波形', ...
    'NumberTitle','off','Color','w','Position',[120 70 1250 850]);
tl = tiledlayout(fdiag,nRuns,2,'TileSpacing','compact');
for k = 1:nRuns
    ax = nexttile(tl);
    plot(ax,y,S.h_true(:,centerCol),'k-', ...
        y,S.hFull{k}(:,centerCol),'-',y,S.hBit{k}(:,centerCol),'--');
    grid(ax,'on');
    xlabel(ax,'距离向/m'); ylabel(ax,'高程/m');
    title(ax,sprintf('中心剖面，seed=%d',S.amplitudeSeeds(k)));
    legend(ax,'真实','全精度','单比特','Location','best');
    ax = nexttile(tl);
    plot(ax,yr,S.profileError(region,k),'-',yr,fitCurves{k},'--', ...
        'LineWidth',1.1);
    grid(ax,'on');
    xlabel(ax,'距离向/m'); ylabel(ax,'单比特高程误差/m');
    title(ax,sprintf('整段正弦拟合 R^2=%.3f',analysisRows(k,3)));
    legend(ax,'实际误差','线性趋势+最优单正弦','Location','best');
end
savefig(fdiag,fullfile(dataDir,'onebit_profile_diagnostics.fig'));
exportgraphics(fdiag,fullfile(dataDir,'onebit_profile_diagnostics.png'), ...
    'Resolution',180);
drawnow;
viewerReady = true;
save(fullfile(dataDir,'onebit_viewer_ready.mat'),'viewerReady','SineAnalysis');
figure(f3d);
fprintf('Interactive figures are ready. MATLAB and figures remain open.\n');

function showSurface(ax,x,y,h,limits)
surf(ax,x,y,h,'EdgeColor','none');
view(ax,45,35);
axis(ax,'tight');
zlim(ax,limits); clim(ax,limits);
colormap(ax,parula); colorbar(ax);
xlabel(ax,'方位向/m'); ylabel(ax,'距离向/m'); zlabel(ax,'高程/m');
grid(ax,'on');
end
