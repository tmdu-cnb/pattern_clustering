%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustering.m
% 
% 【概要】
%   1) presentation.csv を読み込み、
%      - ROIごとの蛍光時系列データ (Frame)
%      - ROIのx,y座標
%      - ROI番号 など
%   2) 「座標ベース」と「発火パターンベース」で別々に階層的クラスタリング
%   3) それぞれの結果を組み合わせて A1, B2... のような複合ラベルを付与
%   4) 画像への描画やCSV出力
%
% 【変更内容】
%  ・最終ラベルに属するROIが1つしか無いクラスター（単独クラスター）を画像表示から除外
%  ・単独クラスターがいくつ存在したかを表示
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1. データ読み込み
csvFile = 'presentation.csv';
T = readtable(csvFile);

nCols = width(T);
timeSeriesData = T{:, 1:(nCols-4)};
roiNumber    = T{:, nCols-3};
is_cell      = T{:, nCols-2};
y_coords     = T{:, nCols-1};
x_coords     = T{:, nCols};

N = size(timeSeriesData, 1);
numFrames = size(timeSeriesData, 2);
fprintf('Loaded %d ROIs, each with %d frames of data.\n', N, numFrames);

%% 2. 座標ベースクラスタリング
roiXY = [x_coords, y_coords];
distCoord = pdist(roiXY, 'euclidean');
Z_coord = linkage(distCoord, 'average');

figure('Name','Dendrogram for Coordinate-based Clustering');
dendrogram(Z_coord, 0);
title('Coordinate-based Dendrogram');
xlabel('ROI index');
ylabel('Distance');

cutoffCoord = 0.3 * max(Z_coord(:,3));
clusterCoordIdx = cluster(Z_coord, 'cutoff', cutoffCoord, 'criterion', 'distance');

uniqueCoordClusters = unique(clusterCoordIdx);
numCoordClusters = length(uniqueCoordClusters);

lettersAll = arrayfun(@(x) char('A' + x - 1), 1:26, 'UniformOutput', false);

coordClusterLabelCell = cell(N,1);
for i = 1:N
    cidx = clusterCoordIdx(i);
    if cidx <= 26
        coordClusterLabelCell{i} = lettersAll{cidx};
    else
        coordClusterLabelCell{i} = ['X' num2str(cidx)];
    end
end

%% 3. 発火パターンベースクラスタリング
corrMat_custom = eye(N);
for i = 1:N
    xi = timeSeriesData(i, :);
    for j = i+1:N
        xj = timeSeriesData(j, :);
        mask = ~((xi == 0) & (xj == 0));  
        xi2 = xi(mask);
        xj2 = xj(mask);
        if length(xi2) < 2
            r = 0;
        else
            tmp = corrcoef(xi2, xj2);
            r = tmp(1,2);
        end
        corrMat_custom(i,j) = r;
        corrMat_custom(j,i) = r;
    end
end

distPatternMat = 1 - corrMat_custom;
distPatternVec = squareform(distPatternMat);
Z_pattern = linkage(distPatternVec, 'ave');

figure('Name','Dendrogram for Pattern-based Clustering');
dendrogram(Z_pattern, 0);
title('Pattern-based Dendrogram (exclude frames both=0)');
xlabel('ROI index');
ylabel('Distance');

cutoffPattern = 0.5;
clusterPatternIdx = cluster(Z_pattern, 'cutoff', cutoffPattern, 'criterion','distance');

uniquePatternClusters = unique(clusterPatternIdx);
numPatternClusters = length(uniquePatternClusters);

patternClusterLabelCell = arrayfun(@(x) num2str(x), clusterPatternIdx, ...
    'UniformOutput', false);

%% 4. 複合ラベル付け (A1, B2, ...)
finalClusterLabel = cell(N,1);
for i = 1:N
    finalClusterLabel{i} = [coordClusterLabelCell{i}, patternClusterLabelCell{i}];
end

% ★変更点1: 最終ラベルの種類数
uniqueFinalLabels = unique(finalClusterLabel);
numFinalClusters = length(uniqueFinalLabels);

%% 5. 結果の可視化・出力

% 5-1) ターミナル表示
disp('=== Final cluster labels (ROI-based) ===');
for i = 1:N
    fprintf('ROI_number=%d -> %s\n', roiNumber(i), finalClusterLabel{i});
end

%--------------------------------------------------------------------------
% 5-2) クラスター数や度数を表示
%--------------------------------------------------------------------------
fprintf('\n=== Cluster Count Info ===\n');
fprintf('座標クラスタの数: %d\n', numCoordClusters);
fprintf('発火パターンクラスタの数: %d\n', numPatternClusters);
fprintf('最終ラベル(組み合わせ)の数: %d\n', numFinalClusters);

% (A) 座標クラスタごとの度数
disp('=== Frequency of each Coordinate Cluster ===');
uCoord = unique(clusterCoordIdx);
for c = 1:length(uCoord)
    cID = uCoord(c);
    freqCount = sum(clusterCoordIdx == cID);
    if cID <= 26
        labelStr = lettersAll{cID};
    else
        labelStr = ['X' num2str(cID)];
    end
    fprintf('  CoordCluster %s (ID=%d): %d ROIs\n', labelStr, cID, freqCount);
end

% (B) 発火パターンクラスタごとの度数
disp('=== Frequency of each Pattern Cluster ===');
uPattern = unique(clusterPatternIdx);
for c = 1:length(uPattern)
    cID = uPattern(c);
    freqCount = sum(clusterPatternIdx == cID);
    fprintf('  PatternCluster %d : %d ROIs\n', cID, freqCount);
end

% (C) 最終ラベル(A1, B2...)ごとの度数
disp('=== Frequency of each Final Cluster ===');
uFinal = unique(finalClusterLabel);
freqFinal = zeros(size(uFinal));

for c = 1:length(uFinal)
    clusterName = uFinal{c};
    freqCount = sum(strcmp(finalClusterLabel, clusterName));
    freqFinal(c) = freqCount;
    fprintf('  FinalLabel %s : %d ROIs\n', clusterName, freqCount);
end

%--------------------------------------------------------------------------
% ★変更点2: 単独クラスター(ROIが1つだけ)の除外ロジック
%--------------------------------------------------------------------------
% freqFinal == 1 のラベルを "単独クラスター" とみなす
singletonMask = (freqFinal == 1);
singletonLabels = uFinal(singletonMask);

% "単独クラスター" がいくつあったか
numSingletonClusters = sum(singletonMask);
fprintf('\n 単独クラスター (ROIが1つしかない複合ラベル) の数: %d\n', numSingletonClusters);

% 画像表示では、単独クラスターに属するROIは除外
% → keepMask = true if ROI belongs to a multi-ROI cluster
keepMask = true(N,1);
for i = 1:N
   if ismember(finalClusterLabel{i}, singletonLabels)
       keepMask(i) = false;
   end
end

%--------------------------------------------------------------------------
% 5-3) 画像上で可視化 (単独クラスターを除いたROIのみプロット)
%--------------------------------------------------------------------------
figure('Name','Clustering result (exclude singleton clusters)');
imshow('Drd1_1_AAV-L7-RFP_L7-j7s_P9.tif');  
hold on;

% 散布図 (座標クラスタは変わらず)
%   → keepMask でフィルタ
keptCoordIdx = clusterCoordIdx(keepMask);
keptX = x_coords(keepMask);
keptY = y_coords(keepMask);
keptFinalLabel = finalClusterLabel(keepMask);

% 既存の座標クラスタ数は numCoordClusters のままだが，
% 単独クラスター削除で実際に何個残るかは再計算してもいい．
% ここではそのまま lines(nCoordClusters) で着色
nCoordClusters = length(uniqueCoordClusters);
cmap = lines(nCoordClusters);

for c = 1:nCoordClusters
    idxC = find(keptCoordIdx == c);
    if isempty(idxC), continue; end
    scatter(keptX(idxC), keptY(idxC), 30, ...
        'MarkerFaceColor', cmap(c,:), 'MarkerEdgeColor','k', ...
        'DisplayName',['CoordCluster ' lettersAll{c}]);
end
legend('show');

% ROIごとの文字ラベル
for i = 1:sum(keepMask)
   text(keptX(i), keptY(i), keptFinalLabel{i}, ...
       'Color','w','FontSize',8,'HorizontalAlignment','center', ...
       'VerticalAlignment','middle','FontWeight','bold');
end

% 複合ラベルごとに輪郭
% 輪郭描画するのは keepMask後のROIのみ
uFinalKept = unique(keptFinalLabel);
cmap2 = lines(length(uFinalKept));

for cIdx = 1:length(uFinalKept)
    thisLabel = uFinalKept{cIdx};
    idxCluster = find(strcmp(keptFinalLabel, thisLabel));
    if length(idxCluster) < 3, continue; end
    
    x_c = keptX(idxCluster);
    y_c = keptY(idxCluster);
    points = [x_c, y_c];
    
    k = boundary(points(:,1), points(:,2), 0.8);
    plotColor = cmap2(cIdx,:);
    patch('XData', points(k,1), 'YData', points(k,2), ...
          'EdgeColor', plotColor, 'FaceColor','none', ...
          'LineWidth', 2, 'LineStyle','-');
      
    cx = mean(points(k,1));
    cy = mean(points(k,2));
    text(cx, cy, thisLabel, 'Color', plotColor, 'FontSize', 10, ...
         'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
end
hold off;

%--------------------------------------------------------------------------
% 5-4) テーブル出力
%  単独クラスターのROIはCSV上から削除するかどうかはお好み
%  → 今回は削除せず、すべて出力(要件に合わせて変更可)
%--------------------------------------------------------------------------
Tout = table(roiNumber, x_coords, y_coords, ...
    clusterCoordIdx, clusterPatternIdx, finalClusterLabel, ...
    'VariableNames',{'ROI_number','x','y','CoordClusterID','PatternClusterID','FinalLabel'});
writetable(Tout, 'clustering_result.csv');
fprintf('Saved clustering_result.csv\n');
