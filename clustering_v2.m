%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustering.m
% 
% 【概要】
%  1) presentation.csv を読み込み、
%     - ROIごとの蛍光時系列 (フレームデータ)
%     - ROIの (x, y) 座標
%     - ROI番号 など
%    を取得
%
%  2) 「座標ベース」と「発火パターン（相関）ベース」で別々に階層的クラスタリングを実施
%     - 座標クラスタリングは従来通り
%     - 発火パターンクラスタリングでは、両ROIの発火が同時に0のフレームを除外して相関を計算
%
%  3) それぞれの結果を組み合わせて、A1, A2, B1, B2, ... のようなラベルを付与
%
%  4) 画像への描画や結果の保存を行う
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1. データ読み込み
csvFile = 'presentation.csv';
T = readtable(csvFile);

% 下記は "Frame_1 ~ Frame_n" が最初～(end-4)列にある想定
nCols = width(T);
timeSeriesData = T{:, 1:(nCols-4)};

roiNumber    = T{:, nCols-3};  % ROI_number 列
is_cell      = T{:, nCols-2};  % is_cell 列
y_coords     = T{:, nCols-1};  % y座標
x_coords     = T{:, nCols};    % x座標

N = size(timeSeriesData, 1);   % ROI数
numFrames = size(timeSeriesData, 2);  % フレーム数

fprintf('Loaded %d ROIs, each with %d frames of data.\n', N, numFrames);

%% 2. 座標ベースのクラスタリング
roiXY = [x_coords, y_coords];  % (N x 2)

distCoord = pdist(roiXY, 'euclidean');
Z_coord = linkage(distCoord, 'average');

% デンドログラム可視化 (座標)
figure('Name','Dendrogram for Coordinate-based Clustering');
dendrogram(Z_coord, 0);
title('Coordinate-based Dendrogram');
xlabel('ROI index');
ylabel('Distance');

% クラスタリング (例: cutoff を 0.3 * max(Z_coord(:,3)) )
cutoffCoord = 0.3 * max(Z_coord(:,3));
clusterCoordIdx = cluster(Z_coord, 'cutoff', cutoffCoord, ...
                          'criterion', 'distance');

% 座標クラスタの種類数（ユニーク数）
uniqueCoordClusters = unique(clusterCoordIdx);
numCoordClusters = length(uniqueCoordClusters);

% A,B,C... ラベル付け (最大26想定)
uniqueCoordClusters = unique(clusterCoordIdx);
nCoordClusters = length(uniqueCoordClusters);
lettersAll = arrayfun(@(x) char('A' + x - 1), 1:26, 'UniformOutput', false);


coordClusterLabelCell = cell(N,1);
for i = 1:N
    cidx = clusterCoordIdx(i);
    if cidx <= 26
        coordClusterLabelCell{i} = lettersAll{cidx};
    else
        coordClusterLabelCell{i} = ['X' num2str(cidx)]; % 暫定
    end
end

%% 3. 発火パターンベースのクラスタリング
% 両方のROIが同時に0 だったフレームは相関計算から除外する

% N x T
corrMat_custom = eye(N);  % 対角(自己相関)は 1 or 0 にしてもOK
for i = 1:N
    xi = timeSeriesData(i, :);
    for j = i+1 : N
        xj = timeSeriesData(j, :);
        
        % 両方が0のフレームを除外
        mask = ~((xi == 0) & (xj == 0));  % true の箇所が「使う」フレーム
        xi2 = xi(mask);
        xj2 = xj(mask);
        
        if length(xi2) < 2
            % 相関定義できないので r=0 とする (お好みで NaN や1でも可)
            r = 0;
        else
            % ここでピアソン相関を計算
            tmp = corrcoef(xi2, xj2);
            r = tmp(1,2);
        end
        
        corrMat_custom(i,j) = r;
        corrMat_custom(j,i) = r;
    end
end

% 距離行列に変換
distPatternMat = 1 - corrMat_custom;  
distPatternVec = squareform(distPatternMat);

% linkage  (average, ward, complete など試せる)
Z_pattern = linkage(distPatternVec, 'ward');

% デンドログラム可視化 (発火パターン)
figure('Name','Dendrogram for Pattern-based Clustering');
dendrogram(Z_pattern, 0);
title('Pattern-based Dendrogram (exclude frames both=0)');
xlabel('ROI index');
ylabel('Distance');

% クラスタリング (cutoffPattern=0.5 等は要調整)
cutoffPattern = 0.5;
clusterPatternIdx = cluster(Z_pattern, 'cutoff', cutoffPattern, ...
                            'criterion', 'distance');

%発火パターンクラスタの種類数
uniquePatternClusters = unique(clusterPatternIdx);
numPatternClusters = length(uniquePatternClusters);


% 数字ラベル (1, 2, 3, ...) → 文字列に
patternClusterLabelCell = arrayfun(@(x) num2str(x), clusterPatternIdx, ...
                                   'UniformOutput', false);

%% 4. 座標クラスタ & 発火パターンクラスタ → 複合ラベル付け
finalClusterLabel = cell(N,1);
for i = 1:N
    finalClusterLabel{i} = [coordClusterLabelCell{i}, patternClusterLabelCell{i}];
end

% (★変更点1) 最終ラベルの種類数を定義
uniqueFinalLabels = unique(finalClusterLabel);
numFinalClusters = length(uniqueFinalLabels);  % ← 追加した行

%% 5. 結果の可視化・出力

% 5-1) ターミナル表示
disp('=== Final cluster labels (ROI-based) ===');
for i = 1:N
    fprintf('ROI_number=%d -> %s\n', roiNumber(i), finalClusterLabel{i});
end

%--------------------------------------------------------------------------
% 5-2) 各クラスタの度数表示
%--------------------------------------------------------------------------

% (A) それぞれの「クラスター数」をターミナル表示
fprintf('\n=== Cluster Count Info ===\n');
fprintf('座標クラスタの数: %d\n', numCoordClusters);
fprintf('発火パターンクラスタの数: %d\n', numPatternClusters);
fprintf('最終ラベル(組み合わせ)の数: %d\n', numFinalClusters);

% (A) 座標クラスタごとの度数
disp('=== Frequency of each Coordinate Cluster ===');
uCoord = unique(clusterCoordIdx);
for c = 1:length(uCoord)
    cID = uCoord(c);  % クラスタ番号
    freqCount = sum(clusterCoordIdx == cID);
    % ラベル(A,B,...)も分かるように対応付け
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

% (C) 最終ラベル(例: A1, B2...)ごとの度数
disp('=== Frequency of each Final Cluster ===');
uFinal = unique(finalClusterLabel);
for c = 1:length(uFinal)
    clusterName = uFinal{c};
    freqCount = sum(strcmp(finalClusterLabel, clusterName));
    fprintf('  FinalLabel %s : %d ROIs\n', clusterName, freqCount);
    fprintf('  最終ラベル %s → %d個\n', clusterName, freqCount);
end






% 5-2) 画像上で可視化 (簡易例)
figure('Name','Clustering result'); 
imshow('Drd1_1_AAV-L7-RFP_L7-j7s_P9.tif');  % 背景画像 (例)
hold on;

% 座標クラスタ (A,B,C...) 単位で散布図
nCoordClusters = length(uniqueCoordClusters);
cmap = lines(nCoordClusters);
for c = 1:nCoordClusters
    idxC = find(clusterCoordIdx == c);
    scatter(x_coords(idxC), y_coords(idxC), 30, ...
        'MarkerFaceColor', cmap(c,:), 'MarkerEdgeColor', 'k', ...
        'DisplayName',['CoordCluster ' lettersAll{c}]);
end
legend('show');

% ROIごとにラベル(A1, B2...)を文字表示
for i = 1:N
   text(x_coords(i), y_coords(i), finalClusterLabel{i}, ...
       'Color','w','FontSize',8,'HorizontalAlignment','center', ...
       'VerticalAlignment','middle','FontWeight','bold');
end

% 5-3) クラスターごとに輪郭(例: boundary)を引く
uniqueFinalLabels = unique(finalClusterLabel);
nLabels = length(uniqueFinalLabels);
cmap2 = lines(nLabels);

for cIdx = 1:nLabels
    thisLabel = uniqueFinalLabels{cIdx};
    idxCluster = find(strcmp(finalClusterLabel, thisLabel));
    
    if length(idxCluster) < 3
        continue;  % ROIが2点以下ならスキップ
    end
    
    x_c = x_coords(idxCluster);
    y_c = y_coords(idxCluster);
    points = [x_c, y_c];
    
    k = boundary(points(:,1), points(:,2), 0.8);
    plotColor = cmap2(cIdx,:);
    patch('XData', points(k,1), 'YData', points(k,2), ...
          'EdgeColor', plotColor, 'FaceColor', 'none', ...
          'LineWidth', 2, 'LineStyle','-');
      
    % 領域の重心付近にラベル表示
    cx = mean(points(k,1));
    cy = mean(points(k,2));
    text(cx, cy, thisLabel, ...
        'Color', plotColor, 'FontSize', 10, 'FontWeight','bold', ...
        'HorizontalAlignment','center','VerticalAlignment','middle');
end
hold off;

% 5-4) 結果をテーブルで出力
Tout = table(roiNumber, x_coords, y_coords, ...
             clusterCoordIdx, clusterPatternIdx, finalClusterLabel, ...
             'VariableNames',{'ROI_number','x','y','CoordClusterID','PatternClusterID','FinalLabel'});
writetable(Tout, 'clustering_result.csv');
fprintf('Saved clustering_result.csv\n');
