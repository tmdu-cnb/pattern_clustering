%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustering.m (イベントタイミング相関バージョン)
%
% 【概要】
%  1) presentation.csv を読み込み
%     - ROIごとの蛍光時系列 (フレームデータ)
%     - ROIの (x, y) 座標
%     - ROI番号 など
%
%  2) 「座標ベース」と「イベントタイミング相関ベース」で別々に階層的クラスタリングを実施
%     - 座標クラスタリングは従来通り
%     - 発火パターンクラスタリングは、閾値を超えたフレームをイベント(=1)とみなし、
%       2つのROI間のJaccard係数を相関指標とする
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

nCols = width(T);
timeSeriesData = T{:, 1:(nCols-4)};

roiNumber    = T{:, nCols-3};  % ROI_number
is_cell      = T{:, nCols-2};  % is_cell
y_coords     = T{:, nCols-1};  % y座標
x_coords     = T{:, nCols};    % x座標

N = size(timeSeriesData, 1);   % ROI数
numFrames = size(timeSeriesData, 2);
fprintf('Loaded %d ROIs, each with %d frames.\n', N, numFrames);

%% 2. 座標ベースのクラスタリング（従来通り）
roiXY = [x_coords, y_coords];
distCoord = pdist(roiXY, 'euclidean');
Z_coord = linkage(distCoord, 'average');

figure('Name','Dendrogram for Coordinate-based Clustering');
dendrogram(Z_coord, 0);
title('Coordinate-based Dendrogram');
xlabel('ROI index');
ylabel('Distance');

% しきい値を適宜設定 (例: 全距離の最大値の 30%に相当する位置)
cutoffCoord = 0.3 * max(Z_coord(:,3));
clusterCoordIdx = cluster(Z_coord, 'cutoff', cutoffCoord, 'criterion', 'distance');

% 座標クラスタの種類数（ユニーク数）
uniqueCoordClusters = unique(clusterCoordIdx);
numCoordClusters = length(uniqueCoordClusters);

% A,B,C... ラベル付け (最大26想定)
uniqueCoordClusters = unique(clusterCoordIdx);
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

%% 3. 発火パターン(イベントタイミング)ベースのクラスタリング

% 3-1) イベント抽出
%     例: 閾値を定義し、timeSeriesData(i,t) > threshold のフレームを 1(イベント) とする
threshold = 5; % 例: 要調整
binaryEvents = timeSeriesData > threshold;  % (N x T) の logical配列

% 3-2) ROIペア間のイベントタイミング相関（ここでは Jaccard係数）を計算
%      Jaccard = (A AND B) / (A OR B)
%      ただし A,B は {0,1}のbinary vector

corrMat_events = eye(N);  % (N x N)
for i = 1:N
    for j = i+1:N
        A = binaryEvents(i,:);
        B = binaryEvents(j,:);
        
        overlap = sum(A & B);       % 両方1のフレーム数
        unionAB = sum(A | B);       % どちらか1のフレーム数
        if unionAB == 0
            r = 0;  % 両方ほぼイベントがないなら相関0とする
        else
            r = overlap / unionAB;  % Jaccard index
        end
        
        corrMat_events(i,j) = r;
        corrMat_events(j,i) = r;
    end
end

% 3-3) 距離行列を作り、階層クラスタリング
distPatternMat = 1 - corrMat_events;      % 類似度が高いほど距離小
distPatternVec = squareform(distPatternMat);
Z_pattern = linkage(distPatternVec, 'ward');

figure('Name','Dendrogram for Event-timing Clustering');
dendrogram(Z_pattern, 0);
title('Event-timing-based Dendrogram (Jaccard)');
xlabel('ROI index');
ylabel('Distance');

% cutoff設定
cutoffPattern = 0.5;
clusterPatternIdx = cluster(Z_pattern, 'cutoff', cutoffPattern, ...
                            'criterion','distance');

%発火パターンクラスタの種類数
uniquePatternClusters = unique(clusterPatternIdx);
numPatternClusters = length(uniquePatternClusters);

patternClusterLabelCell = arrayfun(@(x) num2str(x), clusterPatternIdx, ...
                                   'UniformOutput', false);

%% 4. 座標クラスタ & イベントタイミングクラスタ → 複合ラベル
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
end







% 5-2) 画像上での可視化 (例)
figure('Name','Clustering result'); 
imshow('Drd1_1_AAV-L7-RFP_L7-j7s_P9.tif');  % 例: 背景画像
hold on;

uniqueCoordClusters2 = unique(clusterCoordIdx);
numC = length(uniqueCoordClusters2);
cmap = lines(numC);

for c = 1:numC
    idxC = find(clusterCoordIdx == c);
    scatter(x_coords(idxC), y_coords(idxC), 30, ...
        'MarkerFaceColor', cmap(c,:), 'MarkerEdgeColor','k', ...
        'DisplayName', ['CoordCluster ' lettersAll{c}]);
end
legend('show');

% ROIごとにラベル(A1, B2...)を文字表示
for i = 1:N
    text(x_coords(i), y_coords(i), finalClusterLabel{i}, ...
        'Color','w','FontSize',8,'HorizontalAlignment','center', ...
        'VerticalAlignment','middle','FontWeight','bold');
end

% 5-3) クラスターごとに枠線表示 (例: boundary)
uniqueFinalLabels = unique(finalClusterLabel);
nLabels = length(uniqueFinalLabels);
cmap2 = lines(nLabels);

for cIdx = 1:nLabels
    thisLabel = uniqueFinalLabels{cIdx};
    idxCluster = find(strcmp(finalClusterLabel, thisLabel));
    
    if length(idxCluster) < 3, continue; end
    
    x_c = x_coords(idxCluster);
    y_c = y_coords(idxCluster);
    points = [x_c, y_c];
    k = boundary(points(:,1), points(:,2), 0.8);
    plotColor = cmap2(cIdx,:);
    patch('XData', points(k,1), 'YData', points(k,2), ...
          'EdgeColor', plotColor, 'FaceColor','none', ...
          'LineWidth', 2, 'LineStyle','-');
      
    cx = mean(points(k,1));
    cy = mean(points(k,2));
    text(cx, cy, thisLabel, 'Color', plotColor, 'FontSize', 10, ...
        'FontWeight','bold','HorizontalAlignment','center',...
        'VerticalAlignment','middle');
end
hold off;

% 5-4) 結果のテーブル出力
Tout = table(roiNumber, x_coords, y_coords, ...
             clusterCoordIdx, clusterPatternIdx, finalClusterLabel, ...
             'VariableNames', ...
             {'ROI_number','x','y','CoordClusterID','PatternClusterID','FinalLabel'});
writetable(Tout, 'clustering_result.csv');
fprintf('Saved clustering_result.csv\n');
