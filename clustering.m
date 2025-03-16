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
%
%  3) それぞれの結果を組み合わせて、A1, A2, B1, B2, ... のようなラベルを付与
%
%  4) 画像への描画や結果の保存を行う（必要に応じて stat も参照）
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1. データ読み込み
% presentation.csv の想定フォーマット（列の並び一例）
% [Frame_1, Frame_2, ... , Frame_n, ROI_number, is_cell, y, x]
%  → 具体的には detectROI.m の最後で出力した形に合わせる

csvFile = 'presentation.csv';
T = readtable(csvFile);

% ヘッダーや構成はユーザーの環境によって異なる場合があるので、
% 必要に応じて「何列目が何に当たるか」を確認しながら修正してください。

% --- ROIの時系列部分を取り出す（例: Frame_1 ~ Frame_n 列） ---
% ここでは仮に Frame_1 ~ Frame_n が最初～(end-4)列 に入っている想定とする
nCols = width(T);  % テーブル列数
timeSeriesData = T{:, 1:(nCols-4)};  % numeric配列として抜き出す

% --- ROI番号, 座標などを取り出す ---
roiNumber    = T{:, nCols-3};  % ROI_number 列 (最終4列前)
is_cell      = T{:, nCols-2};  % is_cell 列    (最終3列前)
y_coords     = T{:, nCols-1};  % y座標         (最終2列前)
x_coords     = T{:, nCols};    % x座標         (最終列)

% ROI の数やフレーム数
N = size(timeSeriesData, 1);  % ROI数
numFrames = size(timeSeriesData, 2);  % フレーム数

% 必要なら確認
fprintf('Loaded %d ROIs, each with %d frames of data.\n', N, numFrames);

%% 2. 座標ベースのクラスタリング
% (x, y) を使い、ユークリッド距離で階層的クラスタリング。
roiXY = [x_coords, y_coords];  % Nx2

% 2-1) 距離ベクトルを計算
distCoord = pdist(roiXY, 'euclidean');  % NxN の上三角成分をベクトル化したもの

% 2-2) linkage (階層的クラスタリング)
Z_coord = linkage(distCoord, 'average');  % 手法は "average", "ward", "complete" 等お好みで

% --- (追加) デンドログラムの可視化 ---
figure('Name','Dendrogram for Coordinate-based Clustering');
dendrogram(Z_coord, 0);
title('Coordinate-based Dendrogram');
xlabel('ROI index');
ylabel('Distance');


% 2-3) クラスタリングのしきい値またはクラスター数を決める
%   - しきい値の例 (cutoff)
cutoffCoord = 0.3 * max(Z_coord(:,3));  % とりあえず一例(割合で指定)
%   - 或いはクラスター数 k で指定も可能(例: 5クラスターに固定)
%     clusterCoordIdx = cluster(Z_coord, 'maxclust', 5);

clusterCoordIdx = cluster(Z_coord, 'cutoff', cutoffCoord, ...
                          'criterion', 'distance');

% clusterCoordIdx は [N x 1] のベクトルで，座標クラスタ番号(1,2,3,...)が入る

% 文字ラベル(A,B,C,...)を振りたい(最大26グループ想定)
uniqueCoordClusters = unique(clusterCoordIdx);
nCoordClusters = length(uniqueCoordClusters);
% 上限チェック (ここでは26以上になったらZから始める等、適宜対応)
lettersAll = arrayfun(@(x) char('A' + x - 1), 1:26, 'UniformOutput', false);

coordClusterLabelCell = cell(N,1);
for i = 1:N
    cidx = clusterCoordIdx(i);
    % cidx番目のアルファベットを割り当て (範囲外なら要工夫)
    if cidx <= 26
        coordClusterLabelCell{i} = lettersAll{cidx};
    else
        coordClusterLabelCell{i} = ['X' num2str(cidx)]; % 暫定対処
    end
end

%% 3. 発火パターンベースのクラスタリング
% ROIごとの発火パターン(時系列)間の相関を用いてクラスタリング。
% 相関係数 → 1 - 相関 を距離とする方法が一般的。

% 3-1) 相関行列
% 注意: MATLABの corr() は「列同士の相関」を計算するため、ROIが行なら転置が必要
corrMat = corr(timeSeriesData');  % (N x N)

% 3-2) 距離行列に変換 & pdist形式に
distPatternMat = 1 - corrMat;  % 相関が高いほど距離が小さい
distPatternVec = squareform(distPatternMat);  % pdist形式ベクトルへ

% 3-3) linkage
Z_pattern = linkage(distPatternVec, 'average');

% --- (追加) デンドログラムの可視化 ---
figure('Name','Dendrogram for Pattern-based Clustering');
dendrogram(Z_pattern, 0);  % 0指定で全ての枝を表示
title('Pattern-based Dendrogram');
xlabel('ROI index');
ylabel('Distance');

% 3-4) クラスタリング
cutoffPattern = 0.5;  % ここも適宜
clusterPatternIdx = cluster(Z_pattern, 'cutoff', cutoffPattern, ...
                            'criterion', 'distance');

% 数字ラベル(1, 2, 3, ...)をそのまま使うか、文字列化する
patternClusterLabelCell = arrayfun(@(x) num2str(x), clusterPatternIdx, ...
                                   'UniformOutput', false);

%% 4. 座標クラスタ(A,B,C...) と 発火パターンクラスタ(1,2,3,...) を結合し A1, B2 ... を作る
finalClusterLabel = cell(N,1);
for i = 1:N
    finalClusterLabel{i} = [coordClusterLabelCell{i}, patternClusterLabelCell{i}];
end

%% 5. 結果の可視化・出力
% 例1) ターミナルにクラスタラベルを表示 (ROIごと)
disp('=== Final cluster labels (ROI-based) ===');
for i = 1:N
    fprintf('ROI_number=%d -> %s\n', roiNumber(i), finalClusterLabel{i});
end

% 例2) 画像上にラベルを重ねて表示
%     すでに ROI_in_image.m がある場合、それを応用してもOKです。
%     ここでは簡易的に scatter + text で座標上にラベルを表示。
figure('Name','Clustering result'); 
imshow('Drd1_1_AAV-L7-RFP_L7-j7s_P9.tif'); % 背景に用いる画像
hold on;

% 色分け可視化したい場合
% 各複合クラスターに一意の色を付けるのは大変なので、まずは座標クラスタ (A,B,C...) 単位で色分けする例
nCoordClusters = length(uniqueCoordClusters);
cmap = lines(nCoordClusters);  % カラーマップ

for c = 1:nCoordClusters
    idxC = find(clusterCoordIdx == c);
    scatter(x_coords(idxC), y_coords(idxC), 30, ...
        'MarkerFaceColor', cmap(c,:), 'MarkerEdgeColor', 'k', ...
        'DisplayName',['CoordCluster ' lettersAll{c}]);
end
legend('show');

% さらに ROIごとにラベル(A1, B2...)を文字として重ねる
for i = 1:N
   text(x_coords(i), y_coords(i), finalClusterLabel{i}, ...
       'Color','w','FontSize',8,'HorizontalAlignment','center', ...
       'VerticalAlignment','middle','FontWeight','bold');
end


%--------------------------------------------------------------------------
%「クラスターごとに大きく枠線を囲む」処理
%--------------------------------------------------------------------------
% たとえば finalClusterLabel（A1, A2, B1, B2...）をキーにして領域をまとめます。
uniqueFinalLabels = unique(finalClusterLabel);
nLabels = length(uniqueFinalLabels);

% 輪郭を引くときの色を決める(クラスター数ぶん)
%  - ここではラベルごとに別の色を割り当てる例
cmap2 = lines(nLabels);

for cIdx = 1:nLabels
    thisLabel = uniqueFinalLabels{cIdx};  % 例: "A1"
    idxCluster = find(strcmp(finalClusterLabel, thisLabel));
    
    if length(idxCluster) < 3
        % ROIが2点以下ならconvhull or boundaryは作りづらいのでスキップ
        continue;
    end
    
    % すべてのROI座標をまとめる
    x_c = x_coords(idxCluster);
    y_c = y_coords(idxCluster);
    points = [x_c, y_c];

    % boundaryまたはconvhullを使って領域を求める
    % boundary(points(:,1), points(:,2), alpha) アルファ値は0~1あたりを調整
    % convhull(points(:,1), points(:,2)) でもOK（単純凸包）
    
    k = boundary(points(:,1), points(:,2), 0.8);  % 0.8は例: 曲線の滑らかさ

    % 枠線だけを描画(FaceColor='none'で塗りつぶしなし, EdgeColorで色指定)
    plotColor = cmap2(cIdx,:);
    patch('XData', points(k,1), 'YData', points(k,2), ...
          'EdgeColor', plotColor, 'FaceColor', 'none', ...
          'LineWidth', 2, 'LineStyle','-');
      
    % クラスターラベル名を、領域の中心あたりに一括表示したい場合
    cx = mean(points(k,1));
    cy = mean(points(k,2));
    text(cx, cy, thisLabel, ...
        'Color', plotColor, 'FontSize', 10, 'FontWeight','bold', ...
        'HorizontalAlignment','center','VerticalAlignment','middle');
end


hold off;

% 例3) 結果をテーブルにまとめてCSV出力
%      ROI_number / x / y / coordClusterIdx / patternClusterIdx / finalLabel
Tout = table(roiNumber, x_coords, y_coords, ...
             clusterCoordIdx, clusterPatternIdx, finalClusterLabel, ...
             'VariableNames',{'ROI_number','x','y','CoordClusterID','PatternClusterID','FinalLabel'});
writetable(Tout, 'clustering_result.csv');
fprintf('Saved clustering_result.csv\n');

