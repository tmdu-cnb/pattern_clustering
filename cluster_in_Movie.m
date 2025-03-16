%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cluster_in_Movie.m
%
% 【概要】
%  1) TIF動画を読み込み、指定したフレーム範囲だけ処理
%  2) clustering_result.csv (あるいは clustering.m の出力) を読み込み、
%     - ROIごとの (x, y) 座標
%     - クラスターラベル (A1, B2, ...)
%     を取得
%  3) クラスターラベルごとに ROI 座標をまとめて boundary() し、枠線を取得
%  4) 各フレームに対して、クラスター枠線を `insertShape` で上書きし、新しい TIF に保存
%
%  ※ 既存の ROI_in_Movie.m では "ROI単位" のアウトラインを描画していましたが、
%    ここでは "クラスター単位" でまとめた領域を囲うイメージです。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cluster_in_Movie()

% -------------------------
% 動画ファイル (TIF)
% -------------------------
tifFile = 'Drd1_1_AAV-L7-RFP_L7-j7s_P9.tif';  % 読み込むTIF動画ファイル名
info = imfinfo(tifFile);
numFrames = numel(info);
fprintf('TIF frames: %d\n', numFrames);

% 例: 全フレームの1/3だけ処理など
startFrame = 1;
endFrame   = floor(numFrames / 3);
selectedFrames = startFrame:endFrame;
fprintf('Processing frames %d to %d.\n', startFrame, endFrame);

% -------------------------
% クラスタリング結果のCSVを読み込み
% clustering_result.csv は clustering.m で保存したものを想定
%   Columns例: ROI_number / x / y / CoordClusterID / PatternClusterID / FinalLabel
% -------------------------
clusterCSV = 'clustering_result.csv';
TC = readtable(clusterCSV);
roiNumber = TC.ROI_number;    % ROI番号 (Suite2P上の番号 or ユーザ定義)
x_coords  = TC.x;
y_coords  = TC.y;
labels    = TC.FinalLabel;    % A1, B2, ... など

% ROI数
N = length(roiNumber);

% -------------------------
% クラスターごとの枠線 (boundary) を事前に計算
% -------------------------
uniqueLabels = unique(labels);
nClusters = length(uniqueLabels);

% カラー設定 (cluster毎の色を lines で生成)
cmap = uint8(255 * lines(nClusters));

% クラスター領域のアウトラインを polygon として保存しておく
clusterPolygons = cell(nClusters, 1);

for cIdx = 1:nClusters
    thisLabel = uniqueLabels{cIdx};
    idxC = find(strcmp(labels, thisLabel));
    if length(idxC) < 3
        % ROIが2点以下なら輪郭を作れないのでスキップ
        clusterPolygons{cIdx} = [];
        continue;
    end
    
    % すべてのROI座標をまとめる
    x_c = x_coords(idxC);
    y_c = y_coords(idxC);
    points = [x_c, y_c];
    
    % boundary でアウトラインを取得
    k = boundary(points(:,1), points(:,2), 0.8);  % alpha=0.8 は例
    
    % boundaryで得た頂点座標を保存 (Polygon用に x1,y1,x2,y2,... の形に並べる)
    outlineXY = points(k,:);  % [Nx2]
    
    % insertShape は (x1,y1, x2,y2, x3,y3,...) のフォーマットが必要
    % reshape で 1x(2*N) に変換
    clusterPolygons{cIdx} = reshape(outlineXY', 1, []);
end

% -------------------------
% 出力用TIF
% -------------------------
outputTifFile = 'output_cluster_video.tif';

% フレーム処理時間を記録する配列
frameTimes = zeros(1, length(selectedFrames));

% パラレル処理 (必要ならパラレルプール open)
parfor i = 1:length(selectedFrames)
    frameIdx = selectedFrames(i);
    frameTimer = tic;  % 計測開始
    
    % フレーム読み込み (グレースケール)
    imgGray = imread(tifFile, frameIdx);
    
    % RGB 3チャネルに拡張
    if ndims(imgGray) == 2
        imgRGB = repmat(imgGray, [1,1,3]);
    else
        imgRGB = imgGray;
    end
    
    % クラスター枠線を描画 (各クラスターごとに insertShape)
    for cIdx = 1:nClusters
        polyXY = clusterPolygons{cIdx};
        if isempty(polyXY)
            continue; % ROIが2点以下等で枠がない
        end
        
        colorThisCluster = cmap(cIdx, :);  % [R G B]
        
        imgRGB = insertShape(imgRGB, 'Polygon', polyXY, ...
            'Color', colorThisCluster, 'LineWidth', 2);
    end
    
    % TIFに書き込み
    if i == 1
        imwrite(imgRGB, outputTifFile, 'WriteMode','overwrite','Compression','none');
    else
        imwrite(imgRGB, outputTifFile, 'WriteMode','append','Compression','none');
    end
    
    frameTimes(i) = toc(frameTimer);
    fprintf('Frame %d processed in %.2f sec\n', frameIdx, frameTimes(i));
end

% 処理時間を保存
writematrix(frameTimes, 'frame_processing_times.csv');
disp(['Output TIF: ', outputTifFile]);
disp('Frame times saved to frame_processing_times.csv');

end
