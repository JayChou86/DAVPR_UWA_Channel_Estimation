cc
folder = 'D:\A_Matlab_program\A_Matlab_program\时变信道估计代码\滑动窗分段\五月份实测数据\mean_CIR';              % 放 mat 文件的目录
files = dir(fullfile(folder, '*.mat'));

for k = 1:numel(files)
    fname = files(k).name;
    % 加载mat文件中的所有变量
    dataStruct = load(fullfile(folder, fname));
    
    % 获取所有变量名
    varNames = fieldnames(dataStruct);
    
    if numel(varNames) > 3
        warning('文件 %s 中变量数量 > 3，只读取前3个变量。', fname);
    elseif numel(varNames) < 3
        warning('文件 %s 中变量数量 < 3，将读取所有变量。', fname);
    end
    
    % 遍历所有变量，将它们保存到工作空间
    for i = 1:min(3, numel(varNames))
        % 生成新的变量名：文件名_变量名
        baseName = matlab.lang.makeValidName([erase(fname, '.mat') '_' varNames{i}]);
        
        % 将数据保存到工作空间
        assignin('base', baseName, dataStruct.(varNames{i}));
        
        % 可选：显示读取的变量信息
        fprintf('已加载变量: %s (大小: %s)\n', baseName, mat2str(size(dataStruct.(varNames{i}))));
    end
end

num_segments = size(x______475_h_matrix,1);
step_duration = 0.1;
fs = 16000;
step_length = step_duration * fs;
seg_length = (size(x______475_h_matrix,2) + 1)/2;

% 显示总结信息
fprintf('\n总计读取了 %d 个文件\n', numel(files));
fprintf('每个文件最多读取了 3 个变量\n');

h475 = x______475_h_matrix(:, seg_length:end);
h624 = x______624_h_matrix(:, seg_length:end);
h824 = x______824_h_matrix(:, seg_length:end);
h924 = x______924_h_matrix(:, seg_length:end);
h475_mean = abs(mean(h475,1)).^2 / max(abs(mean(h475,1))).^2;
h624_mean = abs(mean(h624,1)).^2 / max(abs(mean(h624,1))).^2;
h824_mean = abs(mean(h824,1)).^2 / max(abs(mean(h824,1))).^2;
h924_mean = abs(mean(h924,1)).^2 / max(abs(mean(h924,1))).^2;

delay_ms = (0:seg_length-1)' / fs * 1000;


means  = {h475_mean, h624_mean, h824_mean, h924_mean};
labels = {'475.8', '624.8', '824.8', '924.8'};
% colors = lines(numel(means));
colors = [0.0039    0.3373    0.6000; 0.9804    0.7529    0.0588;...
    0.9529    0.4627    0.2902; 0.3725    0.7765    0.7882];
x = 1:numel(h475_mean);    % 样本索引
figure;
hold on;
% 示例代码，假设已经定义了means、delay_ms等变量
for k = 1:numel(means)
    y = k * ones(size(x)); % 给每条线一个固定的"序列编号"高度，便于区分
    
    % 找出delay_ms在0到500ms范围内的索引
    valid_indices = delay_ms >= 0 & delay_ms <= 500;
    
    % 只绘制0到500ms范围内的数据
    plot3(delay_ms(valid_indices), y(valid_indices), means{k}(valid_indices), ...
          'LineWidth', 0.5, 'Color', colors(k,:));
    
    % 如果需要保持x轴范围一致，可以添加hold on
    if k == 1
        hold on;
    end
end

% 设置x轴范围为0-500
xlim([0, 500]);

% 添加坐标轴标签和标题（根据需要）
xlabel('Delay (ms)');
ylabel('Depth (m)');
zlabel('Normalized Power');

% 如果需要网格
grid on;
hold off;
hold off;
grid on;
yticks(1:numel(labels));
yticklabels(labels);
view(45, 25);
% title('四组均值冲激响应');