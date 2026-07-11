cc
folder = 'D:\A_Matlab_program\A_Matlab_program\时变信道估计代码\相关';              % 放 mat 文件的目录
files = dir(fullfile(folder, '*.mat'));
for k = 1:numel(files)
    fname = files(k).name;
    dataStruct = load(fullfile(folder, fname));
    % 假设每个 mat 文件里只有一个变量
    varNames = fieldnames(dataStruct);
    if numel(varNames) ~= 1
        warning('文件 %s 中变量数量 ≠ 1，已跳过。', fname);
        continue
    end
    % 根据文件名生成变量名（去掉扩展名并确保是合法的 MATLAB 标识符）
    baseName = matlab.lang.makeValidName(erase(fname, '.mat'));
    % 将数据保存到新的变量中
    assignin('base', baseName, dataStruct.(varNames{1}));
end

h475 = x______475;
h624 = x______624_2_.h_matrix;
h824 = x______824;
h924 = x______924;
h475_mean = abs(mean(h475,1)).^2 / max(abs(mean(h475,1))).^2;
h624_mean = abs(mean(h624,1)).^2 / max(abs(mean(h624,1))).^2;
h824_mean = abs(mean(h824,1)).^2 / max(abs(mean(h824,1))).^2;
h924_mean = abs(mean(h924,1)).^2 / max(abs(mean(h924,1))).^2;

delay_ms_all = (0:size(h924_mean, 2)-1)' / 16000 * 1000;     % 去除负时延的时延轴


means  = {h475_mean, h624_mean, h824_mean, h924_mean};
labels = {'475.8m', '624.8m', '824.8m', '924.8m'};
colors = lines(numel(means));
x = 1:numel(h475_mean);    % 样本索引
figure;
hold on;
for k = 1:numel(means)
    y = k * ones(size(x)); % 给每条线一个固定的“序列编号”高度，便于区分
    plot3(delay_ms_all, y, means{k}, 'LineWidth', 0.5, 'Color', colors(k,:));
end
hold off;
grid on;
xlabel('\fontname{宋体}时延\fontname{Times new roman}(ms)', ...
    'fontname','宋体', 'fontsize', 10, 'Interpreter', 'tex');


ylabel('深度');
zlabel('归一化幅度');
yticks(1:numel(labels));
yticklabels(labels);
view(45, 25);
% title('四组均值冲激响应');