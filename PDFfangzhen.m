%% 1. 数据准备 (模拟含宽度的多径)
clear; clc; close all;

% 模拟参数
n_delay = 200; 
n_time = 1000;
h_matrix = (randn(n_delay, n_time) + 1j*randn(n_delay, n_time)) * 0.01; % 噪声底

% 模拟一条有“宽度”的主径
% 假设脉冲成型导致能量扩散在 +/- 3 个点
t_pulse = -3:3;
pulse_shape = sinc(t_pulse/1.5); % 模拟成型滤波

% 在第 50 个采样点添加强径，并让它随时间瑞利衰落
for t = 1:n_time
    fading_coeff = (randn + 1j*randn) + 3; % 莱斯衰落系数
    % 将波形叠加进去，形成有"宽度"的径
    h_matrix(50-3:50+3, t) = h_matrix(50-3:50+3, t) + (fading_coeff * pulse_shape).';
end

%% 2. 正确统计方法：局部峰值提取
h_abs = abs(h_matrix);

% 参数设置
threshold_val = max(h_abs(:)) * 10^(-20/20); % -20dB 阈值
min_peak_dist = 3; % 最小峰值距离 (防止把同一个sinc波形的旁瓣当成多径)

valid_peaks = []; % 用于存储提取出来的样本

% --- 核心循环 ---
for t = 1:n_time
    % 取出当前时刻的一列数据
    current_cir = h_abs(:, t);
    
    % 使用 findpeaks 寻找局部极大值
    % 'MinPeakHeight': 过滤噪声
    % 'MinPeakDistance': 过滤同一路径的宽度效应
    [pks, ~] = findpeaks(current_cir, ...
                         'MinPeakHeight', threshold_val, ...
                         'MinPeakDistance', min_peak_dist);
    
    % 将找到的峰值加入统计池
    valid_peaks = [valid_peaks; pks];
end

%% 3. 错误统计方法 (作为对比)
% 错误做法：直接利用阈值截断所有点 (包含波形的旁瓣)
binary_mask = h_abs > threshold_val;
wrong_samples = h_abs(binary_mask); 

%% 4. 绘图对比
figure('Position', [100, 100, 1000, 500]);

% 归一化
data_correct = valid_peaks / sqrt(mean(valid_peaks.^2));
data_wrong = wrong_samples / sqrt(mean(wrong_samples.^2));

% 绘制 PDF 对比
x_axis = linspace(0, 4, 100);

subplot(1, 2, 1);
histogram(data_correct, 50, 'Normalization', 'pdf'); hold on;
% 拟合 Rician
pd_corr = fitdist(data_correct, 'Rician');
plot(x_axis, pdf(pd_corr, x_axis), 'r-', 'LineWidth', 2);
title('方法A：局部峰值提取 (正确)');
subtitle(['拟合 K 因子: s=' num2str(pd_corr.s, '%.2f')]);
legend('样本直方图', 'Rician 拟合');
grid on;

subplot(1, 2, 2);
histogram(data_wrong, 50, 'Normalization', 'pdf'); hold on;
% 拟合 Rician
pd_wrong = fitdist(data_wrong, 'Rician');
plot(x_axis, pdf(pd_wrong, x_axis), 'r-', 'LineWidth', 2);
title('方法B：包含宽度内所有点 (错误)');
subtitle(['拟合 K 因子: s=' num2str(pd_wrong.s, '%.2f')]);
legend('样本直方图', 'Rician 拟合');
grid on;

fprintf('--- 结果分析 ---\n');
fprintf('正确方法的样本数: %d\n', length(valid_peaks));
fprintf('错误方法的样本数: %d (因为包含了波形旁瓣)\n', length(wrong_samples));
