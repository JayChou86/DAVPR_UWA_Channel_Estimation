%% =================================================================
%  visualize_cir.m — 修改后的CIR可视化（自动多径提取对比）
%  展示: 原始CIR热力图 + 自动提取后的有效多径 + PDP峰值标注
% =================================================================
cc
addpath('utils');

depths = [475, 624, 824, 924];
depth_labels = {'Node A (475.8m)', 'Node B (624.8m)', 'Node C (824.8m)', 'Node D (924.8m)'};
data_files = {
    '时变冲激响应475_0.1_3.mat',
    '时变冲激响应624_0.1_3.mat',
    '时变冲激响应824_0.1_3.mat',
    '时变冲激响应924_0.1_3.mat'
    };
fs = 16000;

%% ================================================================
%  Fig 1: 4深度原始CIR热力图对比 (2x2)
% ================================================================
figure('Position', [30, 30, 1100, 800]);

for d = 1:4
    load(data_files{d});
    
    subplot(2, 2, d);
    
    % 去除负时延，限制500ms显示范围
    seg_length = (size(h_matrix,2) + 1)/2;
    h_disp = h_matrix(:, seg_length:end);
    delay_ms = (0:size(h_disp,2)-1)' / fs * 1000;
    time_sec_full = (0:size(h_matrix,1)-1)' * 0.1;
    
    imagesc(delay_ms, time_sec_full, abs(h_disp));
    ylabel('Time (s)'); xlabel('Delay (ms)');
    title(sprintf('%s — Raw CIR', depth_labels{d}), 'FontSize', 11);
    colorbar; axis xy; colormap('jet');
    xlim([0, 500]);
    caxis([0, 0.3]);
end
sgtitle('Measured Time-Varying CIR: Four UASN Deployment Depths', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'cir_fig1_raw_cir_4depth.png');
saveas(gcf, 'cir_fig1_raw_cir_4depth.fig');
fprintf('Fig CIR-1: 原始CIR热力图 → cir_fig1_raw_cir_4depth.png\n');

%% ================================================================
%  Fig 2: 自动多径提取过程展示 (每个深度一行: CIR + PDP + 峰值)
% ================================================================
figure('Position', [30, 30, 1400, 900]);

for d = 1:4
    load(data_files{d});
    
    % ----- 预处理 (与Channel_review一致) -----
    avg_pdp = mean(abs(h_matrix).^2, 1);
    [~, peak_idx] = max(avg_pdp);
    start_idx = max(1, peak_idx - 100);
    end_idx = min(size(h_matrix,2), peak_idx + 2500);
    h_cut1 = h_matrix(:, start_idx:end_idx);
    
    [m, n] = size(h_cut1);
    result_matrix = zeros(m, n);
    for row = 1:m
        row_data = h_cut1(row, 1:min(150, n));
        [~, max_val] = max(row_data);
        shift_amount = 100 - max_val;
        if shift_amount > 0
            result_matrix(row, shift_amount+1:end) = h_cut1(row, 1:end-shift_amount);
        elseif shift_amount < 0
            shift_amount = abs(shift_amount);
            result_matrix(row, 1:end-shift_amount) = h_cut1(row, shift_amount+1:end);
        else
            result_matrix(row, :) = h_cut1(row, :);
        end
    end
    
    % 自动多径提取
    [h_extracted, keep_indices, pdp_norm] = extract_multipath_by_threshold(...
        abs(result_matrix), -20, 5);
    
    time_sec = (0:size(result_matrix,1)-1)' * 0.1;
    delay_full = (0:size(result_matrix,2)-1)' / fs * 1000;
    delay_extr = (0:size(h_extracted,2)-1)' / fs * 1000;
    
    % 列1: 同步后的CIR热力图(截取区域)
    subplot(4, 3, (d-1)*3 + 1);
    imagesc(delay_full, time_sec, abs(result_matrix));
    ylabel(sprintf('Node %s\nTime (s)', char('A'+d-1)), 'FontSize', 9);
    if d == 1, title('Synchronized CIR', 'FontSize', 11); end
    colorbar; axis xy; colormap('jet'); caxis([0, 0.2]);
    
    % 列2: 提取后的有效多径CIR
    subplot(4, 3, (d-1)*3 + 2);
    imagesc(delay_extr, time_sec, abs(h_extracted));
    if d == 1, title(sprintf('Auto-Extracted (%d paths)', length(keep_indices)), 'FontSize', 11); end
    colorbar; axis xy; colormap('jet'); caxis([0, 0.2]);
    
    % 列3: PDP + 检测到的峰值标注
    subplot(4, 3, (d-1)*3 + 3);
    pdp_dB = 10*log10(pdp_norm + eps);
    plot(delay_full, pdp_dB, 'b-', 'LineWidth', 1);
    hold on;
    % 标注检测到的峰值
    plot(delay_full(keep_indices), pdp_dB(keep_indices), 'ro', ...
        'MarkerSize', 6, 'MarkerFaceColor', 'r', 'LineWidth', 1);
    % -20dB 阈值线
    yline(-20, '--r', '-20dB', 'LineWidth', 0.8, 'Alpha', 0.5);
    xlabel('Delay (ms)'); ylabel('Power (dB)');
    if d == 1, title(sprintf('PDP with %d Detected Peaks', length(keep_indices)), 'FontSize', 11); end
    xlim([0, max(delay_full)]);
    ylim([-40, 5]);
    grid on;
    
    fprintf('  Node %s: %d paths detected (from %d columns → %d columns)\n', ...
        char('A'+d-1), length(keep_indices), size(result_matrix,2), size(h_extracted,2));
end

sgtitle('Automatic Multipath Extraction: Before vs After', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'cir_fig2_extraction_process.png');
saveas(gcf, 'cir_fig2_extraction_process.fig');
fprintf('Fig CIR-2: 多径提取过程 → cir_fig2_extraction_process.png\n');

%% ================================================================
%  Fig 3: 新旧方法对比 — 手工选径 vs 自动提取
%  用475m为例说明改进
% ================================================================
load('时变冲激响应475_0.1_3.mat');

% 预处理 (同上)
avg_pdp = mean(abs(h_matrix).^2, 1);
[~, peak_idx] = max(avg_pdp);
start_idx = max(1, peak_idx - 100);
end_idx = min(size(h_matrix,2), peak_idx + 2500);
h_cut1 = h_matrix(:, start_idx:end_idx);
[m, n] = size(h_cut1);
result_matrix = zeros(m, n);
for row = 1:m
    row_data = h_cut1(row, 1:min(150, n));
    [~, max_val] = max(row_data);
    shift_amount = 100 - max_val;
    if shift_amount > 0
        result_matrix(row, shift_amount+1:end) = h_cut1(row, 1:end-shift_amount);
    elseif shift_amount < 0
        shift_amount = abs(shift_amount);
        result_matrix(row, 1:end-shift_amount) = h_cut1(row, shift_amount+1:end);
    else
        result_matrix(row, :) = h_cut1(row, :);
    end
end

% 自动提取
[h_auto, keep_auto, pdp_auto] = extract_multipath_by_threshold(abs(result_matrix), -20, 5);

% 模拟旧方法: 手工选2条径 (原Channel_review_475的做法)
h_manual = zeros(74, 400);
h_manual(:, [10, 361]) = abs(result_matrix(:, [100, 451]));

time_sec = (0:size(result_matrix,1)-1)' * 0.1;

figure('Position', [30, 30, 1200, 450]);

% 旧方法
subplot(1, 3, 1);
imagesc((0:399)/fs*1000, time_sec, abs(h_manual));
ylabel('Time (s)'); xlabel('Delay (ms)');
title(sprintf('OLD: Manual Selection\n(2 paths, hand-picked)'), 'FontSize', 12);
colorbar; axis xy; colormap('jet');

% 新方法
subplot(1, 3, 2);
delay_auto = (0:size(h_auto,2)-1)' / fs * 1000;
imagesc(delay_auto, time_sec, abs(h_auto));
ylabel('Time (s)'); xlabel('Delay (ms)');
title(sprintf('NEW: Auto Extraction\n(%d paths, -20dB threshold)', length(keep_auto)), 'FontSize', 12);
colorbar; axis xy; colormap('jet');

% PDP对比
subplot(1, 3, 3);
pdp_original_dB = 10*log10(pdp_auto + eps);
delay_full = (0:length(pdp_auto)-1)' / fs * 1000;
plot(delay_full, pdp_original_dB, 'b-', 'LineWidth', 1.5); hold on;
% 标注自动提取的峰值
plot(delay_full(keep_auto), pdp_original_dB(keep_auto), 'go', ...
    'MarkerSize', 8, 'MarkerFaceColor', 'g');
% 标注旧方法手工选的2条径位置
old_indices = [100, 451];  % 对应原代码中 h_matrix_cut(:,[10,361]) = abs(result_matrix(:,[100,451]))
plot(delay_full(old_indices), pdp_original_dB(old_indices), 'rx', ...
    'MarkerSize', 12, 'LineWidth', 2);
yline(-20, '--r', '-20dB', 'LineWidth', 0.8);
xlabel('Delay (ms)'); ylabel('Power (dB)');
title('PDP: Hand-Picked vs Auto-Detected', 'FontSize', 12);
legend('PDP', sprintf('%d Auto-Detected', length(keep_auto)), ...
    '2 Hand-Picked (OLD)', '-20dB Threshold', 'Location', 'SouthEast');
grid on; ylim([-40, 5]);

sgtitle('Node A (475.8m): Old vs New Multipath Extraction Method', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'cir_fig3_old_vs_new.png');
saveas(gcf, 'cir_fig3_old_vs_new.fig');
fprintf('Fig CIR-3: 新旧方法对比 → cir_fig3_old_vs_new.png\n');

%% ================================================================
%  Fig 4: 4深度PDP对比 + 带温跃层标注
% ================================================================
figure('Position', [30, 30, 800, 500]);
colors = lines(4);
max_delay_ms = 300;  % 显示0-300ms

hold on;
for d = 1:4
    load(data_files{d});
    
    % 计算归一化PDP
    avg_pdp = mean(abs(h_matrix).^2, 1);
    seg_length = (size(h_matrix,2)+1)/2;
    pdp = avg_pdp(seg_length:end);
    delay_ms_all = (0:length(pdp)-1)' / fs * 1000;
    pdp_dB = 10*log10(pdp / max(pdp) + eps);
    
    % 截取到显示范围
    idx_range = delay_ms_all <= max_delay_ms;
    plot(delay_ms_all(idx_range), pdp_dB(idx_range), ...
        'Color', colors(d,:), 'LineWidth', 1.5);
end

xlabel('Delay (ms)', 'FontSize', 12);
ylabel('Normalized Power (dB)', 'FontSize', 12);
title('Power Delay Profiles Across Thermocline Depths', 'FontSize', 14);
legend(depth_labels, 'Location', 'NorthEast', 'FontSize', 10);
yline(-20, '--k', '-20 dB', 'LineWidth', 0.8, 'Alpha', 0.5);
grid on;
ylim([-40, 5]);
xlim([0, max_delay_ms]);

saveas(gcf, 'cir_fig4_pdp_comparison.png');
saveas(gcf, 'cir_fig4_pdp_comparison.fig');
fprintf('Fig CIR-4: PDP深度对比 → cir_fig4_pdp_comparison.png\n');

%% ================================================================
fprintf('\n========== CIR可视化图表生成完毕 ==========\n');
fprintf('  cir_fig1_raw_cir_4depth.png      - 4深度原始CIR热力图\n');
fprintf('  cir_fig2_extraction_process.png  - 自动多径提取过程\n');
fprintf('  cir_fig3_old_vs_new.png          - 新旧方法对比(475m为例)\n');
fprintf('  cir_fig4_pdp_comparison.png      - PDP深度对比\n');
