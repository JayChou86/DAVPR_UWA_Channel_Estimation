%% =================================================================
%  generate_figures.m — 论文图表批量生成脚本
%  生成 UASN 深度感知接收机预配置方案的全部关键图表
% =================================================================
cc
addpath('utils');

%% 配置
depths = [475, 624, 824, 924];
depth_labels = {'Node A (475.8m)', 'Node B (624.8m)', 'Node C (824.8m)', 'Node D (924.8m)'};
short_labels = {'A:475m', 'B:624m', 'C:824m', 'D:924m'};
thermocline_labels = {'Upper', 'Within', 'Lower*', 'Below'};
data_files = {
    '时变冲激响应475_0.1_3.mat',
    '时变冲激响应624_0.1_3.mat',
    '时变冲激响应824_0.1_3.mat',
    '时变冲激响应924_0.1_3.mat'
    };
colors = [
    0.00 0.45 0.74;  % 蓝
    0.85 0.33 0.10;  % 橙
    0.93 0.69 0.13;  % 黄
    0.49 0.18 0.56;  % 紫
    ];

fs = 16000;
symbol_rate = 4000;
T_step = 0.1;

%% ================================================================
%  阶段1: 自动参数选择 — 收集4个深度的结果
% ================================================================
fprintf('========== 阶段1: 信道统计与自动参数选择 ==========\n');
all_stats = cell(4,1);
all_params = cell(4,1);
all_n_multipath = zeros(4,1);

for d = 1:4
    fprintf('处理深度 %s...\n', depth_labels{d});
    load(data_files{d});

    % 预处理
    avg_pdp = mean(abs(h_matrix).^2, 1);
    [~, peak_idx] = max(avg_pdp);
    start_idx = max(1, peak_idx - 100);
    end_idx = min(size(h_matrix,2), peak_idx + 2500);
    h_cut1 = h_matrix(:, start_idx:end_idx);
    clear h_matrix;

    [m, n] = size(h_cut1);
    result_matrix = zeros(m, n);
    for row = 1:m
        row_data = h_cut1(row, 1:min(150, n));
        [~, max_idx] = max(row_data);
        shift_amount = 100 - max_idx;
        if shift_amount > 0
            result_matrix(row, shift_amount+1:end) = h_cut1(row, 1:end-shift_amount);
        elseif shift_amount < 0
            shift_amount = abs(shift_amount);
            result_matrix(row, 1:end-shift_amount) = h_cut1(row, shift_amount+1:end);
        else
            result_matrix(row, :) = h_cut1(row, :);
        end
    end
    clear h_cut1;

    [h_cut, keep_indices] = extract_multipath_by_threshold(abs(result_matrix), -20, 5);
    clear result_matrix;
    all_n_multipath(d) = length(keep_indices);

    params = auto_parameter_selector(h_cut, fs, symbol_rate, T_step, depth_labels{d});
    all_params{d} = params;
    all_stats{d} = params.stats;
    clear h_cut;
end

%% ================================================================
%  Fig 1: 自动参数对比柱状图 (fb_taps, tau_rms, ff_taps)
% ================================================================
fprintf('\n========== 生成图表 ==========\n');

figure('Position', [50, 50, 900, 600]);

% 子图1: fb_taps + tau_rms
subplot(2,2,1);
fb_values = cellfun(@(p) p.eq_fb_taps, all_params);
tau_values = cellfun(@(s) s.tau_rms, all_stats);
yyaxis left;
b1 = bar(1:4, fb_values, 0.5, 'FaceColor', [0.3 0.5 0.8]);
ylabel('Feedback Taps (N_{fb})', 'FontSize', 11, 'Color', [0.3 0.5 0.8]);
yyaxis right;
plot(1:4, tau_values, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
ylabel('RMS Delay Spread (ms)', 'FontSize', 11, 'Color', 'r');
set(gca, 'XTick', 1:4, 'XTickLabel', short_labels);
xlabel('Sensor Node');
title('Feedback Taps vs. RMS Delay Spread', 'FontSize', 12);
grid on;

% 子图2: ff_taps + multipath count
subplot(2,2,2);
ff_values = cellfun(@(p) p.eq_ff_taps, all_params);
yyaxis left;
b2 = bar(1:4, ff_values, 0.5, 'FaceColor', [0.3 0.7 0.4]);
ylabel('Forward Taps (N_{ff})', 'FontSize', 11, 'Color', [0.3 0.7 0.4]);
ylim([0, 70]);
yyaxis right;
plot(1:4, all_n_multipath, 'm-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'm');
ylabel('Significant Multipaths', 'FontSize', 11, 'Color', 'm');
set(gca, 'XTick', 1:4, 'XTickLabel', short_labels);
xlabel('Sensor Node');
title('Forward Taps vs. Multipath Richness', 'FontSize', 12);
grid on;

% 子图3: lambda + coherence time
subplot(2,2,3);
lambda_values = cellfun(@(p) p.rls_forget, all_params);
tc_values = cellfun(@(s) s.T_c * 1e3, all_stats);
yyaxis left;
plot(1:4, lambda_values, 'b-o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
ylabel('RLS Forgetting Factor \lambda', 'FontSize', 11, 'Color', 'b');
ylim([0.9985, 0.9995]);
yyaxis right;
plot(1:4, tc_values, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
ylabel('Coherence Time T_c (ms)', 'FontSize', 11, 'Color', 'r');
set(gca, 'XTick', 1:4, 'XTickLabel', short_labels);
xlabel('Sensor Node');
title('Forgetting Factor vs. Coherence Time', 'FontSize', 12);
grid on;

% 子图4: 复杂度对比
subplot(2,2,4);
rls_ops = cellfun(@(p) p.rls_flops_per_sym, all_params);
lms_ops = cellfun(@(p) p.lms_flops_per_sym, all_params);
bar_data = [rls_ops(:), lms_ops(:)];
b = bar(bar_data);
b(1).FaceColor = [0.2 0.4 0.8];
b(2).FaceColor = [0.9 0.5 0.3];
set(gca, 'XTick', 1:4, 'XTickLabel', short_labels);
ylabel('Multiplications per Symbol');
title('Computational Complexity: RLS vs LMS', 'FontSize', 12);
legend('RLS-DFE (O(N^2))', 'LMS-DFE (O(N))', 'Location', 'NorthWest');
grid on;

sgtitle('Depth-Aware Parameter Selection for UASN Nodes', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'fig1_parameter_selection.png');
saveas(gcf, 'fig1_parameter_selection.fig');
fprintf('  Fig 1: 参数选择对比图 → fig1_parameter_selection.png\n');

%% ================================================================
%  Fig 2: 信道统计特征雷达图/条形图
% ================================================================
figure('Position', [50, 50, 900, 400]);

subplot(1,3,1);
% 相干时间
tc_data = cellfun(@(s) s.T_c * 1e3, all_stats);
bar(1:4, tc_data, 0.6, 'FaceColor', [0.3 0.6 0.9]);
set(gca, 'XTick', 1:4, 'XTickLabel', short_labels);
ylabel('Coherence Time (ms)');
title('Coherence Time T_c', 'FontSize', 12);
% 标注准静态
yl = ylim;
text(0.5, yl(2)*0.95, 'All nodes: quasi-static (T_c > 3s)', ...
    'FontSize', 9, 'Color', [0.4 0.4 0.4]);
grid on;

subplot(1,3,2);
% 相干带宽
bc_data = cellfun(@(s) s.B_c, all_stats);
bar(1:4, bc_data, 0.6, 'FaceColor', [0.9 0.5 0.2]);
set(gca, 'XTick', 1:4, 'XTickLabel', short_labels);
ylabel('Coherence Bandwidth (Hz)');
title('Coherence Bandwidth B_c', 'FontSize', 12);
text(0.5, max(bc_data)*0.95, 'Narrowband: B_c = 4-6 Hz', ...
    'FontSize', 9, 'Color', [0.4 0.4 0.4]);
grid on;

subplot(1,3,3);
% 多普勒扩展
fd_data = cellfun(@(s) s.f_d_rms, all_stats);
bar(1:4, fd_data, 0.6, 'FaceColor', [0.5 0.3 0.7]);
set(gca, 'XTick', 1:4, 'XTickLabel', short_labels);
ylabel('RMS Doppler Spread (Hz)');
title('Doppler Spread f_d^{rms}', 'FontSize', 12);
text(0.5, max(fd_data)*0.95, 'Minimal Doppler (< 0.6 Hz)', ...
    'FontSize', 9, 'Color', [0.4 0.4 0.4]);
grid on;

sgtitle('Channel Statistics Across Thermocline Depths', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'fig2_channel_statistics.png');
saveas(gcf, 'fig2_channel_statistics.fig');
fprintf('  Fig 2: 信道统计图 → fig2_channel_statistics.png\n');

%% ================================================================
%  Fig 3: 部署配置表 (热力图/矩阵)
% ================================================================
figure('Position', [50, 50, 700, 500]);

% 构建对比矩阵
config_matrix = [
    cellfun(@(p) p.eq_ff_taps, all_params)';
    cellfun(@(p) p.eq_fb_taps, all_params)';
    cellfun(@(s) s.tau_rms, all_stats)';
    cellfun(@(s) s.T_c*1e3, all_stats)';
    cellfun(@(s) s.B_c, all_stats)';
    cellfun(@(s) s.f_d_rms*1000, all_stats)';  % mHz
    cellfun(@(p) p.rls_forget*10000-9990, all_params)';  % scaled
    ];

% 归一化各行用于热力图
config_norm = (config_matrix - min(config_matrix,[],2)) ./ ...
    (max(config_matrix,[],2) - min(config_matrix,[],2) + eps);

imagesc(config_norm);
colormap('parula');
colorbar;
set(gca, 'XTick', 1:4, 'XTickLabel', short_labels);
set(gca, 'YTick', 1:7, 'YTickLabel', ...
    {'N_{ff}', 'N_{fb}', '\tau_{rms} (ms)', 'T_c (ms)', 'B_c (Hz)', 'f_d (mHz)', '\lambda (scaled)'});
title('Normalized Configuration Matrix Across Nodes', 'FontSize', 14);
xlabel('Sensor Node');

% 标注数值
for i = 1:7
    for j = 1:4
        if i <= 2
            text(j, i, sprintf('%d', round(config_matrix(i,j))), ...
                'HorizontalAlign', 'center', 'FontSize', 9, 'FontWeight', 'bold');
        else
            text(j, i, sprintf('%.1f', config_matrix(i,j)), ...
                'HorizontalAlign', 'center', 'FontSize', 9);
        end
    end
end

saveas(gcf, 'fig3_configuration_matrix.png');
saveas(gcf, 'fig3_configuration_matrix.fig');
fprintf('  Fig 3: 配置矩阵热力图 → fig3_configuration_matrix.png\n');

%% ================================================================
%  Fig 4: 深度 vs 信道特征 + 温跃层标注
% ================================================================
figure('Position', [50, 50, 800, 550]);

subplot(2,2,1);
plot(depths, tau_values, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 12, 'MarkerFaceColor', 'b');
xlabel('Depth (m)'); ylabel('RMS Delay Spread (ms)');
title('RMS Delay Spread vs. Depth');
grid on; xlim([400, 1000]);
% 标注温跃层区间
hold on; yl = ylim;
fill([600 850 850 600], [yl(1) yl(1) yl(2) yl(2)], ...
    [1 0.9 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
text(725, yl(2)*0.95, 'Thermocline', 'HorizontalAlign', 'center', 'FontSize', 9);

subplot(2,2,2);
plot(depths, tc_data, 'r-s', 'LineWidth', 2.5, 'MarkerSize', 12, 'MarkerFaceColor', 'r');
xlabel('Depth (m)'); ylabel('Coherence Time (ms)');
title('Coherence Time vs. Depth');
grid on; xlim([400, 1000]);
hold on; yl = ylim;
fill([600 850 850 600], [yl(1) yl(1) yl(2) yl(2)], ...
    [1 0.9 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

subplot(2,2,3);
plot(depths, bc_data, 'g-^', 'LineWidth', 2.5, 'MarkerSize', 12, 'MarkerFaceColor', 'g');
xlabel('Depth (m)'); ylabel('Coherence Bandwidth (Hz)');
title('Coherence Bandwidth vs. Depth');
grid on; xlim([400, 1000]);
hold on; yl = ylim;
fill([600 850 850 600], [yl(1) yl(1) yl(2) yl(2)], ...
    [1 0.9 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

subplot(2,2,4);
plot(depths, fb_values, 'm-d', 'LineWidth', 2.5, 'MarkerSize', 12, 'MarkerFaceColor', 'm');
xlabel('Depth (m)'); ylabel('Optimal Feedback Taps');
title('Auto-Selected Feedback Taps vs. Depth');
grid on; xlim([400, 1000]);
hold on; yl = ylim;
fill([600 850 850 600], [yl(1) yl(1) yl(2) yl(2)], ...
    [1 0.9 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

sgtitle('Channel Characteristics Across Thermocline Depths', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'fig4_depth_analysis.png');
saveas(gcf, 'fig4_depth_analysis.fig');
fprintf('  Fig 4: 深度分析图 → fig4_depth_analysis.png\n');

%% ================================================================
%  Fig 5: 自动推导公式图 (概念示意图)
% ================================================================
figure('Position', [50, 50, 800, 400]);

% 左侧: 参数推导流程示意
subplot(1,2,1);
axis off; hold on;
% 用文字块示意
text(0.5, 0.95, 'Parameter Derivation Pipeline', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlign', 'center');

% CIR → Stats
rectangle('Position', [0.1, 0.7, 0.8, 0.15], 'Curvature', 0.1, ...
    'FaceColor', [0.8 0.9 1], 'EdgeColor', [0.2 0.4 0.8], 'LineWidth', 1.5);
text(0.5, 0.775, 'Measured CIR h(\tau, t)', 'FontSize', 12, 'HorizontalAlign', 'center');

arrow_y = 0.68;
plot([0.5, 0.5], [arrow_y, arrow_y-0.05], 'k-', 'LineWidth', 1.5);
plot([0.48, 0.5, 0.52], [arrow_y-0.03, arrow_y, arrow_y-0.03], 'k-', 'LineWidth', 1.5);

% Stats
rectangle('Position', [0.1, 0.48, 0.8, 0.15], 'Curvature', 0.1, ...
    'FaceColor', [0.9 1 0.8], 'EdgeColor', [0.2 0.7 0.3], 'LineWidth', 1.5);
text(0.5, 0.555, 'Channel Statistics (\tau_{rms}, T_c, B_c, f_d)', 'FontSize', 11, 'HorizontalAlign', 'center');

plot([0.5, 0.5], [0.46, 0.41], 'k-', 'LineWidth', 1.5);
plot([0.48, 0.5, 0.52], [0.43, 0.46, 0.43], 'k-', 'LineWidth', 1.5);

% Rules
rectangle('Position', [0.05, 0.25, 0.4, 0.15], 'Curvature', 0.1, ...
    'FaceColor', [1 0.95 0.7], 'EdgeColor', [0.9 0.6 0.1], 'LineWidth', 1.5);
text(0.25, 0.325, 'N_{fb} = f(\tau_{rms})', 'FontSize', 10, 'HorizontalAlign', 'center');

rectangle('Position', [0.55, 0.25, 0.4, 0.15], 'Curvature', 0.1, ...
    'FaceColor', [1 0.85 0.85], 'EdgeColor', [0.9 0.3 0.2], 'LineWidth', 1.5);
text(0.75, 0.325, '\lambda = f(T_c, f_d)', 'FontSize', 10, 'HorizontalAlign', 'center');

plot([0.25, 0.25], [0.23, 0.18], 'k-', 'LineWidth', 1.5);
plot([0.75, 0.75], [0.23, 0.18], 'k-', 'LineWidth', 1.5);
plot([0.23, 0.25, 0.27], [0.20, 0.23, 0.20], 'k-', 'LineWidth', 1.5);
plot([0.73, 0.75, 0.77], [0.20, 0.23, 0.20], 'k-', 'LineWidth', 1.5);

% Output
rectangle('Position', [0.1, 0.02, 0.8, 0.15], 'Curvature', 0.1, ...
    'FaceColor', [0.7 0.9 0.7], 'EdgeColor', [0.1 0.6 0.2], 'LineWidth', 1.5);
text(0.5, 0.095, 'Deployment Config: (N_{ff}, N_{fb}, \lambda)_{per node}', ...
    'FontSize', 12, 'HorizontalAlign', 'center', 'FontWeight', 'bold');

% 右侧: 公式
subplot(1,2,2);
axis off; hold on;
text(0.5, 0.95, 'Derivation Rules', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlign', 'center');

text(0.1, 0.78, 'Feedback Taps:', 'FontSize', 12, 'FontWeight', 'bold');
text(0.15, 0.70, 'N_{fb} = \lceil 3 \cdot \tau_{rms} \cdot R_s \rceil', ...
    'FontSize', 13, 'Color', [0.2 0.4 0.8]);
text(0.15, 0.64, '\tau_{rms}: RMS delay spread (s)', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);
text(0.15, 0.60, 'R_s: symbol rate (sps)', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);

text(0.1, 0.48, 'Forgetting Factor:', 'FontSize', 12, 'FontWeight', 'bold');
text(0.15, 0.40, '\lambda = exp(-T_{sym} / (10 \cdot T_c))', ...
    'FontSize', 13, 'Color', [0.8 0.3 0.2]);
text(0.15, 0.34, 'T_{sym}: symbol duration (s)', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);
text(0.15, 0.30, 'T_c: coherence time (s)', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);

text(0.1, 0.18, 'Forward Taps:', 'FontSize', 12, 'FontWeight', 'bold');
text(0.15, 0.10, 'N_{ff} = 3 \cdot \tau_{pre} \cdot R_s', ...
    'FontSize', 13, 'Color', [0.3 0.7 0.3]);
text(0.15, 0.04, '\tau_{pre}: pre-cursor delay @ -10dB (s)', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);

saveas(gcf, 'fig5_derivation_rules.png');
saveas(gcf, 'fig5_derivation_rules.fig');
fprintf('  Fig 5: 推导规则图 → fig5_derivation_rules.png\n');

%% ================================================================
%  Fig 6: 汇总部署配置表 (文本方式呈现)
% ================================================================
figure('Position', [50, 50, 750, 300]);
axis off; hold on;

text(0.5, 0.95, 'Deployment Configuration Table for UASN Nodes', ...
    'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlign', 'center');

% 表头
headers = {'Depth(m)', 'Layer', '\tau_{rms}(ms)', 'T_c(s)', 'B_c(Hz)', ...
    'f_d(Hz)', 'N_{ff}', 'N_{fb}', '\lambda', '#Paths'};
col_widths = [0.09, 0.09, 0.11, 0.08, 0.08, 0.08, 0.07, 0.07, 0.09, 0.07];
col_starts = cumsum([0.02, col_widths(1:end-1)]);

% 画表头
y_header = 0.82;
for c = 1:length(headers)
    text(col_starts(c)+col_widths(c)/2, y_header, headers{c}, ...
        'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlign', 'center');
end
plot([0.02, 0.98], [y_header-0.04, y_header-0.04], 'k-', 'LineWidth', 1.5);

% 画数据行
for d = 1:4
    y_row = 0.78 - d*0.12;
    row_color = colors(d,:);
    row_data = {sprintf('%d', depths(d)), thermocline_labels{d}, ...
        sprintf('%.1f', all_stats{d}.tau_rms), ...
        sprintf('%.1f', all_stats{d}.T_c), ...
        sprintf('%.1f', all_stats{d}.B_c), ...
        sprintf('%.3f', all_stats{d}.f_d_rms), ...
        sprintf('%d', all_params{d}.eq_ff_taps), ...
        sprintf('%d', all_params{d}.eq_fb_taps), ...
        sprintf('%.4f', all_params{d}.rls_forget), ...
        sprintf('%d', all_n_multipath(d))};

    for c = 1:length(row_data)
        text(col_starts(c)+col_widths(c)/2, y_row, row_data{c}, ...
            'FontSize', 9, 'HorizontalAlign', 'center', 'Color', row_color);
    end
    % 标注瓶颈行
    if depths(d) == 824
        text(0.99, y_row, '← Bottleneck', 'FontSize', 9, 'Color', 'r', ...
            'FontWeight', 'bold', 'HorizontalAlign', 'right');
    end
end

saveas(gcf, 'fig6_deployment_table.png');
saveas(gcf, 'fig6_deployment_table.fig');
fprintf('  Fig 6: 部署配置表 → fig6_deployment_table.png\n');

%% ================================================================
fprintf('\n========== 全部图表生成完毕 ==========\n');
fprintf('生成文件:\n');
fprintf('  fig1_parameter_selection.png    - 参数选择对比\n');
fprintf('  fig2_channel_statistics.png     - 信道统计特征\n');
fprintf('  fig3_configuration_matrix.png   - 配置矩阵热力图\n');
fprintf('  fig4_depth_analysis.png         - 深度分析\n');
fprintf('  fig5_derivation_rules.png       - 推导规则\n');
fprintf('  fig6_deployment_table.png       - 部署配置表\n');

%% 保存数据供 network_ber_comparison 使用
save('figures_data.mat', 'all_params', 'all_stats', 'all_n_multipath', ...
    'depths', 'depth_labels', 'short_labels', 'thermocline_labels');
disp('数据已保存到 figures_data.mat');
