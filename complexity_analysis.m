%% =================================================================
%  complexity_analysis.m — 计算复杂度量化对比
%  目的: 对比 RLS-DFE vs LMS-DFE vs 连续信道估计 的计算开销
%  回答审稿人R2#4: "能效优势未量化，无计算复杂度或能耗实测"
% =================================================================
cc

%% 配置
fs = 16000;
symbol_rate = 4000;
T_sym = 1 / symbol_rate;

% 典型参数范围 (来自各深度的自动推导结果)
scenarios = {
    % {label,          eq_ff, eq_fb,  forget}
    'Node A (475m)',     60,    90,   0.995;
    'Node B (624m)',     60,   400,   0.998;
    'Node C (824m)',     60,   450,   0.998;
    'Node D (924m)',     60,   200,   0.996;
    };

% LMS 通常需要更多抽头来匹配 RLS 性能
lms_tap_factor = 2;  % LMS 需要约2倍抽头数

%% 计算每种场景的复杂度
fprintf('========== 计算复杂度对比分析 ==========\n\n');
fprintf('%-20s | %-8s | %-8s | %-8s | %-8s | %-8s\n', ...
    'Scenario', 'N_RLS', 'N_LMS', 'RLS ops', 'LMS ops', 'Ratio');
fprintf('%-20s-+-%-8s-+-%-8s-+-%-8s-+-%-8s-+-%-8s\n', ...
    repmat('-',1,20), repmat('-',1,8), repmat('-',1,8), ...
    repmat('-',1,8), repmat('-',1,8), repmat('-',1,8));

results = [];
for s = 1:size(scenarios, 1)
    label = scenarios{s, 1};
    ff = scenarios{s, 2};
    fb = scenarios{s, 3};
    lambda = scenarios{s, 4};

    % RLS-DFE 复杂度
    N_rls = ff + fb;
    rls_mul = 4 * N_rls^2 + 8 * N_rls;     % 每符号乘法
    rls_add = 4 * N_rls^2 + 4 * N_rls;      % 每符号加法

    % LMS-DFE 复杂度 (2倍抽头)
    N_lms = (ff + fb) * lms_tap_factor;
    lms_mul = 2 * N_lms + 2;
    lms_add = 2 * N_lms + 1;

    ratio = lms_mul / rls_mul * 100;

    fprintf('%-20s | %-8d | %-8d | %-8.0f | %-8.0f | %-6.1f%%\n', ...
        label, N_rls, N_lms, rls_mul, lms_mul, ratio);

    results(s).label = label;
    results(s).N_rls = N_rls;
    results(s).N_lms = N_lms;
    results(s).rls_mul = rls_mul;
    results(s).lms_mul = lms_mul;
    results(s).ratio = ratio;
    results(s).ff = ff;
    results(s).fb = fb;
    results(s).lambda = lambda;
end

%% ===== 预配置 vs 传统自适应接收机的资源对比 =====
% 传统方案: 周期性发送训练序列进行信道估计 + 自适应跟踪
% 预配置方案: 部署前一次估计 → 固定参数运行, 零训练开销
%
% 关键论证: 不是 CPU 周期节省 (那只有 ~2%), 而是:
%   1. 带宽效率: 无训练序列 → 100% 帧用于数据传输
%   2. 能量效率: 无估计计算 + 无导频发射
%   3. 延迟: 无周期性估计中断

fprintf('\n========== 预配置 vs 传统自适应接收机 ==========\n\n');

% --- 带宽效率 ---
% 典型水声通信: 每帧 20% 为导频/训练符号 (保守估计)
pilot_overhead_ratio = 0.20;  % 20% 带宽用于训练
T_frame = 2.0;                % 典型帧长 (秒)
T_train = T_frame * pilot_overhead_ratio;  % 训练时长

fprintf('--- 带宽效率分析 ---\n');
fprintf('传统方案训练开销:   %.0f%% 帧长 (= %.1f s 每 %.1f s 帧)\n', ...
    pilot_overhead_ratio*100, T_train, T_frame);
fprintf('预配置方案训练开销: 0%% (一次部署前估计, 零运行时开销)\n');
fprintf('有效数据吞吐量提升: %.0f%%\n\n', pilot_overhead_ratio*100);

% --- 能量效率 ---
% 传感器节点主要能耗: 发射功率放大 > 接收处理 > DSP计算
% 节省训练序列发射 → 节省发射能量
tx_power_ratio = 0.70;  % 发射占总能耗 70% (水声通信典型值)
energy_saved = pilot_overhead_ratio * tx_power_ratio * 100;
fprintf('--- 能量效率分析 ---\n');
fprintf('发射能耗占比:        %.0f%% (水声传感器典型值)\n', tx_power_ratio*100);
fprintf('节省发射能量:        %.0f%% (无需周期发射训练序列)\n', energy_saved);
fprintf('额外节省:            接收端无需信道估计算法运行\n\n');

% --- 信道估计计算开销 (作为参考) ---
T_update = 2.0;  % 每2秒重新估计
N_fft_ce = 2^nextpow2(2 * 3 * fs - 1);
ce_mul_per_est = 3 * (N_fft_ce / 2) * log2(N_fft_ce);
fprintf('--- 信道估计计算开销 (参考) ---\n');
fprintf('每次信道估计乘法:    %.0f (FFT点数 %d, 3s窗)\n', ce_mul_per_est, N_fft_ce);
fprintf('每秒信道估计乘法:    %.0f (更新间隔 %.1fs)\n', ce_mul_per_est / T_update, T_update);
fprintf('预配置方案:          0 (仅部署时运行一次)\n\n');

%% 绘图1: RLS vs LMS 乘法复杂度柱状图
figure('Position', [100, 100, 600, 400]);
labels = {results.label};
rls_data = [results.rls_mul];
lms_data = [results.lms_mul];

bar_data = [rls_data(:), lms_data(:)];
b = bar(bar_data);
b(1).FaceColor = [0.2 0.4 0.8];
b(2).FaceColor = [0.8 0.4 0.2];
set(gca, 'XTickLabel', labels, 'FontSize', 9);
xtickangle(45);
ylabel('Multiplications per Symbol', 'FontSize', 12);
title('Computational Complexity: RLS-DFE vs LMS-DFE', 'FontSize', 13);
legend('RLS-DFE (O(N^2))', 'LMS-DFE (O(N))', 'Location', 'NorthWest');
grid on;

%% 绘图2: 带宽效率对比
figure('Position', [100, 100, 700, 450]);
subplot(1, 2, 1);
% 饼图: 传统方案的帧结构
pie_data_trad = [1-pilot_overhead_ratio, pilot_overhead_ratio];
pie(pie_data_trad, {'Data (80%)', 'Training (20%)'});
title('Traditional: Per-Frame Structure', 'FontSize', 12);
colormap([0.2 0.6 0.4; 0.9 0.4 0.3]);

subplot(1, 2, 2);
pie_data_pre = [1, 0];
pie(pie_data_pre, {'Data (100%)', ''});
title('Pre-Configured: Per-Frame Structure', 'FontSize', 12);
colormap([0.2 0.6 0.4; 1 1 1]);

sgtitle('Bandwidth Efficiency: Training Overhead Elimination', 'FontSize', 14);

%% 绘图3: 复杂度-信道参数权衡
figure('Position', [100, 100, 600, 400]);
fb_taps = [results.fb];
yyaxis left;
plot(fb_taps, rls_data, 'b-o', 'LineWidth', 2, 'MarkerSize', 10);
ylabel('RLS Multiplications/Symbol', 'Color', 'b');
yyaxis right;
plot(fb_taps, [results.lambda], 'r-s', 'LineWidth', 2, 'MarkerSize', 10);
ylabel('RLS Forgetting Factor \lambda', 'Color', 'r');
xlabel('Feedback Taps (fb)');
title('Complexity vs. Channel Variability Trade-off', 'FontSize', 13);
grid on;

%% 输出汇总
fprintf('========== 部署建议 ==========\n');
fprintf('对于固定深度传感器网络节点 (准静态信道, T_c = 3-4s):\n');
fprintf('  1. 预配置策略消除训练开销:  带宽效率 +%.0f%%, 发射能量 -%.0f%%\n', ...
    pilot_overhead_ratio*100, energy_saved);
fprintf('  2. 计算资源充裕节点 → RLS-DFE (O(N^2), 更优BER)\n');
fprintf('  3. 功耗敏感节点 → LMS-DFE (%.1f%% of RLS 乘法, O(N))\n', ...
    mean([results.ratio]));
fprintf('  4. 824.8m (下温跃层) 需最多反馈抽头, 为计算最密集节点\n\n');

save('complexity_analysis_results.mat', 'results', 'pilot_overhead_ratio', ...
    'energy_saved');
disp('结果已保存到 complexity_analysis_results.mat');
