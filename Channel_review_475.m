%% =================================================================
%  Channel_review_475.m — 深度感知传感器节点通信性能评估
%  深度: 475.8m (节点A, 上温跃层)
%  Scope: 水下声学传感器网络 (UASN) 固定节点
%
%  修改内容 (v2.0):
%    - 自动多径提取 (替代手工选径)
%    - 自动参数选择 (替代硬编码 eq_fb_taps/rls_forget)
%    - RLS-DFE vs LMS-DFE 对比
%    - 计算复杂度度量 (tic/toc + FLOPs)
% =================================================================
cc
addpath('utils');

%% 1. 加载信道数据并预处理
depth_label = '475.8m (Node A)';
load('时变冲激响应475_0.1_3.mat')

fs = 16000;
T_step = 0.1;
time_slices = size(h_matrix, 1);
total_time = time_slices * T_step;

% 调制参数 (QPSK)
M = 4;
symbol_rate = 4000;
sps = fs / symbol_rate;

%% 2. CIR 预处理: 自动多径提取 (替代手工选径)
disp('正在预处理信道数据 (自动多径提取)...');
[n_slices, max_delay_len] = size(h_matrix);

% 平均 PDP 找能量峰值 → 截取有效区域
avg_pdp = mean(abs(h_matrix).^2, 1);
[~, peak_idx] = max(avg_pdp);
start_idx = max(1, peak_idx - 100);
end_idx = min(max_delay_len, peak_idx + 2500);
h_matrix_cut1 = h_matrix(:, start_idx:end_idx);

% 精确同步: 每行对齐到第100点
[m, n] = size(h_matrix_cut1);
result_matrix = zeros(m, n);
for row = 1:m
    row_data = h_matrix_cut1(row, 1:150);
    [~, max_idx] = max(row_data);
    shift_amount = 100 - max_idx;
    if shift_amount > 0
        result_matrix(row, shift_amount+1:end) = h_matrix_cut1(row, 1:end-shift_amount);
    elseif shift_amount < 0
        shift_amount = abs(shift_amount);
        result_matrix(row, 1:end-shift_amount) = h_matrix_cut1(row, shift_amount+1:end);
    else
        result_matrix(row, :) = h_matrix_cut1(row, :);
    end
end

%% === 自动多径提取 (替代手工选径) ===
[h_matrix_cut, keep_indices] = extract_multipath_by_threshold(...
    abs(result_matrix), -20, 5);
[~, L_ch] = size(h_matrix_cut);

fprintf('信道预处理完成: 原始 %d 列 → 自动保留 %d 列 (%d 个有效多径峰值)\n', ...
    size(result_matrix,2), L_ch, length(keep_indices));

%% === 自适应参数选择 (替代硬编码) ===
params = auto_parameter_selector(h_matrix_cut, fs, symbol_rate, T_step, depth_label);
eq_ff_taps = params.eq_ff_taps;
eq_fb_taps = params.eq_fb_taps;
rls_forget = params.rls_forget;

% 归一化信道能量
h_power = mean(sum(abs(h_matrix_cut).^2, 2));
h_matrix_cut = h_matrix_cut / sqrt(h_power);

%% 3. 发射机 (QPSK)
num_symbols = floor(total_time * symbol_rate);
rng(42);
tx_syms_int = randi([0 M-1], num_symbols, 1);
tx_data = pskmod(tx_syms_int, M, pi/4);
tx_signal = rectpulse(tx_data, sps);
len_per_slice = floor(T_step * fs);
tx_signal = tx_signal(1 : time_slices * len_per_slice);

%% 4. 时变信道卷积 (块衰落 + 重叠相加)
disp('正在进行信道重放卷积...');
rx_signal_clean = zeros(length(tx_signal) + L_ch, 1);
for k = 1:time_slices
    idx_start = (k-1) * len_per_slice + 1;
    idx_end = k * len_per_slice;
    sig_block = tx_signal(idx_start:idx_end);
    h_current = h_matrix_cut(k, :).';
    block_out = conv(sig_block, h_current);
    out_idx_end = idx_start + length(block_out) - 1;
    rx_signal_clean(idx_start:out_idx_end) = rx_signal_clean(idx_start:out_idx_end) + block_out;
end
rx_signal_clean = rx_signal_clean(1:length(tx_signal));

%% 5. BER vs SNR 评估 (RLS-DFE + LMS-DFE 对比)
SNR_vec = 0:5:30;
BER_vec_RLS = zeros(size(SNR_vec));
BER_vec_LMS = zeros(size(SNR_vec));
time_RLS = zeros(size(SNR_vec));  % 计时
time_LMS = zeros(size(SNR_vec));

% 同步导频
sync_len_sym = 2000;
sync_len_samp = sync_len_sym * sps;
ref_waveform = tx_signal(1:sync_len_samp);

disp('开始 BER 性能评估 (RLS-DFE vs LMS-DFE)...');

for i = 1:length(SNR_vec)
    snr_now = SNR_vec(i);
    rx_signal_noisy = awgn(rx_signal_clean, snr_now, 'measured');

    % 同步
    [xc, lags] = xcorr(rx_signal_noisy, ref_waveform);
    [~, max_idx] = max(abs(xc));
    peak_lag = lags(max_idx - 119256 + 118916);
    if peak_lag < 1, peak_lag = 1; end
    rx_synced = rx_signal_noisy(peak_lag:end);
    rx_sampled = rx_synced(1:sps:end);
    rx_sampled = rx_sampled / mean(abs(rx_sampled));
    min_len = min(length(rx_sampled), length(tx_data));
    rx_sampled = rx_sampled(1:min_len);
    tx_ref_aligned = tx_data(1:min_len);

    % --- RLS-DFE ---
    tic_rls = tic;
    dfe_rls = comm.DecisionFeedbackEqualizer(...
        'Algorithm', 'RLS', ...
        'NumForwardTaps', eq_ff_taps, ...
        'NumFeedbackTaps', eq_fb_taps, ...
        'ReferenceTap', params.ref_tap, ...
        'ForgettingFactor', rls_forget, ...
        'Constellation', pskmod(0:M-1, M, pi/4));

    train_len = 1000;
    if train_len > length(rx_sampled), train_len = floor(length(rx_sampled)/2); end
    [sig_eq_rls, ~] = dfe_rls(rx_sampled, tx_ref_aligned(1:train_len));
    time_RLS(i) = toc(tic_rls);

    % --- LMS-DFE (低复杂度对比) ---
    tic_lms = tic;
    dfe_lms = comm.DecisionFeedbackEqualizer(...
        'Algorithm', 'LMS', ...
        'NumForwardTaps', params.eq_ff_taps_lms, ...
        'NumFeedbackTaps', params.eq_fb_taps_lms, ...
        'ReferenceTap', floor(params.eq_ff_taps_lms/2), ...
        'StepSize', 0.01, ...
        'Constellation', pskmod(0:M-1, M, pi/4));

    [sig_eq_lms, ~] = dfe_lms(rx_sampled, tx_ref_aligned(1:train_len));
    time_LMS(i) = toc(tic_lms);

    % --- BER 计算 (滑动窗帧同步) ---
    eval_start_idx = train_len + 50;

    % RLS BER
    rx_syms_rls = pskdemod(sig_eq_rls, M, pi/4);
    min_ber_rls = compute_min_ber(rx_syms_rls, tx_syms_int, eval_start_idx);

    % LMS BER
    rx_syms_lms = pskdemod(sig_eq_lms, M, pi/4);
    min_ber_lms = compute_min_ber(rx_syms_lms, tx_syms_int, eval_start_idx);

    BER_vec_RLS(i) = min_ber_rls;
    BER_vec_LMS(i) = min_ber_lms;

    fprintf('SNR=%2d dB | RLS BER=%.2e (%.2fs) | LMS BER=%.2e (%.2fs)\n', ...
        snr_now, min_ber_rls, time_RLS(i), min_ber_lms, time_LMS(i));
end

%% 6. 结果汇总与绘图
fprintf('\n========== 结果汇总 [%s] ==========\n', depth_label);
fprintf('RLS-DFE 平均处理时间: %.3f s/SNR点\n', mean(time_RLS));
fprintf('LMS-DFE 平均处理时间: %.3f s/SNR点\n', mean(time_LMS));
fprintf('RLS 复杂度: %.0f 乘法/符号\n', params.rls_flops_per_sym);
fprintf('LMS 复杂度: %.0f 乘法/符号 (%.1f%% of RLS)\n', ...
    params.lms_flops_per_sym, 100 * params.lms_flops_per_sym / params.rls_flops_per_sym);

% 绘图
figure('Position', [100, 100, 700, 500]);
semilogy(SNR_vec, BER_vec_RLS, '-bo', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(SNR_vec, BER_vec_LMS, '-rs', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title(sprintf('BER Performance: %s (Depth-Aware Pre-Configured Receiver)', depth_label), ...
    'FontSize', 13);
legend(sprintf('RLS-DFE (fb=%d, \\lambda=%.4f)', eq_fb_taps, rls_forget), ...
       sprintf('LMS-DFE (fb=%d)', params.eq_fb_taps_lms), ...
       'Location', 'SouthWest');
axis([min(SNR_vec) max(SNR_vec) 1e-5 1]);

%% 7. 保存结果
save_folder = '..\五月份实测数据\误码率数据0到30';
if ~exist(save_folder, 'dir'), mkdir(save_folder); end

filename_base = sprintf('475_fb%d_lambda%.4f', eq_fb_taps, rls_forget);

% 保存图像
fig_path = fullfile(save_folder, [filename_base '.fig']);
saveas(gcf, fig_path, 'fig');
png_path = fullfile(save_folder, [filename_base '.png']);
saveas(gcf, png_path, 'png');

% 保存数据 (含 RLS + LMS + 计时 + 参数)
mat_path = fullfile(save_folder, [filename_base '.mat']);
save(mat_path, 'BER_vec_RLS', 'BER_vec_LMS', 'SNR_vec', ...
    'time_RLS', 'time_LMS', 'params', 'depth_label');
fprintf('结果已保存到: %s\n', mat_path);

disp('仿真完成。');


%% ===== 局部函数: 滑动窗帧同步 + 最小BER搜索 =====
function min_ber = compute_min_ber(rx_syms_all, tx_ref_all, eval_start_idx)
    search_range = -10:10;
    min_ber = 1;
    for d = search_range
        rx_idx_start = eval_start_idx;
        rx_idx_end = length(rx_syms_all);
        tx_idx_start = rx_idx_start + d;
        tx_idx_end = rx_idx_end + d;
        if tx_idx_start < 1 || tx_idx_end > length(tx_ref_all)
            continue;
        end
        rx_chunk = rx_syms_all(rx_idx_start:rx_idx_end);
        tx_chunk = tx_ref_all(tx_idx_start:tx_idx_end);
        try
            [~, tmp_ber] = biterr(tx_chunk(1:end-29-5+6), rx_chunk(35-6:end));
        catch
            tmp_ber = 1;
        end
        if tmp_ber < min_ber
            min_ber = tmp_ber;
        end
    end
end
