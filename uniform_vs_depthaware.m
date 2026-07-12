%% =================================================================
%  uniform_vs_depthaware.m — 统一接收机 vs 深度感知接收机对比
%  目的: 量化 "一刀切" 部署策略的性能损失
%  回答审稿人: "查表策略太简单? 统一接收机在不同深度性能如何?"
% =================================================================
cc
addpath('utils');

%% 配置
depths = [475, 624, 824, 924];
depth_labels = {'475.8m (Node A)', '624.8m (Node B)', '824.8m (Node C)', '924.8m (Node D)'};
data_files = {
    '时变冲激响应475_0.1_3.mat',
    '时变冲激响应624_0.1_3.mat',
    '时变冲激响应824_0.1_3.mat',
    '时变冲激响应924_0.1_3.mat'
    };

fs = 16000;
symbol_rate = 4000;
sps = fs / symbol_rate;
M = 4;
T_step = 0.1;
SNR_vec = 0:5:30;

%% 第一步: 为每个深度自动推导最优参数
fprintf('========== 阶段1: 各深度最优参数推导 ==========\n');
all_params = cell(4, 1);
all_h = cell(4, 1);

for d = 1:4
    fprintf('\n--- 深度: %s ---\n', depth_labels{d});
    load(data_files{d});

    % 预处理 (与 Channel_review 一致)
    avg_pdp = mean(abs(h_matrix).^2, 1);
    [~, peak_idx] = max(avg_pdp);
    start_idx = max(1, peak_idx - 100);
    end_idx = min(size(h_matrix,2), peak_idx + 2500);
    h_cut1 = h_matrix(:, start_idx:end_idx);

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

    [h_cut, ~] = extract_multipath_by_threshold(abs(result_matrix), -20, 5);
    h_power = mean(sum(abs(h_cut).^2, 2));
    h_cut = h_cut / sqrt(h_power);
    all_h{d} = h_cut;

    params = auto_parameter_selector(h_cut, fs, symbol_rate, T_step, depth_labels{d});
    all_params{d} = params;
end

%% 第二步: 统一接收机 vs 深度感知接收机 BER 对比
fprintf('\n========== 阶段2: BER 对比 ==========\n');

% 统一接收机: 使用475m的最优参数 (最"乐观"的参数)
uniform_params = all_params{1};
% 深度感知接收机: 每深度用自己的最优参数
depthaware_params = all_params;

% 每个深度跑两组BER: uniform 和 depth-aware
BER_uniform = zeros(4, length(SNR_vec));
BER_depthaware = zeros(4, length(SNR_vec));

for d = 1:4
    fprintf('\n深度 %s:\n', depth_labels{d});
    h_cut = all_h{d};
    time_slices = size(h_cut, 1);
    total_time = time_slices * T_step;
    L_ch = size(h_cut, 2);

    % 发射机
    num_symbols = floor(total_time * symbol_rate);
    rng(42);
    tx_syms_int = randi([0 M-1], num_symbols, 1);
    tx_data = pskmod(tx_syms_int, M, pi/4);
    tx_signal = rectpulse(tx_data, sps);
    len_per_slice = floor(T_step * fs);
    tx_signal = tx_signal(1 : time_slices * len_per_slice);

    % 信道卷积
    rx_signal_clean = zeros(length(tx_signal) + L_ch, 1);
    for k = 1:time_slices
        idx_start = (k-1) * len_per_slice + 1;
        idx_end = k * len_per_slice;
        sig_block = tx_signal(idx_start:idx_end);
        block_out = conv(sig_block, h_cut(k, :).');
        out_idx_end = idx_start + length(block_out) - 1;
        rx_signal_clean(idx_start:out_idx_end) = rx_signal_clean(idx_start:out_idx_end) + block_out;
    end
    rx_signal_clean = rx_signal_clean(1:length(tx_signal));

    % 导频
    sync_len_sym = 2000;
    ref_waveform = tx_signal(1:sync_len_sym * sps);

    for i = 1:length(SNR_vec)
        snr_now = SNR_vec(i);
        rx_noisy = awgn(rx_signal_clean, snr_now, 'measured');

        [xc, lags] = xcorr(rx_noisy, ref_waveform);
        [~, max_idx] = max(abs(xc));
        peak_lag = lags(max_idx - 119256 + 118916);
        if peak_lag < 1, peak_lag = 1; end
        rx_synced = rx_noisy(peak_lag:end);
        rx_sampled = rx_synced(1:sps:end);
        rx_sampled = rx_sampled / mean(abs(rx_sampled));
        min_len = min(length(rx_sampled), length(tx_data));
        rx_sampled = rx_sampled(1:min_len);
        tx_ref_aligned = tx_data(1:min_len);
        train_len = min(1000, floor(length(rx_sampled)/2));

        % --- 统一接收机 (uniform) ---
        p = uniform_params;
        dfe_uniform = comm.DecisionFeedbackEqualizer(...
            'Algorithm', 'RLS', ...
            'NumForwardTaps', p.eq_ff_taps, ...
            'NumFeedbackTaps', p.eq_fb_taps, ...
            'ReferenceTap', p.ref_tap, ...
            'ForgettingFactor', p.rls_forget, ...
            'Constellation', pskmod(0:M-1, M, pi/4));
        [sig_eq_u, ~] = dfe_uniform(rx_sampled, tx_ref_aligned(1:train_len));
        rx_syms_u = pskdemod(sig_eq_u, M, pi/4);
        BER_uniform(d, i) = compute_min_ber(rx_syms_u, tx_syms_int, train_len+50);

        % --- 深度感知接收机 (depth-aware) ---
        p = depthaware_params{d};
        dfe_da = comm.DecisionFeedbackEqualizer(...
            'Algorithm', 'RLS', ...
            'NumForwardTaps', p.eq_ff_taps, ...
            'NumFeedbackTaps', p.eq_fb_taps, ...
            'ReferenceTap', p.ref_tap, ...
            'ForgettingFactor', p.rls_forget, ...
            'Constellation', pskmod(0:M-1, M, pi/4));
        [sig_eq_da, ~] = dfe_da(rx_sampled, tx_ref_aligned(1:train_len));
        rx_syms_da = pskdemod(sig_eq_da, M, pi/4);
        BER_depthaware(d, i) = compute_min_ber(rx_syms_da, tx_syms_int, train_len+50);
    end
    fprintf('  完成 (深度 %d)\n', depths(d));
end

%% 第三步: 绘图 — 统一 vs 深度感知对比
figure('Position', [100, 100, 1200, 500]);
colors = lines(4);
markers = {'o', 's', '^', 'd'};

for d = 1:4
    subplot(1, 4, d);
    semilogy(SNR_vec, BER_uniform(d,:), '-o', 'Color', [0.7 0.2 0.2], ...
        'LineWidth', 1.5, 'MarkerSize', 6); hold on;
    semilogy(SNR_vec, BER_depthaware(d,:), '-s', 'Color', [0.2 0.5 0.2], ...
        'LineWidth', 1.5, 'MarkerSize', 6);
    grid on;
    xlabel('SNR (dB)'); ylabel('BER');
    title(depth_labels{d}, 'FontSize', 10);
    legend('Uniform (475m params)', 'Depth-Aware', 'Location', 'SouthWest', 'FontSize', 8);
    axis([0 30 1e-5 1]);
end
sgtitle('Uniform Receiver vs Depth-Aware Receiver: Per-Node BER Comparison', 'FontSize', 14);

%% 第四步: 汇总图 — 全部深度在一张图上
figure('Position', [100, 100, 800, 550]);

% 深度感知 (实线)
for d = 1:4
    semilogy(SNR_vec, BER_depthaware(d,:), ['-' markers{d}], ...
        'Color', colors(d,:), 'LineWidth', 2, 'MarkerSize', 8); hold on;
end
% 统一接收机 (虚线)
for d = 1:4
    semilogy(SNR_vec, BER_uniform(d,:), ['--' markers{d}], ...
        'Color', colors(d,:) * 0.6, 'LineWidth', 1, 'MarkerSize', 6);
end

grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title('Network-Wide BER: Uniform vs Depth-Aware Receiver Deployment', 'FontSize', 14);
leg_entries = cell(8, 1);
for d = 1:4
    leg_entries{2*d-1} = sprintf('%s (Depth-Aware, fb=%d)', ...
        depth_labels{d}, depthaware_params{d}.eq_fb_taps);
    leg_entries{2*d} = sprintf('%s (Uniform, fb=%d)', ...
        depth_labels{d}, uniform_params.eq_fb_taps);
end
legend(leg_entries, 'Location', 'SouthWest', 'FontSize', 8);
axis([0 30 1e-5 1]);

%% 第五步: 性能恶化量化
fprintf('\n========== 统一接收机性能恶化分析 ==========\n');
fprintf('基准参数 (来自475m): fb=%d, lambda=%.4f\n', ...
    uniform_params.eq_fb_taps, uniform_params.rls_forget);
for d = 2:4
    ber_diff_15dB = BER_uniform(d, 4) / max(1e-6, BER_depthaware(d, 4));
    fprintf('  %s: SNR=15dB 时统一接收机 BER 恶化 %.1f 倍\n', ...
        depth_labels{d}, ber_diff_15dB);
end

%% 保存
save('uniform_vs_depthaware_results.mat', ...
    'BER_uniform', 'BER_depthaware', 'SNR_vec', ...
    'all_params', 'depth_labels', 'depths');
disp('结果已保存到 uniform_vs_depthaware_results.mat');


%% 局部函数
function min_ber = compute_min_ber(rx_syms_all, tx_ref_all, eval_start_idx)
    search_range = -10:10;
    min_ber = 1;
    for d = search_range
        rx_idx_start = eval_start_idx;
        rx_idx_end = length(rx_syms_all);
        tx_idx_start = rx_idx_start + d;
        tx_idx_end = rx_idx_end + d;
        if tx_idx_start < 1 || tx_idx_end > length(tx_ref_all), continue; end
        rx_chunk = rx_syms_all(rx_idx_start:rx_idx_end);
        tx_chunk = tx_ref_all(tx_idx_start:tx_idx_end);
        try
            [~, tmp_ber] = biterr(tx_chunk(1:end-28), rx_chunk(29:end));
        catch, tmp_ber = 1;
        end
        if tmp_ber < min_ber, min_ber = tmp_ber; end
    end
end
