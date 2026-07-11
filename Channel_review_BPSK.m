%% ==============================================================================
%  Deep Sea Acoustic Channel Replay Simulation for TVT
%  Author: Assistant
%  Date: 2026-01-26
%  Description: Evaluates BER performance using measured time-varying CIR.
% ==============================================================================

clc; close all;
load('时变冲激响应924_0.1_3.mat')

%% 1. 参数设置 (System Parameters)
fs = 16000;              % 采样率 (Hz)
fc = 800;                % 中心频率 (Hz) - 仅用于说明，基带仿真中体现为带宽
sym_rate = 1000;         % 符号率 (Symbols/s) - 建议根据实际带宽调整 (如 B=fc/2)
sps = fs / sym_rate;     % 每个符号的采样点数 (Samples per Symbol)
mod_order = 2;           % BPSK Modulation
slice_time = 0.1;        % 信道切片时间长度 (s)
slice_len = round(slice_time * fs); % 每个切片的采样点数 (1600点)

% 检查输入数据
if ~exist('h_matrix', 'var')
    error('请先加载你的信道数据 h_matrix (74x95999)');
end
[num_slices, max_delay_taps] = size(h_matrix);

%% 2. 信道数据预处理 (Channel Pre-processing) [修正版]
fprintf('正在预处理信道数据...\n');
% 计算平均功率延迟谱
avg_pdp = mean(abs(h_matrix).^2, 1);

% --- 修正：使用“以峰值为中心”的固定窗口截取 ---
[~, max_idx] = max(avg_pdp); % 找到能量最强的点

% 设定保留范围：峰值前 50ms (800点) + 峰值后 200ms (3200点)
% 这个范围 (250ms) 足够覆盖深海多径，同时去除了前面几秒的空白时延
pre_cursor = 800;  
post_cursor = 3200; 

start_idx = max(1, max_idx - pre_cursor);
end_idx = min(max_delay_taps, max_idx + post_cursor);

% 截取并归一化
h_eff = h_matrix(:, start_idx:end_idx); 
[~, L_ch] = size(h_eff);
avg_power = mean(sum(abs(h_eff).^2, 2));
h_eff = h_eff / sqrt(avg_power);

fprintf('信道已优化截断: 峰值位置 %d, 保留长度 %d (%.1f ms)\n', ...
    max_idx, L_ch, L_ch/fs*1000);

%% 3. 发射机 (Transmitter)
% 计算总仿真时间
total_time = num_slices * slice_time;
num_syms = floor(total_time * sym_rate) - 200; % 预留一点尾部空间

% 生成随机比特
tx_bits = randi([0 1], num_syms, 1);

% BPSK 调制
tx_syms = pskmod(tx_bits, mod_order);

% 脉冲成型 (RRC Filter)
rolloff = 0.5; span = 10;
rrcFilter = rcosdesign(rolloff, span, sps);
tx_signal = upfirdn(tx_syms, rrcFilter, sps);

% 补零以匹配分段长度
total_len_needed = num_slices * slice_len;
tx_signal_padded = [tx_signal; zeros(total_len_needed - length(tx_signal), 1)];

%% 4. 时变信道重放 (Time-Varying Channel Replay)
% 核心步骤：分段卷积 + 重叠相加 (Overlap-Add)

fprintf('正在执行时变信道卷积...\n');
rx_signal_clean = zeros(length(tx_signal_padded) + L_ch, 1);

for i = 1:num_slices
    % 1. 提取当前时间窗的发射信号
    idx_start = (i-1)*slice_len + 1;
    idx_end = i*slice_len;
    
    if idx_end > length(tx_signal_padded), break; end
    sig_segment = tx_signal_padded(idx_start:idx_end);
    
    % 2. 提取当前时刻的 CIR
    h_current = h_eff(i, :).'; % 转置为列向量
    
    % 3. 卷积 (Convolution)
    rx_segment = conv(sig_segment, h_current);
    
    % 4. 重叠相加放入接收buffer
    rx_idx_end = idx_start + length(rx_segment) - 1;
    rx_signal_clean(idx_start:rx_idx_end) = ...
        rx_signal_clean(idx_start:rx_idx_end) + rx_segment;
end

% 截掉多余的尾部，只保留对应发射长度
rx_signal_clean = rx_signal_clean(1:length(tx_signal_padded));

%% 5. 接收机与性能评估 (Rx & Evaluation) [同步逻辑修复版]

SNR_vec = 0:3:21; 
BER_vec = zeros(size(SNR_vec));

% 均衡器参数
nFwd = 60;  
nFb = 30;   
rls_forget = 0.995; 
filter_delay = span; 

fprintf('正在执行全局同步校准 (基于复数序列)...\n');

% --- 步骤 A: 全局同步 ---
% 1. 匹配滤波 + 下采样
rx_clean_syms_all = upfirdn(rx_signal_clean, rrcFilter, 1, sps);

% 2. 去头
if length(rx_clean_syms_all) > filter_delay
    rx_clean_syms = rx_clean_syms_all(filter_delay+1 : end);
else
    rx_clean_syms = rx_clean_syms_all;
end

% 3. 互相关同步 (关键修改：不要对信号取 abs !)
L_corr = min(length(rx_clean_syms), length(tx_syms));
% 使用原始复数信号进行互相关，寻找波形相似度
[xc, lags] = xcorr(rx_clean_syms(1:L_corr), tx_syms(1:L_corr));

% 对“互相关结果”取模，找到能量峰值
[max_val, max_idx] = max(abs(xc));
best_lag = lags(max_idx);

% 计算初始相位偏差 (Phase Offset)
peak_complex_val = xc(max_idx);
initial_phase_rot = angle(peak_complex_val); 

fprintf('同步完成: lag = %d, 峰值强度 = %.2e, 初始相位偏差 = %.2f rad\n', ...
    best_lag, max_val, initial_phase_rot);

% --- 步骤 B: SNR 循环 ---
fprintf('开始 BER 性能扫描...\n');

for k = 1:length(SNR_vec)
    snr_db = SNR_vec(k);
    
    % 1. 加噪
    rx_signal_noisy = awgn(rx_signal_clean, snr_db, 'measured');
    
    % 2. 匹配滤波
    rx_syms_all = upfirdn(rx_signal_noisy, rrcFilter, 1, sps);
    
    % 3. 去滤波器延迟
    if length(rx_syms_all) > filter_delay
        rx_raw_syms = rx_syms_all(filter_delay+1 : end);
    else
        rx_raw_syms = rx_syms_all;
    end
    
    % 4. 对齐数据
    if best_lag >= 0
        s_idx = best_lag + 1;
        valid_len = min(length(rx_raw_syms)-s_idx+1, length(tx_syms));
        rx_input = rx_raw_syms(s_idx : s_idx+valid_len-1);
        ref_sig = tx_syms(1 : valid_len);
    else
        s_idx = -best_lag + 1;
        valid_len = min(length(rx_raw_syms), length(tx_syms)-s_idx+1);
        rx_input = rx_raw_syms(1 : valid_len);
        ref_sig = tx_syms(s_idx : s_idx+valid_len-1);
    end
    
    % 再次截断确保长度一致
    L_proc = min(length(rx_input), length(ref_sig));
    rx_input = rx_input(1:L_proc);
    ref_sig = ref_sig(1:L_proc);
    
    % --- 关键修改：功率归一化 (Power Normalization) ---
    % RLS 对幅度敏感，将接收信号归一化到单位功率
    rx_power = mean(abs(rx_input).^2);
    if rx_power > 0
        rx_input = rx_input / sqrt(rx_power);
    end
    
    % --- 关键修改：相位粗补偿 (Phase Derotation) ---
    % 利用同步时计算的相位偏差，先转回来一点，减轻均衡器压力
    rx_input = rx_input * exp(-1j * initial_phase_rot);
    
    % 5. 均衡器
    dfeObj = comm.DecisionFeedbackEqualizer(...
        'Algorithm', 'RLS', ...
        'NumForwardTaps', nFwd, ...
        'NumFeedbackTaps', nFb, ...
        'ReferenceTap', floor(nFwd/2), ...
        'ForgettingFactor', rls_forget, ...
        'Constellation', complex([-1 1])); 
    
    [rx_syms_eq, ~] = dfeObj(rx_input, ref_sig);
    
    % 6. BER 计算
    converge_idx = 500; 
    
    if length(rx_syms_eq) > converge_idx + 100
        rx_bits_dec = pskdemod(rx_syms_eq(converge_idx:end), mod_order);
        ref_bits = pskdemod(ref_sig(converge_idx:end), mod_order);
        
        [numErrors, ber] = biterr(ref_bits, rx_bits_dec);
        BER_vec(k) = ber;
    else
        numErrors = -1;
        BER_vec(k) = 0.5;
    end
    
    fprintf('SNR: %2d dB | BER: %.2e | Errors: %d\n', snr_db, BER_vec(k), numErrors);
end

%% 6. 绘图
if ishandle(1), close(1); end
figure('Color', 'w');

%% 6. 绘图 (Plotting)
if ishandle(1), close(1); end
figure('Color', 'w');
semilogy(SNR_vec, BER_vec, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title(sprintf('Replay Simulation (Lag: %d, Symbols: %d)', best_lag, length(rx_input)), 'FontSize', 12);
axis([min(SNR_vec) max(SNR_vec) 1e-5 1]);

hold on;
semilogy(SNR_vec, berawgn(SNR_vec,'psk',2,'nondiff'), 'k--', 'LineWidth', 1.5);
legend('RLS-DFE Performance', 'Theory (AWGN)', 'Location', 'SouthWest');
%% 6. 绘图 (Plotting)
if ishandle(1), close(1); end
figure('Color', 'w');
semilogy(SNR_vec, BER_vec, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title(sprintf('Performance at Selected Depth (Lag: %d)', best_lag), 'FontSize', 14);
axis([min(SNR_vec) max(SNR_vec) 1e-5 1]);

% 理论参考
hold on;
semilogy(SNR_vec, berawgn(SNR_vec,'psk',2,'nondiff'), 'k--', 'LineWidth', 1.5);
legend('Simulated (RLS-DFE)', 'Theoretical BPSK', 'Location', 'SouthWest');

%% 6. 绘图 (Plotting for Paper)
figure('Color', 'w');
semilogy(SNR_vec, BER_vec, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title(['BPSK Performance over Replayed Channel (Depth: Unknown)'], 'FontSize', 14);
axis([min(SNR_vec) max(SNR_vec) 1e-5 1]);

% 添加理论 BPSK 曲线对比
hold on;
ber_theo = berawgn(SNR_vec, 'psk', 2, 'nondiff');
semilogy(SNR_vec, ber_theo, 'k--', 'LineWidth', 1.5);
legend('Replayed Channel (RLS-DFE)', 'AWGN Theory', 'Location', 'SouthWest');

fprintf('仿真完成。\n');