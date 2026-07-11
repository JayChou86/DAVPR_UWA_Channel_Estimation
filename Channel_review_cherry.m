%% =================================================================
%  Project: 水声时变信道重放与通信性能评估 (Channel Replay Simulation)
%  Target: IEEE TVT Trans.
%  Author: Assistant
%  =================================================================

clear; close all; clc;
load('时变冲激响应475_0.1_3.mat')
%% 1. 参数设置 (System Parameters)
fs = 16000;                 % 采样率 16kHz
T_slice = 0.1;              % 信道时间切片分辨率 (s)
[num_slices, max_delay_len] = deal(74, 95999); % 原始信道矩阵维度

% 通信参数
symbol_rate = 2000;         % 符号率 (Baud)，需小于 fs/2，建议 2k-4k
sps = fs / symbol_rate;     % 每个符号的采样点数 (Oversampling factor)
mod_order = 4;              % QPSK
total_time = num_slices * T_slice; % 总时长 7.4s
num_symbols = ceil(total_time * symbol_rate);

% 接收机参数 (DFE)
df_taps_ff = 60;            % 前馈滤波器抽头数 (覆盖主径附近)
df_taps_fb = 120;            % 反馈滤波器抽头数 (消除后向多径)
rls_forget_factor = 0.995;  % RLS 遗忘因子 (针对时变信道)

%% 2. 加载与预处理信道 (Channel Pre-processing)
% =============================================================
% 注意：此处假设 h_matrix 已经加载到工作区。
% 如果没有，请取消下面两行注释生成模拟数据用于测试代码
if ~exist('h_matrix', 'var')
    disp('正在生成模拟信道数据...');
    h_matrix = (randn(74, 95999) + 1j*randn(74, 95999)) * 0.01;
    % 模拟一个移动的主径
    for i = 1:74
        pos = 5000 + i*10; % 模拟多普勒/时延漂移
        h_matrix(i, pos:pos+50) = exp(-0.1*(0:50)); % 主径能量
    end
end

% --- 关键步骤：信道对齐与截断 (Windowing) ---
% 原始 95999 长度包含绝对时延，直接卷积会导致计算量巨大且接收机无法同步。
% 我们需要找到能量中心，截取有效部分。

% 计算平均功率时延谱
pdp = mean(abs(h_matrix).^2, 1);
[~, max_loc] = max(pdp);

% 设定截取窗口 (例如：主径前 50ms, 主径后 200ms)
pre_cursor_s = 0.05; 
post_cursor_s = 0.20;
win_start = max(1, max_loc - round(pre_cursor_s * fs));
win_end = min(max_delay_len, max_loc + round(post_cursor_s * fs));

h_effective1 = h_matrix(:, win_start:win_end); % 截断后的有效信道
[n_slices, n_taps_eff] = size(h_effective1);

h_effective = zeros(size(h_effective1));
imagesc(abs(h_effective1))
h_effective(:,[800,1151]) = abs(h_effective1(:,[800,1151]));
h_effective(:,800) = 1;




% 归一化信道能量 (为了保证 SNR 设置准确)
h_power = mean(sum(abs(h_effective).^2, 2));
h_effective = h_effective / sqrt(h_power);

fprintf('信道预处理完成: 原始长度 %d -> 有效长度 %d (%.2f ms)\n', ...
    max_delay_len, n_taps_eff, n_taps_eff/fs*1000);

figure; imagesc(abs(h_effective)); title('预处理后的有效时变信道 (幅度)');
xlabel('Delay Taps'); ylabel('Time Slice'); colormap jet; colorbar;

%% 3. 发射机 (Transmitter)
% =============================================================
disp('生成发射信号...');
data_bits = randi([0 1], num_symbols * log2(mod_order), 1);

% QPSK 调制
modulator = comm.QPSKModulator('BitInput', true);
tx_symbols = modulator(data_bits);

% 根升余弦成型滤波 (RRC Pulse Shaping)
rrcFilter = comm.RaisedCosineTransmitFilter(...
    'Shape', 'Square root', ...
    'RolloffFactor', 0.25, ...
    'FilterSpanInSymbols', 10, ...
    'OutputSamplesPerSymbol', sps);

tx_waveform = rrcFilter(tx_symbols);
tx_len = length(tx_waveform);

%% 4. 时变信道重放仿真 (Time-Varying Channel Replay)
% =============================================================
% 核心逻辑：在 0.1s 的切片之间进行线性插值，模拟连续变化的水声信道
disp('正在进行时变信道卷积 (这可能需要一点时间)...');

rx_waveform_clean = zeros(tx_len + n_taps_eff, 1);
samples_per_slice = round(T_slice * fs);

% 为了速度，使用重叠相加法 (Overlap-Add) 或者逐样本插值
% 这里为了演示清晰，采用简化的"分段卷积+平滑过渡"逻辑

current_idx = 1;
h_prev = h_effective(1, :).';

for i = 1:n_slices
    % 确定当前时间段的信号索引
    start_idx = (i-1) * samples_per_slice + 1;
    end_idx = min(start_idx + samples_per_slice - 1, tx_len);
    
    if start_idx > tx_len, break; end
    
    sig_segment = tx_waveform(start_idx:end_idx);
    
    % 获取当前时刻的 CIR
    h_curr = h_effective(i, :).';
    
    % --- 简单重放：直接卷积当前切片 (Block Fading assumption within 0.1s) ---
    % 如果想更精细，可以对 h_prev 和 h_curr 做 interpolation
    % 但考虑到计算量，0.1s 分辨率下直接分段卷积是可接受的近似
    
    seg_out = conv(sig_segment, h_curr);
    
    % 叠加到输出信号 (Overlap-Add)
    L_seg = length(seg_out);
    rx_waveform_clean(start_idx : start_idx + L_seg - 1) = ...
        rx_waveform_clean(start_idx : start_idx + L_seg - 1) + seg_out;
end

% 截断多余尾部
rx_waveform_clean = rx_waveform_clean(1:length(tx_waveform));

%% 5. 信道添加噪声 (Add AWGN)
% =============================================================
target_snr_db = 20; % 设定目标信噪比
rx_waveform_noisy = awgn(rx_waveform_clean, target_snr_db, 'measured');

%% 6. 接收机：同步与均衡 (Receiver & Equalization)
% =============================================================
disp('接收机处理 (RLS-DFE)...');

% 1. 简单的粗同步 (基于能量最大值，实际应使用 Preamble)
% 由于我们做的是 replay，且已经截断了信道，信号起始点基本在 0 附近
% 这里做一个简单的相关峰寻找以微调
lag_search_range = 0:500;
corrs = xcorr(rx_waveform_noisy(1:1000), tx_waveform(1:1000), 500);
[~, peak_idx] = max(abs(corrs));
delay_est = peak_idx - 500; 
if delay_est < 0, delay_est = 0; end

rx_sync = rx_waveform_noisy(delay_est+1:end);

% 2. 降采样与匹配滤波 (此处简化，直接抽取，实际应用应过匹配滤波)
% 注意：为了 DFE 工作更好，通常在分数间隔(Fractionally Spaced)下工作
% 这里为了代码简洁，先降采样到符号率
rx_symbols_uneq = rx_sync(1:sps:end);
rx_symbols_uneq = rx_symbols_uneq(1:min(length(rx_symbols_uneq), length(tx_symbols)));
ref_symbols = tx_symbols(1:length(rx_symbols_uneq));

% 3. DFE 均衡 (使用 MATLAB 通信工具箱)
% RLS 算法适合快时变信道
dfe = comm.DecisionFeedbackEqualizer(...
    'Algorithm', 'RLS', ...
    'NumForwardTaps', df_taps_ff, ...
    'NumFeedbackTaps', df_taps_fb, ...
    'ReferenceTap', floor(df_taps_ff/2), ...
    'ForgettingFactor', rls_forget_factor, ...
    'Constellation', complex([-1-1i, -1+1i, 1-1i, 1+1i]/sqrt(2))); 

% 训练模式：前 1000 个符号作为训练序列 (Training Sequence)
train_len = 1000;
[rx_symbols_eq, error_sq] = dfe(rx_symbols_uneq, ref_symbols(1:train_len));

%% 7. 性能评估与可视化 (智能对齐版)
% =============================================================

% --- A. 准备数据 ---
% 取出均衡后的有效数据（去除训练序列）
rx_sym_payload = rx_symbols_eq(train_len+1:end);
% 取出对应的发射参考数据（为了计算方便，取足够长一段）
tx_sym_ref = tx_symbols(train_len+1:end); 

% 确保长度一致，以较短者为准
L_eval = min(length(rx_sym_payload), length(tx_sym_ref));
rx_sym_eval = rx_sym_payload(1:L_eval);
tx_sym_eval = tx_sym_ref(1:L_eval);

% --- B. 解决相位模糊与延时对齐 (关键步骤) ---
% QPSK 有 4 种可能的相位旋转：0, 90, 180, 270 度
phase_rotations = [1, 1j, -1, -1j]; 
min_ber = 1;
best_phase_idx = 0;
best_delay = 0;

% 实例化解调器
demodulator = comm.QPSKDemodulator('BitOutput', true);

fprintf('正在搜索最佳相位与时延对齐...\n');

for p = 1:4
    % 1. 尝试旋转接收符号
    rx_rotated = rx_sym_eval * phase_rotations(p);
    
    % 2. 使用互相关寻找最佳时延 (Symbol level sync)
    % 计算这一小段的互相关，找到对齐点
    % 注意：为了速度，只用前 2000 个符号计算相关性
    check_len = min(2000, L_eval);
    [xc, lags] = xcorr(rx_rotated(1:check_len), tx_sym_eval(1:check_len));
    [~, max_idx] = max(abs(xc));
    delay_lag = lags(max_idx); % rx 相对于 tx 的滞后
    
    % 3. 根据时延对齐数据
    if delay_lag >= 0
        % Rx 比 Tx 滞后 (正常情况，因为有滤波器延迟)
        % Rx: [discard_delay ... data ...]
        % Tx: [data ...]
        rx_aligned = rx_rotated(delay_lag+1:end);
        tx_aligned = tx_sym_eval(1:length(rx_aligned));
    else
        % Rx 比 Tx 超前 (很少见，除非截取时出错)
        rx_aligned = rx_rotated(1:end+delay_lag);
        tx_aligned = tx_sym_eval(1-delay_lag:end);
    end
    
    % 4. 解调并计算 BER
        % 4. 解调并计算 BER
    % -----------------------------------------------------------
    
    % A. 解调接收到的符号 (得到接收比特)
    rx_bits_temp = demodulator(rx_aligned);
    
    % B. 解调发射端的参考符号 (得到标准答案比特)
    % 这里的逻辑是：为了确保映射关系一致，我们把本来就知道的发射符号
    % 也通过解调器过一遍，作为比较的基准 (Ground Truth)
    ref_bits_temp = demodulator(tx_aligned); 
    
    % C. 计算误码率
    [~, current_ber] = biterr(ref_bits_temp, rx_bits_temp);
    
    % -----------------------------------------------------------

    
    % 5. 记录最佳结果
    if current_ber < min_ber
        min_ber = current_ber;
        best_phase_idx = p;
        best_delay = delay_lag;
        final_rx_bits = rx_bits_temp;
        final_ref_bits = ref_bits_temp;
    end
end

% --- C. 输出最终结果 ---
fprintf('\n=== 仿真结果 (优化后) ===\n');
fprintf('SNR: %d dB\n', target_snr_db);
fprintf('最佳相位旋转: %d 度\n', angle(phase_rotations(best_phase_idx))*180/pi);
fprintf('系统处理时延: %d symbols\n', best_delay);
fprintf('最终 BER: %.5f\n', min_ber);

% --- D. 绘图 ---
figure('Position', [100, 100, 1000, 600]);

subplot(2,2,1);
plot(abs(rx_waveform_clean)); title('信道重放信号 (幅度)');
xlabel('Samples'); grid on;

subplot(2,2,2);
plot(10*log10(abs(error_sq))); title('DFE 收敛曲线 (MSE)');
xlabel('Symbol Index'); ylabel('MSE (dB)'); grid on;
ylim([-30 10]);

subplot(2,2,3);
% 画出相位校正后的星座图
rx_final_plot = rx_sym_eval * phase_rotations(best_phase_idx);
plot(real(rx_final_plot(500:end)), imag(rx_final_plot(500:end)), 'b.'); hold on;
plot(real(tx_sym_eval(1:100)), imag(tx_sym_eval(1:100)), 'rx', 'MarkerSize', 10);
legend('接收符号 (校正后)', '标准星座点');
title(sprintf('星座图 (BER=%.4f)', min_ber)); 
axis square; grid on; xlim([-2 2]); ylim([-2 2]);

subplot(2,2,4);
imagesc(abs(h_effective)); 
title('截断后的有效信道冲激响应');
xlabel('Delay Taps'); ylabel('Time Slice');
colormap jet;

figure
plot(ref_bits_temp)
hold on
plot(rx_bits_temp)
xlim([1,50])
ylim([-1,2])