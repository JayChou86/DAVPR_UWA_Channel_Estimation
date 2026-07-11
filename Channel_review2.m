%% =================================================================
%  Project: 针对“长时延拓展”信道的 TR-DFE (重型版)
%  Target: 解决能量弥散导致的星座图发散
%  代码不好，放弃了，AIstdio后半部分生成的1.26
%  =================================================================
clc; clear variables; close all;

%% 1. 加载/模拟数据
if ~exist('h_matrix', 'var')
    % 模拟长时延拓展信道 (两个宽山包)
    h_static = zeros(1, 60000); 
    % 第一径：47000处，弥散 100 点
    idx1 = 47000:47100;
    h_static(idx1) = (0.3+0.3j) * hann(length(idx1))'; 
    % 第二径：47400处，弥散 150 点 (更强、更宽)
    idx2 = 47400:47550;
    h_static(idx2) = (1.0-0.5j) * hann(length(idx2))';
    disp('【模拟】生成带“长时延拓展”的复数信道...');
else
    h_static = h_matrix(40, :); 
    disp('【实测】使用 h_matrix 第 40 行...');
end
% 归一化
h_static = h_static / norm(h_static);

%% 2. 发射机
M = 4; symbol_rate = 4000; fs = 16000; sps = 4;
num_symbols = 5000; 

tx_syms_int = randi([0 M-1], num_symbols, 1);
tx_data = pskmod(tx_syms_int, M, pi/4);
tx_signal = rectpulse(tx_data, sps);

%% 3. 过信道
rx_signal = conv(tx_signal, h_static.'); 
rx_signal = awgn(rx_signal, 20, 'measured');

%% 4. TR 处理 (Conjugate + Flip)
disp('执行 TR 处理...');
rx_tr = conj(flipud(rx_signal)); 
tx_ref_complex = conj(flipud(tx_data));

% --- 【优化】能量重心同步 (Energy Search) ---
% 面对弥散信道，max(abs) 不可靠，使用滑动能量窗
window_len = sps * 10; % 积分窗
energy_profile = movsum(abs(rx_tr).^2, window_len);
[~, center_idx] = max(energy_profile);

% 粗略定位到能量团的开始位置
start_idx = max(1, center_idx - 20 * sps); 

rx_sync = rx_tr(start_idx : end);
rx_samp = downsample(rx_sync, sps);
rx_samp = rx_samp / mean(abs(rx_samp));

% 对齐长度
L_proc = min(length(rx_samp), length(tx_ref_complex));
rx_proc = rx_samp(1:L_proc);
tx_proc = tx_ref_complex(1:L_proc);

%% 5. TR-DFE 均衡 (参数扩容)
% 针对“长时延拓展”的关键修改：
% 1. 前馈抽头加倍：收集反转后甩到前面的弥散能量
% 2. 参考抽头后移：让EQ有时间处理
eq_ff_taps = 150;    % 原来50 -> 现在100 (覆盖前驱弥散)
eq_fb_taps = 150;    % 覆盖两条径之间的距离
ref_tap    = 80;     % 放在前馈的后部 (关键!)
rls_forget = 0.99;  % 稍微降低遗忘因子，适应弥散导致的快变相位

dfe = comm.DecisionFeedbackEqualizer(...
    'Algorithm', 'RLS', ...
    'NumForwardTaps', eq_ff_taps, ...
    'NumFeedbackTaps', eq_fb_taps, ...
    'ReferenceTap', ref_tap, ...
    'ForgettingFactor', rls_forget, ...
    'Constellation', pskmod(0:M-1, M, pi/4));

% 训练
train_len = 3000; % 增加训练长度
[sig_eq, ~] = dfe(rx_proc, tx_proc(1:train_len));

%% 6. 结果评估
if length(sig_eq) > train_len + 200
    steady_sig = sig_eq(train_len+100:end);
    
    % --- 真值生成 ---
    tx_ref_truth_int = pskdemod(tx_proc(train_len+100:end), M, pi/4);
    rx_demod = pskdemod(steady_sig, M, pi/4);
    
    % --- 强力同步搜索 ---
    probe_len = min(800, length(rx_demod)); % 增加探针长度
    probe_seq = steady_sig(1:probe_len); 
    ref_seq = tx_proc(train_len+100:end); 
    
    [xc, lags] = xcorr(ref_seq, probe_seq);
    [max_val, max_idx] = max(abs(xc));
    best_lag = lags(max_idx);
    
    % 绘图
    figure('Position',[100 100 1000 400]);
    subplot(1,3,1); 
    plot(abs(h_static)); title('信道 (注意宽度)'); xlim([46000 48000]);
    subplot(1,3,2); 
    plot(abs(sig_eq)); title('均衡器输出幅度 (检查收敛)');
    xlabel('Symbols'); ylabel('Amp');
    subplot(1,3,3); 
    scatterplot(steady_sig); title(['星座图 (Lag=' num2str(best_lag) ')']);
    
    % 计算 BER
    rx_start_in_ref = best_lag + 1;
    if rx_start_in_ref > 0 && rx_start_in_ref + length(rx_demod) - 1 <= length(tx_ref_truth_int)
        truth_chunk = tx_ref_truth_int(rx_start_in_ref : rx_start_in_ref + length(rx_demod) - 1);
        [~, final_ber] = biterr(truth_chunk, rx_demod);
    else
        L_common = min(length(rx_demod), length(tx_ref_truth_int));
        [~, final_ber] = biterr(tx_ref_truth_int(1:L_common), rx_demod(1:L_common));
    end
    
    fprintf('\n========================================\n');
    fprintf('1. 均衡器收敛状态:    %s\n', mat2str(mean(abs(steady_sig(end-100:end))) > 0.5));
    fprintf('2. 复数相关峰值:      %.2f (目标 > 400)\n', max_val);
    fprintf('3. 最终 BER:          %.6f\n', final_ber);
    fprintf('========================================\n');
    
    if max_val < 100
        disp('❌ 失败提示：相关峰值太低，说明星座图依然发散。');
        disp('   建议：继续增大 eq_ff_taps (如 120) 或 减小 rls_forget (如 0.99)');
    end
else
    disp('数据过短');
end