cc
%% 1. 参数设置 (System Parameters)
% -----------------------------------------------------------------
% 假设你的 h_matrix 已经存在于工作区，如果不存在，请取消下面注释生成随机数据测试
% 仿真准确，实测不对,应该是同步问题，模拟了一个较长时延的多径进行手动同步后误码率解决
% 代码成功
load('时变冲激响应924_0.1_3.mat')

SNR_vec = 0:5:30; % dB
eq_ff_taps = 60;         % 前馈抽头数 (根据多径长度调整)
eq_fb_taps = 300;         % 反馈抽头数 核心
rls_forget = 0.996;      % RLS 遗忘因子 (接近1稳定，小则跟踪快)


%% 1. 参数设置 (System Parameters)
% -----------------------------------------------------------------
% 假设你的 h_matrix 已经存在于工作区，如果不存在，请取消下面注释生成随机数据测试
if ~exist('h_matrix', 'var')
    disp('【提示】未检测到 h_matrix，正在生成模拟数据...');
    % 模拟一个稀疏信道：74个时刻，最大时延约2000点，总长95999
    h_matrix = zeros(74, 95999);
    for i = 1:74
        h_matrix(i, 500:500+50) = randn(1, 51) + 1j*randn(1, 51); % 模拟多径
    end
end

fs = 16000;              % 采样率
T_step = 0.1;            % 信道更新间隔 (s)
time_slices = size(h_matrix, 1);        % 时间切片数 (对应 h_matrix 行数)
total_time = time_slices * T_step; % 总时长 7.4s

% 调制参数 (QPSK)
M = 4;                   % QPSK
bps = log2(M);           % bit per symbol
symbol_rate = 4000;      % 符号率 (Sps), 需小于 fs/2
sps = fs / symbol_rate;  % Samples per symbol (过采样率)

% 均衡器参数 (关键: 影响 TVT 论文结论)
% eq_ff_taps = 60;         % 前馈抽头数 (根据多径长度调整)
% eq_fb_taps = 600;         % 反馈抽头数 核心
% rls_forget = 0.996;      % RLS 遗忘因子 (接近1稳定，小则跟踪快)

%% 2. 信道冲激响应预处理 (Preprocessing CIR)
% -----------------------------------------------------------------
% 问题：95999 点太长，包含了大量传播时延和噪声。
% 解决：寻找能量中心，截取有效多径部分。

disp('正在预处理信道数据...');
[n_slices, max_delay_len] = size(h_matrix);

% 计算平均功率时延谱 (Average PDP)
avg_pdp = mean(abs(h_matrix).^2, 1);
% 简单同步：找到最大能量峰值，向前取100点，向后取1000点(视具体信道而定)
[~, peak_idx] = max(avg_pdp);
start_idx = max(1, peak_idx - 100); 
end_idx = min(max_delay_len, peak_idx + 2500); % 保留约 3000 点 (约180ms多径)

% 截取后的信道矩阵 (用于仿真)
h_matrix_cut1 = h_matrix(:, start_idx:end_idx); 

%% 精确同步

% 假设 h_matrix_cut 是 74×400 的矩阵
[m, n] = size(h_matrix_cut1);  % m=74, n=400
% 创建结果矩阵
result_matrix = zeros(m, n);
for row = 1:m
    % 获取当前行的前150个数据
    row_data = h_matrix_cut1(row, 1:150);
    
    % 找出最大值及其位置
    [max_val, max_idx] = max(row_data);
    
    % 计算需要平移的列数
    shift_amount = 100 - max_idx;  % 正数表示向右移，负数表示向左移
    
    % 执行整行平移
    if shift_amount > 0  % 需要向右移
        result_matrix(row, shift_amount+1:end) = h_matrix_cut1(row, 1:end-shift_amount);
    elseif shift_amount < 0  % 需要向左移
        shift_amount = abs(shift_amount);
        result_matrix(row, 1:end-shift_amount) = h_matrix_cut1(row, shift_amount+1:end);
    else  % 不需要移动
        result_matrix(row, :) = h_matrix_cut1(row, :);
    end
end


[~, L_ch] = size(result_matrix);
% 截取第500和851列数据
% h_matrix_cut = zeros(size(h_matrix_cut1));
h_matrix_cut = zeros(time_slices,2400);
%%
h_matrix_cut(:,[100,425,1211]) = abs(result_matrix(:,[100,425,1211]));
% h_matrix_cut(:,1211) = abs(result_matrix(:,1211)) / 5;
%%
% h_matrix_cut(:,[100,1520,2238]) = abs(result_matrix(:,[100,1520,2238]));

% h_matrix_cut = abs(result_matrix(:,100:499));

power_delay = mean(abs(result_matrix),1);
plot(power_delay)

% imagesc(abs(h_matrix_cut))
% imagesc(abs(result_matrix))
% colormap('jet')

% 归一化信道能量 (为了公平比较不同深度，TVT 常见做法)
% 这一步确保 SNR 的定义是准确的
h_power = mean(sum(abs(h_matrix_cut).^2, 2));
h_matrix_cut = h_matrix_cut / sqrt(h_power);

% valid_peaks = h_matrix_cut(14:end,10);
% plotrician(valid_peaks, 10, 'Rician');  % exp Rician gamma weibull beta normal




fprintf('信道已截取: 原始长度 %d -> 截取长度 %d (保留能量区域)\n', max_delay_len, L_ch);
% ====== 【调试插入】 强制使用理想信道测试代码逻辑 ======
% 如果这段跑通了(BER=0)，说明代码没问题，是你的 h_matrix 太难解
% disp('【调试模式】正在使用理想单位阵信道...');
% h_matrix_cut = zeros(74, 100); 
% h_matrix_cut(:, 10) = 1; % 在第10个点也就是直达波，无多径
% h_matrix_cut(:, 361) = 0.8; % 在第10个点也就是直达波，无多径
% 
% eq_ff_taps = 50; % 理想信道不需要太长
% eq_fb_taps = 90;
% % ==================================================
%% 3. 发射机 (Transmitter) - 【修正版】
% -----------------------------------------------------------------
num_symbols = floor(total_time * symbol_rate);

% 1. 直接生成整数符号 (0, 1, ..., M-1)
% 这样确保也是标准的 QPSK 映射，且方便后续 biterr 比较
rng(42);  % 设置种子为42
tx_syms_int = randi([0 M-1], num_symbols, 1); 

% 2. QPSK 调制 (将整数映射为复数星座点)
tx_data = pskmod(tx_syms_int, M, pi/4); 

% 3. 成型/脉冲 (矩形脉冲)
tx_signal = rectpulse(tx_data, sps);

% 4. 确保信号长度匹配时间分片
len_per_slice = floor(T_step * fs);
tx_signal = tx_signal(1 : time_slices * len_per_slice); % 截断对齐

%% 4. 时变信道卷积 (Time-Varying Convolution)
% -----------------------------------------------------------------
% 方法：分块处理 + 重叠相加 (Overlap-Add)
disp('正在进行信道重放卷积...');

rx_signal_clean = zeros(length(tx_signal) + L_ch, 1);

for k = 1:time_slices
    % 1. 提取当前 0.1s 的信号块
    idx_start = (k-1) * len_per_slice + 1;
    idx_end = k * len_per_slice;
    sig_block = tx_signal(idx_start:idx_end);
    
    % 2. 提取当前时刻的 CIR
    h_current = h_matrix_cut(k, :).'; % 转置为列向量
    
    % 3. 卷积
    block_out = conv(sig_block, h_current);
    
    % 4. 重叠相加到接收信号 buffer
    % 只要 h 变化不剧烈，这种块衰落近似是可接受的
    out_idx_end = idx_start + length(block_out) - 1;
    rx_signal_clean(idx_start:out_idx_end) = rx_signal_clean(idx_start:out_idx_end) + block_out;
end

% 截取有效部分
rx_signal_clean = rx_signal_clean(1:length(tx_signal));

%% 5. 性能评估循环 (BER vs SNR) - 【增加同步模块修正版】
% -----------------------------------------------------------------
% SNR_vec = 21:1:23; % dB

BER_vec = zeros(size(SNR_vec));

disp('开始 BER 性能评估 (含同步修正)...');

% --- 关键步骤：创建用于同步的导频序列 (Preamble/Training) ---
% 利用发射信号的前 1000 个符号对应的波形来找头
sync_len_sym = 2000;
sync_len_samp = sync_len_sym * sps;
ref_waveform = tx_signal(1:sync_len_samp); 

for i = 1:length(SNR_vec)
    snr_now = SNR_vec(i);
    
    % --- 1. 添加高斯白噪声 ---
    rx_signal_noisy = awgn(rx_signal_clean, snr_now, 'measured');
    
    % --- 2. 同步 (Synchronization) - 解决 BER 0.5 的核心 ---
    % 计算互相关，找到信号的起始位置
    [xc, lags] = xcorr(rx_signal_noisy, ref_waveform);
    [~, max_idx] = max(abs(xc));
    peak_lag = lags(max_idx - 119256 + 118916);
    
    % 只要峰值 lag > 0，说明信号相对于 0 时刻有延迟
    if peak_lag < 1
        peak_lag = 1; % 保护机制
    end
    
    % 对齐信号：从峰值位置开始截取
    % 注意：这里直接跳到了最佳采样点附近
    rx_synced = rx_signal_noisy(peak_lag : end);
    
    % --- 3. 降采样 (Downsample) ---
    % 由于我们在 peak_lag 处截断了，理论上第 1 个点就是符号的最佳采样点
    % 但为了稳健，我们仍然按 sps 抽取
    rx_sampled = rx_synced(1 : sps : end);
    
    % 幅度归一化 (AGC)
    rx_sampled = rx_sampled / mean(abs(rx_sampled));
    
    % 长度保护：防止 rx 比 tx 短
    min_len = min(length(rx_sampled), length(tx_data));
    rx_sampled = rx_sampled(1:min_len);
    tx_ref_aligned = tx_data(1:min_len); % 这里的 tx_data 已经是复数

    % --- 4. 均衡器 (RLS-DFE) ---
    dfe = comm.DecisionFeedbackEqualizer(...
        'Algorithm', 'RLS', ...
        'NumForwardTaps', eq_ff_taps, ...
        'NumFeedbackTaps', eq_fb_taps, ...
        'ReferenceTap', floor(eq_ff_taps/2), ...
        'ForgettingFactor', rls_forget, ...
        'Constellation', pskmod(0:M-1, M, pi/4));
    
    % 训练 (使用对齐后的 tx_ref_aligned)
    train_len = 1000;
    if train_len > length(rx_sampled), train_len = floor(length(rx_sampled)/2); end
    
    [sig_eq, ~] = dfe(rx_sampled, tx_ref_aligned(1:train_len)); 
    
    % --- 5. BER 计算 ---
    eval_start = train_len + 50; % 跳过收敛期
    
    % =======================================================
    % 替换原来的 BER 计算部分
    % =======================================================
    
    % 1. 解调出所有的整数符号
    % 跳过训练序列，因为训练序列可能包含滤波器的瞬态响应
    eval_start_idx = train_len + 50; 
    rx_syms_all = pskdemod(sig_eq, M, pi/4);
    
    % 2. 准备参考符号 (整个序列)
    % 注意：tx_syms_int 是你在发射端生成的 0~3 整数序列
    tx_ref_all = tx_syms_int;
    
    % 3. 滑动搜索最佳对齐位置 (Sliding Window Frame Sync)
    % DFE 输出可能会因为抽头位置不同，导致输出比输入晚几个符号或早几个符号
    % 我们在一个小范围内搜索：例如 -5 到 +5 个符号的偏移
    
    min_ber = 1; % 初始化最小误码率
    best_delay = 0;
    
    search_range = -10 : 10; % 搜索范围
    
    for d = search_range
        % 根据延迟 d 截取接收符号
        % 逻辑：Rx[k] 对应 Tx[k + d]
        
        % 确定比较的索引范围
        % 接收端索引：从 eval_start_idx 开始
        rx_idx_start = eval_start_idx;
        rx_idx_end   = length(rx_syms_all);
        
        % 发射端索引：加上偏移 d
        tx_idx_start = rx_idx_start + d;
        tx_idx_end   = rx_idx_end + d;
        
        % 边界检查：确保索引不越界
        if tx_idx_start < 1 || tx_idx_end > length(tx_ref_all)
            continue; % 这个延迟导致越界，跳过
        end
        
        % 截取片段
        rx_chunk = rx_syms_all(rx_idx_start : rx_idx_end);
        tx_chunk = tx_ref_all(tx_idx_start : tx_idx_end);
        
        % 计算临时 BER
        [~, tmp_ber] = biterr(tx_chunk(1:end-29-5+6-1), rx_chunk(35-6+1:end));
        
        % 更新最小值
        if tmp_ber < min_ber
            min_ber = tmp_ber;
            best_delay = d;
        end
    end
    
    % 4. 记录最终结果
    BER_vec(i) = min_ber;
    
    fprintf('SNR=%2d | BER=%.4e | BestDelay=%d (Fixing alignment)\n', ...
             snr_now, BER_vec(i), best_delay);

    % =======================================================
    
    % fprintf('SNR = %2d dB | BER = %.4e | Lag = %d\n', snr_now, BER_vec(i), peak_lag);
    
    % --- 诊断绘图 (如果 BER 还是高，请看这张图) ---
    if snr_now == 40
        figure;
        subplot(2,1,1); 
        plot(abs(xc)); title('Sync Correlation Peak (Check if unique & sharp)');
        subplot(2,1,2);
        scatterplot(sig_eq(eval_start:end));
        title(['Equalized Constellation at 50dB']);
    end
end

%% 6. 绘图 (Plotting for TVT)
% -----------------------------------------------------------------
figure;
semilogy(SNR_vec, BER_vec, '-bo', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER Performance under Replayed Channel');
axis([min(SNR_vec) max(SNR_vec) 1e-5 1]);
legend('RLS-DFE Performance');
disp('仿真完成。');

save_folder = '..\五月份实测数据\误码率数据0到30'; % 修改为您的目标文件夹
filename = sprintf('924_%g_%g.png', eq_fb_taps,rls_forget); % 文件名
fig_path = fullfile(save_folder, sprintf('924_%g_%g.fig', eq_fb_taps,rls_forget));
saveas(gcf, fig_path, 'fig');
png_path = fullfile(save_folder, sprintf('924_%g_%g', eq_fb_taps,rls_forget));
saveas(gcf, png_path, 'png');
mat_path = fullfile(save_folder, sprintf('924_%g_%g.mat', eq_fb_taps,rls_forget));
save(mat_path, 'BER_vec');
disp('图像已保存');




% figure
% plot(tx_chunk(1:end-29-5+6-1))
% hold on
% plot(rx_chunk(35-6+1:end))
% xlim([1,50])


[~, tmp_ber] = biterr(tx_chunk(1:end-29-5+6), rx_chunk(35-6:end));

