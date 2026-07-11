cc
% 代码汇总 624.8m 同步有问题 接近
%% 同步
max_window_duration = 3; % 窗长度, 单位: 秒
step_duration = 0.1;   % 滑动步长 (每隔多少秒滑动一次), 单位: 秒
file_name_ori='D:\A_Matlab_program\A_Matlab_program\五月份实验数据处理\signal48k2.wav';   %发射原始信号
file_name_rev='D:\A_Matlab_program\A_Matlab_program\时变信道估计代码\接收单道片段-0528\11-624.8m\624.8_0011_0366_367_32K_240528_000010_0a_01.wav';    %接收信号
file_name_bd='D:\A_Matlab_program\A_Matlab_program\五月份实验数据处理\0007_0000_0_32K_240528_004922_71a_01.wav';   % 标定信号
[sig_ori,fs1]=audioread(file_name_ori);
[sig_rev,fs2]=audioread(file_name_rev);
[sig_bd,fs3]=audioread(file_name_bd);
gcd_ab = gcd(fs1, fs2);     % 最大公约数
gcd_abc = gcd(gcd_ab, fs3);
fs = gcd_abc;
dt = 1/fs;

% 下采样
sig_ori = downsample(sig_ori, fs1 / fs);
sig_rev = downsample(sig_rev, fs2 / fs);
sig_bd = downsample(sig_bd, fs3 / fs);


fl = 600;
fh = 1000;
Wn = [fl fh] / (fs/2);
[b, a] = butter(4, Wn, 'bandpass');
sig_rev = filtfilt(b, a, sig_rev);
sig_bd = filtfilt(b, a, sig_bd);

% 截取信号长度，确保二者长度一致
min_length = min(length(sig_rev), length(sig_bd));
sig_rev = sig_rev(end-min_length+1:end);
sig_bd = sig_bd(end-min_length+1:end);

% 滤波
fl = 600;
fh = 1000;
Wn = [fl fh] / (fs/2);
[b, a] = butter(4, Wn, 'bandpass');
sig_rev = filtfilt(b, a, sig_rev);
sig_bd = filtfilt(b, a, sig_bd);



% 计算互相关
[rxy, lags] = xcorr(sig_rev, sig_bd);
[~, max_idx1] = max(abs(rxy));
time_shift_samples_rev_bd = lags(max_idx1); % 接收信号相对于发射信号的时间延迟

% 根据时间延迟调整接收信号和标定信号
if time_shift_samples_rev_bd > 0
    sig_rev_sync = [sig_rev(time_shift_samples_rev_bd + 1:end); zeros(time_shift_samples_rev_bd, 1)];
elseif time_shift_samples_rev_bd < 0
    sig_rev_sync = [zeros(-time_shift_samples_rev_bd, 1); sig_rev(1:end + time_shift_samples_rev_bd)];
else
    sig_rev_sync = sig_rev;
end
sig_bd_sync = sig_bd;

% plotfft(sig_rev, 1)
% plotfft(sig_bd, 1)
% plotfft(sig_rev_sync, 1)
% plotfft(sig_bd_sync, 1)



% 截出m序列
sig_bd_m = sig_bd_sync(1037410:1202230);     
sig_rev_m = sig_rev_sync(1037410:1202230);     

% cor(sig_rev_m, sig_bd_m, 1)

%% 时变信道估计

% plotfft(sig_rev,fs)

for window_duration = 1:max_window_duration
    window_length = round(window_duration * fs); % 窗长度对应的采样点数
    step_length = round(step_duration * fs);     % 每次滑动的采样点数 (也叫 hop size)

    % 为保证两个信号分段后片段数量一致，可以取两者中较短的长度
    trc_len = min(length(sig_bd_m), length(sig_rev_m));
    sig_bd_re = sig_bd_m(1:trc_len);
    sig_rev_re = sig_rev_m(1:trc_len);

    % 对 sig_bd_re 信号进行滑动窗分段，计算所有窗口的起始点
    start_indices = 1:step_length:(length(sig_bd_re) - window_length + 1);
    seg_num = length(start_indices); % 计算总共能产生多少个片段

    % 初始化结果矩阵 (预分配内存以提高效率)，每列是一个片段
    sig_bd_resp = zeros(window_length, seg_num);

    % 循环提取每个片段
    for i = 1:seg_num
        start_idx = start_indices(i);
        end_idx = start_idx + window_length - 1;
        sig_bd_resp(:, i) = sig_bd_re(start_idx:end_idx);
    end

    % 由于长度已经对齐，可以直接使用相同的参数和循环
    sig_rev_resp = zeros(window_length, seg_num);
    for i = 1:seg_num
        start_idx = start_indices(i);
        end_idx = start_idx + window_length - 1;
        sig_rev_resp(:, i) = sig_rev_re(start_idx:end_idx);
    end


    % 使用基于FFT的方法高效计算所有片段的互相关
    N_fft = 2^nextpow2(2 * window_length - 1);

    % 对所有片段一次性进行FFT (沿列操作)
    F_rev = fft(sig_rev_resp, N_fft, 1);
    F_bd = fft(sig_bd_resp, N_fft, 1);

    % 计算互相关的频域表示 (点乘共轭)
    Rxy_freq = F_rev .* conj(F_bd);

    % 将结果转换回时域，得到互相关矩阵
    rxy_matrix_fft = ifft(Rxy_freq, [], 1);

    % xcorr的结果需要fftshift来将零延迟点移动到中心，并截取有效部分
    rxy_matrix = fftshift(rxy_matrix_fft, 1);
    rxy_matrix = rxy_matrix( (N_fft/2 - (window_length-1)) : (N_fft/2 + (window_length-1)) , : );
    rxy_matrix = rxy_matrix / window_length;        % 幅度

    % 注意：MATLAB的xcorr输出是行向量，而我们这里是列向量，需要转置
    h_matrix = rxy_matrix';

    % 修正时间轴: 每个窗口的中心时间点
    num_segments = size(sig_bd_resp, 2);
    time_sec = (0:num_segments-1) * step_duration; % 使用滑动步长时间

    % 时延轴 (单位: ms)
    delay_lags = -(window_length-1):(window_length-1);
    delay_ms_all = delay_lags' / fs * 1000;

    % 归一化处理 (选择一种即可) 将整个矩阵的最大绝对值归一化为1，便于观察整体响应
    h_matrix = h_matrix / max(abs(h_matrix(:)));

    % 绘制热力图
    figure;
    imagesc(delay_ms_all, time_sec, abs(h_matrix));
    ylabel('时间 (s)');
    xlabel('时延 (ms)');
    title(sprintf('滑动步长为%g,窗长度为%g的时变冲激响应', step_duration,window_duration));colorbar;
    colorbar;
    axis xy;
    colormap('jet');
    xlim([-50,500])

    % save_folder = '..\五月份实测数据\图像\624'; % 修改为您的目标文件夹
    % filename = sprintf('滑动步长为%g_窗长度为%g.png', step_duration,window_duration); % 文件名
    % fig_path = fullfile(save_folder, sprintf('滑动步长为%g_窗长度为%g.fig', step_duration,window_duration));
    % saveas(gcf, fig_path, 'fig');
    % png_path = fullfile(save_folder, sprintf('滑动步长为%g_窗长度为%g.png', step_duration,window_duration));
    % saveas(gcf, png_path, 'png');
    % disp('图像已保存');
end



% %% 时间自相关函数
% maxLag = 20;                                   % 时间维度最大滞后点
% Ndelay = size(h_matrix,2);
% R_tau = zeros(maxLag+1, Ndelay);
% R_tau_all = zeros(2*maxLag+1, Ndelay);
% for l = 1:Ndelay
%     [c,~] = xcorr(h_matrix(:,l), maxLag, 'biased');
%     R_tau(:,l) = c(maxLag+1:end);              % 非负滞后
%     R_tau_all(:,l) = c; 
% end
% R_tau_avg = mean(R_tau, 2);
% rho_tau = R_tau_avg / R_tau_avg(1);
% dt_slow = seg_length/fs;                       % 慢时间采样间隔
% tau_axis_slow = (0:maxLag) * dt_slow;          % 慢时间滞后轴（秒）
% 
% % 频域互相关
% Hfresp = fft(h_matrix, [], 2);                % H(f,t)
% maxFreqLag = 250;
% R_f = zeros(maxFreqLag+1, size(h_matrix,1));
% for t = 1:size(h_matrix,1)
%     [c,~] = xcorr(Hfresp(t,:), maxFreqLag, 'biased');
%     R_f(:,t) = c(maxFreqLag+1:end);
% end
% df = fs / (2*seg_length);                     % 频率间隔（Hz）
% freqLagAxis = (0:maxFreqLag) * df;
% R_f_avg = mean(R_f,2);
% rho_f = real(R_f_avg / R_f_avg(1));
% 
% % 时间相关函数热力图
% tau0 = 0;
% t_axis = (-maxLag:maxLag) * dt_slow;
% tau_axis_full = tau0 + (0:Ndelay-1) * dt;
% figure;
% imagesc(tau_axis_full * 1000, t_axis, abs(R_tau_all));
% set(gca, 'YDir', 'normal');
% xlabel('时延 (ms)');
% ylabel('时间 (s)');
% % title('R_\tau heatmap');
% colorbar;
% 
% % 相干时间
% target = 0.5;
% idx_hi = find(rho_tau <= target, 1);
% if isempty(idx_hi) || idx_hi == 1
%     T_c = NaN;
% else
%     idx_lo = idx_hi - 1;
%     x = tau_axis_slow([idx_lo idx_hi]);
%     y = rho_tau([idx_lo idx_hi]);
%     T_c = interp1(y, x, target, 'linear');
% end
% fprintf('相干时间 = %.2f ms\n', T_c*1e3);
% figure;
% plot(tau_axis_slow*1e3, rho_tau, 'LineWidth', 1.5); grid on;
% xlabel('时间 (ms)');
% ylabel('归一化幅度');
% % title('Channel time autocorrelation');
% yline(exp(-1), '--k', 'e^{-1}');
% yline(0.5, '--k', '0.5');
% yline(0.9, '--k', '0.9');
% save('624time_cor.mat', 'tau_axis_slow','rho_tau');
% 
% 
% % 相干带宽
% targetf = exp(-1);
% idx_hif = find(rho_f <= targetf, 1);
% if isempty(idx_hif) || idx_hif == 1
%     B_c = NaN;
% else
%     idx_lof = idx_hif - 1;
%     B_c = interp1(rho_f([idx_lof idx_hif]), freqLagAxis([idx_lof idx_hif]), targetf);
% end
% fprintf('相干带宽 = %.2f Hz\n', B_c);
% figure;
% plot(freqLagAxis, rho_f, 'LineWidth', 1.5);
% hold on;
% yline(exp(-1), '--k', 'e^{-1}');
% yline(0.5, '--k', '0.5');
% yline(0.9, '--k', '0.9');
% xlabel('频率 (Hz)');
% ylabel('归一化幅度');
% % title('平均频率自相关函数');
% grid on;
% save('624band_cor.mat', 'freqLagAxis','rho_f');
% 
% 
% %% === 多普勒功率谱与扩展 ===
% Nslow = size(h_matrix,1);
% window = hamming(Nslow);                      % 减小频谱泄漏
% Hf = fftshift(fft(h_matrix .* window, [], 1), 1);  % 先按慢时间方向FFT
% Doppler_axis = (-Nslow/2:Nslow/2-1) * (fs/seg_length/Nslow);  % 多普勒频率轴 (Hz)
% S_fd = mean(abs(Hf).^2, 2);                   % 各时延抽头平均后的多普勒功率谱
% 
% figure;
% plot(Doppler_axis, 10*log10(S_fd/max(S_fd)), 'LineWidth', 1.5);
% xlabel('多普勒频率 (Hz)');
% ylabel('归一化功率 (dB)');
% % title('多普勒功率谱');
% grid on;
% 
% % 多普勒扩展统计量
% fd_mean = sum(Doppler_axis(:) .* S_fd(:)) / sum(S_fd);                % 多普勒质心
% fd_rms = sqrt(sum(((Doppler_axis(:) - fd_mean).^2).*S_fd(:)) / sum(S_fd)); % RMS 多普勒宽度
% fd_rms_rad = 2*pi*fd_rms;                                            % 转换为 rad/s
% fprintf('多普勒质心 = %.3f Hz, RMS 多普勒扩展 = %.3f Hz (%.3f rad/s)\n', ...
%         fd_mean, fd_rms, fd_rms_rad);
% 
% %% 多普勒时延谱
% % === 多普勒-时延功率谱 (Delay–Doppler Map) ===
% % Hf 已在慢时间方向 FFT 过：Hf ≈ H(f_d, τ)
% S_tau_fd = abs(Hf).^2;                        % 功率
% S_tau_fd = S_tau_fd ./ max(S_tau_fd(:));      % 归一化
% 
% figure;
% imagesc(delay_axis*1e3, Doppler_axis, (S_tau_fd));
% axis xy;
% xlabel('时延 (ms)');
% ylabel('多普勒频率 (Hz)');
% % title('多普勒-时延功率谱');
% c = colorbar;
% c.Label.String = '功率 (dB)';
% % colormap('jet');
% 
% %% === 功率时延谱 (PDP) ===
% PDP_delay = mean(abs(h_matrix).^2, 1);        % 对慢时间取平均
% figure;
% plot(delay_axis(1:numel(PDP_delay))*1e3, PDP_delay / max(PDP_delay), 'LineWidth', 1.5);
% xlabel('时延 (ms)');
% ylabel('功率');
% % title('功率时延谱 (线性)');
% grid on;
% 
% % figure;
% % plot(delay_axis(1:numel(PDP_delay))*1e3, 10*log10(PDP_delay/max(PDP_delay)), 'LineWidth', 1.5);
% % xlabel('时延 (ms)');
% % ylabel('功率 (dB)');
% % title('功率时延谱 (归一化 dB)');
% % grid on;
% 
% % 平均时延与RMS时延扩展
% tau_mean = sum(delay_axis(1:numel(PDP_delay')) .* PDP_delay') / sum(PDP_delay');
% tau_rms = sqrt(sum(((delay_axis(1:numel(PDP_delay')) - tau_mean).^2) .* PDP_delay') / sum(PDP_delay'));
% fprintf('平均超前时延 = %.3f ms, RMS 时延扩展 = %.3f ms\n', tau_mean*1e3, tau_rms*1e3);
% 
% % 基于 RMS 时延估算的相干带宽
% Bc_tau = sqrt(1/(2*(pi^2)*(tau_rms^2)));
% fprintf('基于RMS时延扩展估算的相干带宽 Bc_tau = %.2f Hz\n', Bc_tau);