cc
%% 475、624、824、924
load('时变冲激响应924_0.1_3.mat');
save_folder = '..\五月份实测数据\图像\时变函数图像\924'; % 修改目标文件夹
mat_name = 'TimeVery_924.mat';

num_segments = size(h_matrix,1);
step_duration = 0.1;
fs = 16000;
step_length = step_duration * fs;
seg_length = (size(h_matrix,2) + 1)/2;
h = h_matrix(:, seg_length:end);     % 去除负时延
delay_ms = (0:seg_length-1)' / fs * 1000;

figure;
imagesc(delay_ms, time_sec, abs(h));
ylabel('时间 (s)');
xlabel('时延 (ms)');
title('时变冲激响应');
colorbar;
axis xy;
colormap('jet');
xlim([0,500])

maxLag = 50;                                   % 时间维度最大滞后点
Ndelay = size(h_matrix,2);
R_tau = zeros(maxLag+1, Ndelay);
R_tau_all = zeros(2*maxLag+1, Ndelay);
for l = 1:Ndelay
    [c,~] = xcorr(h_matrix(:,l), maxLag, 'biased');
    R_tau(:,l) = c(maxLag+1:end);              % 非负滞后
    R_tau_all(:,l) = c;
end
R_tau_avg = mean(R_tau, 2);
rho_tau = R_tau_avg / R_tau_avg(1);
dt_slow = step_duration;                       % 慢时间采样间隔
tau_axis_slow = (0:maxLag) * dt_slow;          % 慢时间滞后轴（秒）

% 频域互相关
Hfresp = fft(h, [], 2);                % H(f,t)
maxFreqLag = 250;
R_f = zeros(maxFreqLag+1, size(h,1));
for t = 1:size(h,1)
    [c,~] = xcorr(Hfresp(t,:), maxFreqLag, 'biased');
    R_f(:,t) = c(maxFreqLag+1:end);
end
df = fs / (2*seg_length);                     % 频率间隔（Hz）
freqLagAxis = (0:maxFreqLag) * df;
R_f_avg = mean(R_f,2);
rho_f = real(R_f_avg / R_f_avg(1));

% 时间相关函数热力图
dt = 1/fs;
tau0 = 0;
t_axis = (-maxLag:maxLag) * dt_slow;
tau_axis_full = tau0 + (-Ndelay:Ndelay-1) * dt;
figure;
imagesc(tau_axis_full*1000, t_axis, abs(R_tau_all));
set(gca, 'YDir', 'normal');
xlabel('时延 (ms)');
ylabel('时间(s)');
title('时间相关函数热力图');
colorbar;
colormap('jet');
xlim([-20,500]);

% filename = sprintf('时间相关函数热力图.png'); % 文件名
% fig_path = fullfile(save_folder, sprintf('时间相关函数热力图.fig'));
% saveas(gcf, fig_path, 'fig');
% png_path = fullfile(save_folder, sprintf('时间相关函数热力图.fig.png'));
% saveas(gcf, png_path, 'png');
% disp('图像已保存');

% 相干时间
target = exp(-1);
idx_hi = find(rho_tau <= target, 1);
if isempty(idx_hi) || idx_hi == 1
    T_c = NaN;
else
    idx_lo = idx_hi - 1;
    x = tau_axis_slow([idx_lo idx_hi]);
    y = rho_tau([idx_lo idx_hi]);
    T_c = interp1(y, x, target, 'linear');
end
fprintf('相干时间 = %.2f ms\n', T_c*1e3);
figure;
plot(tau_axis_slow*1e3, rho_tau, 'LineWidth', 1.5); grid on;
xlabel('时间(ms)');
ylabel('归一化幅度');
title('平均时间自相关函数');
yline(exp(-1), '--k', 'e^{-1}');
yline(0.5, '--k', '0.5');
yline(0.9, '--k', '0.9');

% filename = sprintf('平均时间自相关函数.png'); % 文件名
% fig_path = fullfile(save_folder, sprintf('平均时间自相关函数.fig'));
% saveas(gcf, fig_path, 'fig');
% png_path = fullfile(save_folder, sprintf('平均时间自相关函数.png'));
% saveas(gcf, png_path, 'png');
% save('475time_cor.mat', 'tau_axis_slow','rho_tau');

% 相干带宽
targetf = exp(-1);
idx_hif = find(rho_f <= targetf, 1);
if isempty(idx_hif) || idx_hif == 1
    B_c = NaN;
else
    idx_lof = idx_hif - 1;
    B_c = interp1(rho_f([idx_lof idx_hif]), freqLagAxis([idx_lof idx_hif]), targetf);
end
fprintf('相干带宽 = %.2f Hz\n', B_c);
figure;
plot(freqLagAxis, rho_f, 'LineWidth', 1.5);
hold on;
yline(exp(-1), '--k', 'e^{-1}');
yline(0.5, '--k', '0.5');
yline(0.9, '--k', '0.9');
xlabel('频率 (Hz)');
ylabel('归一化幅度');
title('平均频率自相关函数');
grid on;

% filename = sprintf('平均频率自相关函数.png'); % 文件名
% fig_path = fullfile(save_folder, sprintf('平均频率自相关函数.fig'));
% saveas(gcf, fig_path, 'fig');
% png_path = fullfile(save_folder, sprintf('平均频率自相关函数.png'));
% saveas(gcf, png_path, 'png');
% save('475band_cor.mat', 'freqLagAxis','rho_f');


%% === 多普勒功率谱与扩展 ===
Nslow = size(h,1);
window = hamming(Nslow);                      % 减小频谱泄漏
Hf = fftshift(fft(h .* window, [], 1), 1);  % 先按慢时间方向FFT
% Doppler_axis = (-Nslow/2:Nslow/2-1) * (fs/seg_length/Nslow);  % 多普勒频率轴 (Hz)
fs_slow = 1 / step_duration; % 慢时间采样率 (Hz)
Doppler_axis = (-Nslow/2:Nslow/2-1) * (fs_slow / (step_duration * fs)); % 多普勒频率轴
S_fd = mean(abs(Hf).^2, 2);       % 各时延抽头平均后的多普勒功率谱

figure;
plot(Doppler_axis, 10*log10(S_fd/max(S_fd)), 'LineWidth', 1.5);
xlabel('多普勒频率 (Hz)');
ylabel('归一化功率 (dB)');
title('多普勒功率谱');
grid on;

% filename = sprintf('多普勒功率谱.png'); % 文件名
% fig_path = fullfile(save_folder, sprintf('多普勒功率谱.fig'));
% saveas(gcf, fig_path, 'fig');
% png_path = fullfile(save_folder, sprintf('多普勒功率谱.png'));
% saveas(gcf, png_path, 'png');

% 多普勒扩展统计量
fd_mean = sum(Doppler_axis(:) .* S_fd(:)) / sum(S_fd);                % 多普勒质心
fd_rms = sqrt(sum(((Doppler_axis(:) - fd_mean).^2).*S_fd(:)) / sum(S_fd)); % RMS 多普勒宽度
fd_rms_rad = 2*pi*fd_rms;                                            % 转换为 rad/s
fprintf('多普勒质心 = %.3f Hz, RMS 多普勒扩展 = %.3f Hz (%.3f rad/s)\n', ...
        fd_mean, fd_rms, fd_rms_rad);

%% 多普勒时延谱
% === 多普勒-时延功率谱 (Delay–Doppler Map) ===
% Hf 已在慢时间方向 FFT 过：Hf ≈ H(f_d, τ)
delay_axis = (0:seg_length)' * dt;                    % 秒
S_tau_fd = abs(Hf).^2;                        % 功率
S_tau_fd = S_tau_fd ./ max(S_tau_fd(:));      % 归一化

figure;
imagesc(delay_axis*1e3, Doppler_axis, (S_tau_fd));
axis xy;
xlabel('时延 (ms)');
ylabel('多普勒频率 (Hz)');
title('多普勒-时延功率谱');
c = colorbar;
c.Label.String = '功率 (dB)';
xlim([0,500])
colormap('jet');

% filename = sprintf('多普勒-时延功率谱.png'); % 文件名
% fig_path = fullfile(save_folder, sprintf('多普勒-时延功率谱.fig'));
% saveas(gcf, fig_path, 'fig');
% png_path = fullfile(save_folder, sprintf('多普勒-时延功率谱.png'));
% saveas(gcf, png_path, 'png');

%% === 功率时延谱 (PDP) ===
PDP_delay = mean(abs(h).^2, 1);        % 对慢时间取平均
figure;
plot(delay_axis(1:numel(PDP_delay))*1e3, PDP_delay/max(PDP_delay), 'LineWidth', 1.5);
xlabel('时延 (ms)');
ylabel('归一化功率');
title('功率时延谱');
xlim([0,500])
grid on;

% filename = sprintf('功率时延谱.png'); % 文件名
% fig_path = fullfile(save_folder, sprintf('功率时延谱.fig'));
% saveas(gcf, fig_path, 'fig');
% png_path = fullfile(save_folder, sprintf('功率时延谱.png'));
% saveas(gcf, png_path, 'png');
% figure;
% plot(delay_axis(1:numel(PDP_delay))*1e3, 10*log10(PDP_delay/max(PDP_delay)), 'LineWidth', 1.5);
% xlabel('时延 (ms)');
% ylabel('功率 (dB)');
% title('功率时延谱 (归一化 dB)');
% grid on;

% 平均时延与RMS时延扩展
tau_mean = sum(delay_ms(1:500*fs/1000 + 1) .* PDP_delay(1:500*fs/1000 + 1)') / sum(PDP_delay(1:500*fs/1000 + 1));
tau_rms = sqrt(sum(((delay_ms(1:500*fs/1000 + 1) - tau_mean).^2) .* PDP_delay(1:500*fs/1000 + 1)') / sum(PDP_delay(1:500*fs/1000 + 1)'));
 
fprintf('平均超前时延 = %.3f ms, RMS 时延扩展 = %.3f ms\n', tau_mean, tau_rms);

% 基于 RMS 时延估算的相干带宽
Bc_tau = sqrt(1/(2*(pi^2)*((tau_rms/1000)^2)));
fprintf('基于RMS时延扩展估算的相干带宽 Bc_tau = %.2f Hz\n', Bc_tau);


% save(mat_name);














