function stats = computeChannelStats(h_matrix, step_duration, fs)
%COMPUTECHANNELSTATS 从时变CIR矩阵计算信道统计特性
%  可复用的信道统计量计算函数，供 auto_parameter_selector.m、
%  TimeVeryCharater.m 等脚本调用。
%
%  Inputs:
%    h_matrix      - 时变CIR矩阵 [n_time × n_delay]
%    step_duration - 慢时间采样间隔 (秒)
%    fs            - 采样率 (Hz)
%
%  Outputs (struct):
%    .tau_rms      - RMS 时延扩展 (ms)
%    .tau_mean     - 平均时延 (ms)
%    .T_c          - 相干时间 (s), e^{-1} 准则
%    .B_c          - 相干带宽 (Hz), e^{-1} 准则
%    .B_c_tau      - 基于 tau_rms 估算的相干带宽 (Hz)
%    .f_d_mean     - 多普勒质心 (Hz)
%    .f_d_rms      - RMS 多普勒扩展 (Hz)
%    .rho_tau      - 时间自相关函数
%    .tau_axis_slow- 慢时间滞后轴 (s)
%    .rho_f        - 频率自相关函数
%    .freqLagAxis  - 频率滞后轴 (Hz)
%    .Doppler_axis - 多普勒频率轴 (Hz)
%    .S_fd         - 多普勒功率谱
%    .PDP_delay    - 功率时延谱
%    .delay_ms     - 时延轴 (ms)
%
%  See also: auto_parameter_selector, TimeVeryCharater

    stats = struct();

    %% 去除负时延
    seg_length = (size(h_matrix, 2) + 1) / 2;
    h_full = h_matrix(:, seg_length:end);
    Ndelay = size(h_full, 2);
    delay_ms = (0:Ndelay-1)' / fs * 1000;
    stats.delay_ms = delay_ms;

    %% 时间自相关 → 相干时间
    maxLag = min(50, size(h_full, 1) - 1);
    R_tau = zeros(maxLag + 1, Ndelay);
    for l = 1:Ndelay
        [c, ~] = xcorr(h_full(:, l), maxLag, 'biased');
        R_tau(:, l) = c(maxLag + 1:end);
    end
    R_tau_avg = mean(R_tau, 2);
    rho_tau = R_tau_avg / R_tau_avg(1);
    tau_axis_slow = (0:maxLag) * step_duration;

    target = exp(-1);
    idx_hi = find(rho_tau <= target, 1);
    if isempty(idx_hi) || idx_hi == 1
        T_c = NaN;
    else
        idx_lo = idx_hi - 1;
        T_c = interp1(rho_tau([idx_lo idx_hi]), tau_axis_slow([idx_lo idx_hi]), target, 'linear');
    end
    stats.rho_tau = rho_tau;
    stats.tau_axis_slow = tau_axis_slow;
    stats.T_c = T_c;

    %% 频率自相关 → 相干带宽
    Hfresp = fft(h_full, [], 2);
    maxFreqLag = min(250, Ndelay - 1);
    R_f = zeros(maxFreqLag + 1, size(h_full, 1));
    for t = 1:size(h_full, 1)
        [c, ~] = xcorr(Hfresp(t, :), maxFreqLag, 'biased');
        R_f(:, t) = c(maxFreqLag + 1:end);
    end
    df = fs / (2 * seg_length);
    freqLagAxis = (0:maxFreqLag) * df;
    R_f_avg = mean(R_f, 2);
    rho_f = real(R_f_avg / R_f_avg(1));

    idx_hif = find(rho_f <= target, 1);
    if isempty(idx_hif) || idx_hif == 1
        B_c = NaN;
    else
        idx_lof = idx_hif - 1;
        B_c = interp1(rho_f([idx_lof idx_hif]), freqLagAxis([idx_lof idx_hif]), target);
    end
    stats.rho_f = rho_f;
    stats.freqLagAxis = freqLagAxis;
    stats.B_c = B_c;

    %% 多普勒功率谱 → 多普勒扩展
    Nslow = size(h_full, 1);
    window = hamming(Nslow);
    Hf = fftshift(fft(h_full .* window, [], 1), 1);
    fs_slow = 1 / step_duration;
    Doppler_axis = (-Nslow / 2 : Nslow / 2 - 1) * (fs_slow / Nslow);
    S_fd = mean(abs(Hf).^2, 2);

    f_d_mean = sum(Doppler_axis(:) .* S_fd(:)) / sum(S_fd(:));
    f_d_rms = sqrt(sum(((Doppler_axis(:) - f_d_mean).^2) .* S_fd(:)) / sum(S_fd(:)));

    stats.Doppler_axis = Doppler_axis;
    stats.S_fd = S_fd;
    stats.f_d_mean = f_d_mean;
    stats.f_d_rms = f_d_rms;

    %% 功率时延谱 → 时延扩展
    PDP_delay = mean(abs(h_full).^2, 1);
    % 限制在 0-500ms 范围内
    delay_limit = min(floor(0.5 * fs) + 1, length(PDP_delay));
    PDP_trunc = PDP_delay(1:delay_limit);
    delay_trunc = delay_ms(1:delay_limit);

    tau_mean = sum(delay_trunc(:) .* PDP_trunc(:)) / sum(PDP_trunc(:));
    tau_rms = sqrt(sum(((delay_trunc(:) - tau_mean).^2) .* PDP_trunc(:)) / sum(PDP_trunc(:)));
    B_c_tau = sqrt(1 / (2 * pi^2 * (tau_rms / 1000)^2));

    stats.PDP_delay = PDP_delay;
    stats.tau_mean = tau_mean;
    stats.tau_rms = tau_rms;
    stats.B_c_tau = B_c_tau;

end
