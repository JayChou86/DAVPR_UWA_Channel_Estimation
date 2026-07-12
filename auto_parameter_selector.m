function params = auto_parameter_selector(h_matrix, fs, symbol_rate, step_duration, depth_label)
%AUTO_PARAMETER_SELECTOR  从信道统计量自动推导最优 DFE 均衡器参数
%
%  根据实测 CIR 矩阵计算信道统计量（时延扩展、相干时间、多普勒扩展等），
%  并基于物理推导规则自动给出 RLS-DFE 的最优参数建议。
%  替代之前的手工调参 + 深度查表方式。
%
%  Inputs:
%    h_matrix      - 时变 CIR 矩阵 [n_time × n_delay]
%    fs            - 采样率 (Hz)，默认 16000
%    symbol_rate   - 符号率 (Sps)，默认 4000
%    step_duration - CIR 时间分辨率 (s)，默认 0.1
%    depth_label   - 深度标签 (字符串)，用于输出显示，默认 ''
%
%  Output (struct):
%    .eq_ff_taps       - 推荐前馈抽头数
%    .eq_fb_taps       - 推荐反馈抽头数
%    .rls_forget       - 推荐 RLS 遗忘因子
%    .ref_tap          - 推荐参考抽头位置
%    .eq_ff_taps_lms   - LMS 等效前馈抽头数
%    .eq_fb_taps_lms   - LMS 等效反馈抽头数
%    .stats            - 信道统计量结构体 (来自 computeChannelStats)
%    .rls_flops_per_sym  - RLS-DFE 每符号乘法次数
%    .lms_flops_per_sym  - LMS-DFE 每符号乘法次数
%    .derivation_notes - 参数推导说明文本
%
%  See also: computeChannelStats, extract_multipath_by_threshold

    if nargin < 2 || isempty(fs),           fs = 16000;       end
    if nargin < 3 || isempty(symbol_rate),  symbol_rate = 4000; end
    if nargin < 4 || isempty(step_duration),step_duration = 0.1; end
    if nargin < 5 || isempty(depth_label),  depth_label = ''; end

    params = struct();
    T_sym = 1 / symbol_rate;

    %% 1. 计算信道统计量
    stats = computeChannelStats(h_matrix, step_duration, fs);
    params.stats = stats;

    %% 2. 推导 eq_fb_taps (反馈抽头数)
    % 反馈抽头需覆盖 RMS 时延扩展的主要多径能量
    % 公式: N_fb = ceil(tau_rms / T_sym) + margin
    % 其中 tau_rms 转换为秒: tau_rms 单位为 ms
    tau_rms_sec = stats.tau_rms / 1000;  % ms → s
    n_fb_cover = ceil(tau_rms_sec / T_sym);
    % 安全裕度: 至少 3 倍覆盖, 最少 30 抽头
    eq_fb_taps = max(30, n_fb_cover * 3);
    % 上限: 避免计算量过大
    eq_fb_taps = min(eq_fb_taps, 800);
    params.eq_fb_taps = eq_fb_taps;

    %% 3. 推导 rls_forget (RLS 遗忘因子)
    % λ 控制算法对信道变化的跟踪速度
    % 对于准静态信道 (T_c >> T_sym): λ → 1 (接近 1)
    % 对于快时变信道: λ 需降低以跟踪变化
    %
    % 经验公式: λ = exp(-T_sym / (alpha * T_c))
    %   其中 T_c 为相干时间 (s), alpha 为调节因子
    %   当 T_c 很大 (准静态), λ ≈ 1
    %   当 T_c 很小 (快变), λ 降低

    if isnan(stats.T_c) || isinf(stats.T_c)
        % 无法估计相干时间, 默认准静态
        rls_forget = 0.998;
    else
        alpha = 10;  % 调节因子: 越大则 λ 越接近 1 (更保守)
        rls_forget = exp(-T_sym / (alpha * stats.T_c));
        rls_forget = max(0.97, min(0.999, rls_forget));  % 约束在合理范围
    end
    params.rls_forget = rls_forget;

    %% 4. 推导 eq_ff_taps (前馈抽头数)
    % 前馈抽头主要覆盖: 主径到达之前的前驱能量 + 主径本身
    % 反馈抽头才负责消除后向多径 (已由 eq_fb_taps 覆盖)
    %
    % 策略: 找到 PDP 主峰, 向前搜索能量下降至-10dB的位置作为前驱范围

    pdp = stats.PDP_delay;
    delay_ms = stats.delay_ms;

    % 找到主峰位置
    [~, main_peak_idx] = max(pdp);
    main_peak_power = pdp(main_peak_idx);
    threshold_10dB = main_peak_power * 10^(-10/20);  % -10dB

    % 从主峰向左搜索, 找到能量降到-10dB以下的位置
    precursor_idx = main_peak_idx;
    for idx = main_peak_idx:-1:2
        if pdp(idx) < threshold_10dB
            precursor_idx = idx;
            break;
        end
    end

    % 前驱时延 (ms)
    precursor_delay_ms = delay_ms(main_peak_idx) - delay_ms(precursor_idx);

    % 前馈抽头数: 覆盖前驱 + 主峰宽度 (约2倍前驱作为安全裕度)
    n_ff_cover = ceil((precursor_delay_ms / 1000) / T_sym);
    eq_ff_taps = max(30, n_ff_cover * 3);          % 3倍裕度
    eq_ff_taps = min(eq_ff_taps, 120);              % 上限
    params.eq_ff_taps = eq_ff_taps;

    %% 5. 推导 ReferenceTap (参考抽头)
    % 参考抽头放在前馈滤波器中部偏后
    params.ref_tap = floor(eq_ff_taps / 2);

    %% 6. 复杂度估算
    N_rls = eq_ff_taps + eq_fb_taps;
    params.rls_flops_per_sym = 4 * N_rls^2 + 8 * N_rls;  % RLS ~O(N²)

    % LMS 需要更多抽头来达到相同性能
    eq_ff_taps_lms = eq_ff_taps * 2;
    eq_fb_taps_lms = eq_fb_taps * 2;
    N_lms = eq_ff_taps_lms + eq_fb_taps_lms;
    params.lms_flops_per_sym = 2 * N_lms + 2;            % LMS ~O(N)
    params.eq_ff_taps_lms = eq_ff_taps_lms;
    params.eq_fb_taps_lms = eq_fb_taps_lms;

    %% 7. 推导说明
    params.derivation_notes = sprintf([
        '[%s] 参数推导依据:\n', ...
        '  tau_rms = %.2f ms → eq_fb_taps = %d (覆盖 %.1f 倍 RMS 时延扩展)\n', ...
        '  T_c = %.1f ms → rls_forget = %.4f (准静态信道, lambda ≈ 1)\n', ...
        '  主峰前驱时延 = %.1f ms → eq_ff_taps = %d (覆盖前驱+主峰, -10dB准则)\n', ...
        '  RLS 复杂度: %.0f 乘法/符号 | LMS 复杂度: %.0f 乘法/符号\n', ...
        '  复杂度比 (LMS/RLS): %.2f%%\n'], ...
        depth_label, stats.tau_rms, eq_fb_taps, eq_fb_taps / max(1, n_fb_cover), ...
        stats.T_c * 1e3, rls_forget, ...
        precursor_delay_ms, eq_ff_taps, ...
        params.rls_flops_per_sym, params.lms_flops_per_sym, ...
        100 * params.lms_flops_per_sym / params.rls_flops_per_sym);

    %% 8. 输出显示
    fprintf('\n========== 自动参数选择结果 ==========\n');
    if ~isempty(depth_label)
        fprintf('深度: %s\n', depth_label);
    end
    fprintf('信道统计:\n');
    fprintf('  RMS 时延扩展: %.2f ms\n', stats.tau_rms);
    fprintf('  相干时间:     %.1f ms\n', stats.T_c * 1e3);
    fprintf('  相干带宽:     %.1f Hz\n', stats.B_c);
    fprintf('  RMS 多普勒:   %.3f Hz\n', stats.f_d_rms);
    fprintf('推荐参数:\n');
    fprintf('  前馈抽头:   %d\n', eq_ff_taps);
    fprintf('  反馈抽头:   %d\n', eq_fb_taps);
    fprintf('  参考抽头:   %d\n', params.ref_tap);
    fprintf('  遗忘因子:   %.4f\n', rls_forget);
    fprintf('计算复杂度:\n');
    fprintf('  RLS-DFE: %.0f 乘法/符号\n', params.rls_flops_per_sym);
    fprintf('  LMS-DFE: %.0f 乘法/符号 (%.1f%% of RLS)\n', ...
        params.lms_flops_per_sym, 100 * params.lms_flops_per_sym / params.rls_flops_per_sym);
    fprintf('========================================\n\n');

end
