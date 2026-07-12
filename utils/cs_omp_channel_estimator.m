function [h_omp, support_info] = cs_omp_channel_estimator(sig_bd_resp, sig_rev_resp, varargin)
%CS_OMP_CHANNEL_ESTIMATOR  基于 OMP 压缩感知的稀疏信道估计
%
%  利用 m序列 的优良自相关特性作为测量矩阵，通过 OMP 算法逐个窗口
%  估计稀疏水声信道冲激响应。相比传统的匹配滤波（FFT互相关），OMP能：
%    - 提供超分辨率（分辨间距小于 1/BW 的多径）
%    - 天然降噪（仅保留 K 个最强原子）
%    - 给出更清晰的 PDP（零旁瓣）
%
%  Inputs:
%    sig_bd_resp  - 发射信号窗口矩阵 [L_window × N_segments]
%    sig_rev_resp - 接收信号窗口矩阵 [L_window × N_segments]
%
%  Options (Name-Value pairs):
%    'K_max'      - 最大稀疏度 (每窗口最大径数), 默认 80
%    'stop_residual' - 残差能量比停止阈值, 默认 0.05 (5%)
%    'verbose'    - 是否显示进度, 默认 true
%
%  Outputs:
%    h_omp        - OMP 估计的 CIR 矩阵 [N_segments × L_window]
%    support_info - 结构体: .K_per_window, .residual_ratio, .runtime
%
%  Algorithm:
%    接收信号模型: y = S * h + n
%      其中 S 是发射信号的 Toeplitz 卷积矩阵
%    OMP 迭代选择与残差最相关的 S 的列（即某时延处的 m 序列副本）
%    每次迭代通过最小二乘更新已选原子对应的系数
%
%  Reference:
%    Berger, Zhou, et al. "Sparse channel estimation for multicarrier
%    underwater acoustic communication." IEEE JOE, 2010.
%
%  See also: window_475fenduan, compare_mf_vs_omp

%% Parse inputs
p = inputParser;
p.addParameter('K_max', 80, @(x) isscalar(x) && x > 0);
p.addParameter('stop_residual', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
p.addParameter('verbose', true, @islogical);
p.parse(varargin{:});
opts = p.Results;

K_max = opts.K_max;
stop_ratio = opts.stop_residual;
verbose = opts.verbose;

[L_window, N_seg] = size(sig_bd_resp);
h_omp = zeros(N_seg, L_window);
K_per_window = zeros(N_seg, 1);
residual_ratios = zeros(N_seg, 1);

if verbose
    fprintf('[OMP] 窗口长度=%d, 片段数=%d, K_max=%d\n', L_window, N_seg, K_max);
end

tic_total = tic;

%% 逐窗口 OMP 估计
for seg = 1:N_seg
    y = sig_rev_resp(:, seg);     % 接收信号 (观测)
    x = sig_bd_resp(:, seg);      % 发射信号 (m序列)

    y_norm = norm(y);
    if y_norm < 1e-12
        K_per_window(seg) = 0;
        continue;
    end

    % --- OMP 迭代 ---
    residual = y;
    support = [];                  % 已选时延索引
    h_sparse = zeros(L_window, 1); % 稀疏 CIR

    for k = 1:K_max
        % 步骤1: 计算所有原子与残差的相关性
        % 原子 i = x 循环移位 i (即 x shifted by i)
        % corr_i = <atom_i, residual>
        % 通过 FFT 高效计算:
        L_fft = 2^nextpow2(2 * L_window - 1);
        corr_full = ifft(fft(residual, L_fft) .* conj(fft(x, L_fft)));
        corr = corr_full(1:L_window);  % 正时延部分

        % 步骤2: 选择最相关的原子 (排除已选项)
        [~, corr_abs_sorted_idx] = sort(abs(corr), 'descend');
        new_idx = [];
        for c = 1:length(corr_abs_sorted_idx)
            if ~ismember(corr_abs_sorted_idx(c), support)
                new_idx = corr_abs_sorted_idx(c);
                break;
            end
        end

        if isempty(new_idx)
            break;  % 所有可能时延都已选
        end

        support = [support; new_idx];
        support = sort(support);
        K = length(support);

        % 步骤3: 构建测量子矩阵并最小二乘求解
        % X_support 的每一列是 x 移位相应时延后的向量
        X_sub = zeros(L_window, K);
        for j = 1:K
            tau = support(j);  % 时延 (MATLAB 1-indexed)
            if tau == 1
                X_sub(:, j) = x;
            elseif tau <= L_window
                X_sub(:, j) = [x(tau:end); zeros(tau-1, 1)];
            else
                X_sub(:, j) = zeros(L_window, 1);
            end
        end

        % 最小二乘: h_s = (X'X)^{-1} X' y
        XtX = X_sub' * X_sub;
        Xty = X_sub' * y;

        % 正则化: 当 X'X 接近奇异时加小对角线扰动
        if rcond(XtX) < 1e-10
            h_support = (XtX + 1e-6 * eye(K)) \ Xty;
        else
            h_support = XtX \ Xty;
        end

        % 步骤4: 更新残差
        residual = y - X_sub * h_support;

        % 步骤5: 收敛检查
        residual_ratio = norm(residual) / y_norm;
        if residual_ratio < stop_ratio
            break;
        end
    end

    % 存储结果
    for j = 1:length(support)
        h_sparse(support(j)) = h_support(j);
    end

    h_omp(seg, :) = h_sparse;
    K_per_window(seg) = length(support);
    residual_ratios(seg) = norm(residual) / y_norm;

    % 进度
    if verbose && (mod(seg, ceil(N_seg/10)) == 0 || seg == N_seg)
        fprintf('  [OMP] 进度: %d/%d (%.0f%%), 平均径数=%.1f\n', ...
            seg, N_seg, seg/N_seg*100, mean(K_per_window(1:seg)));
    end
end

runtime = toc(tic_total);

%% 输出
support_info = struct();
support_info.K_per_window = K_per_window;
support_info.K_mean = mean(K_per_window);
support_info.K_std = std(K_per_window);
support_info.residual_ratio_mean = mean(residual_ratios);
support_info.residual_ratios = residual_ratios;
support_info.runtime = runtime;

if verbose
    fprintf('[OMP] 完成! 总耗时=%.2fs, 平均径数=%.1f±%.1f, 残差比=%.3f\n', ...
        runtime, support_info.K_mean, support_info.K_std, ...
        support_info.residual_ratio_mean);
    fprintf('[OMP] 压缩比: %.1f%% (平均保留 %.0f/%d 非零系数)\n', ...
        support_info.K_mean / L_window * 100, support_info.K_mean, L_window);
end

end
