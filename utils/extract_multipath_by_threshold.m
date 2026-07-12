function [h_extracted, keep_indices, pdp_norm] = extract_multipath_by_threshold(h_matrix, threshold_dB, min_distance)
%EXTRACT_MULTIPATH_BY_THRESHOLD 基于能量阈值自动提取有效多径抽头
%   替代手工挑选特定列，基于平均功率时延谱保留-20dB以上的所有有效多径
%
%   Inputs:
%     h_matrix      - 时变CIR矩阵 [n_time × n_delay]，已经过精确同步对齐
%     threshold_dB  - 能量阈值 (dB)，默认 -20 (即保留峰值-20dB以上的抽头)
%     min_distance   - 最小抽头间距 (样本点)，默认 5 (合并相邻径)
%
%   Outputs:
%     h_extracted   - 截取后的CIR矩阵 [n_time × n_kept]，仅含有效多径
%     keep_indices  - 保留的原始列索引
%     pdp_norm       - 归一化平均功率时延谱 [1 × n_delay]
%
%   Usage:
%     [h_cut, idx] = extract_multipath_by_threshold(result_matrix, -20);
%
%   See also: auto_parameter_selector, Channel_review_475

if nargin < 2 || isempty(threshold_dB)
    threshold_dB = -20;  % 默认 -20 dB 阈值
end
if nargin < 3 || isempty(min_distance)
    min_distance = 5;    % 默认最小间距 5 个采样点
end

[n_time, n_delay] = size(h_matrix);

%% 1. 计算归一化平均功率时延谱
avg_pdp = mean(abs(h_matrix).^2, 1);
pdp_norm = avg_pdp / max(avg_pdp);  % 归一化到 [0,1]
pdp_dB = 10 * log10(pdp_norm + eps);  % 转 dB

%% 2. 找出高于阈值的区域
above_threshold = pdp_dB >= threshold_dB;

% 使用形态学闭运算合并间距过小的区域
if min_distance > 1
    se = ones(1, min_distance);
    above_threshold = imclose(double(above_threshold), se) > 0;
end

%% 3. 提取连续区域并找到每个区域内的峰值
% 找到阈值以上的所有连续段
d = diff([0, above_threshold, 0]);
start_segments = find(d == 1);
end_segments = find(d == -1) - 1;

keep_indices = [];
for s = 1:length(start_segments)
    seg_start = start_segments(s);
    seg_end = end_segments(s);
    segment_pdp = pdp_norm(seg_start:seg_end);
    % 在该段内找到峰值位置
    [~, local_peak] = max(segment_pdp);
    peak_idx = seg_start + local_peak - 1;
    % 保留该峰值
    if ~ismember(peak_idx, keep_indices)
        keep_indices(end+1) = peak_idx;
    end
end

%% 4. 如果未找到任何有效径，回退到最强的3个
if isempty(keep_indices)
    warning('extract_multipath:noMultipath', ...
        '阈值 %.1f dB 下未找到有效多径，使用前3强峰值。', threshold_dB);
    [~, sorted_idx] = sort(pdp_norm, 'descend');
    keep_indices = sorted_idx(1:min(3, n_delay));
end

%% 5. 排序并保留
keep_indices = sort(keep_indices);

%% 6. 提取子矩阵
% 为了保持时延连续性，也保留峰值之间的区域
% 这里保留从第一个峰值到最后一个峰值的全部区域
first_keep = keep_indices(1);
last_keep = keep_indices(end);

% 保留区间并向外扩展一点边界
margin = 10;
keep_start = max(1, first_keep - margin);
keep_end = min(n_delay, last_keep + margin);
keep_range = keep_start:keep_end;

h_extracted = h_matrix(:, keep_range);

fprintf('[extract_multipath] 阈值 %.1f dB: %d 个有效径 | 从 %d 列截取到 %d 列 (保留 %d 个多径峰值)\n', ...
    threshold_dB, length(keep_indices), n_delay, length(keep_range), length(keep_indices));

end
