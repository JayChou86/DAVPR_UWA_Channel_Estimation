file_name_ori='D:\A_Matlab_program\A_Matlab_program\五月份实验数据处理\signal48k2.wav';   %发射原始信号
file_name_rev='D:\A_Matlab_program\A_Matlab_program\时变信道估计代码\接收单道片段-0528\34-475.8m\475.8_0034_0368_442_32K_240528_000015_0s_01.wav';    %接收信号
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

plot(sig_ori)
[x, y] = ginput(2);

valid_peaks = abs(sig_ori(x(1):x(2))); %(15:56,351)
valid_peaks = abssig_ori(x(1):x(2)); %(15:56,351)

plotrician(valid_peaks, 10, 'Rician');  % exp Rician gamma weibull beta normal
