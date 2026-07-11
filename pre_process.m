cc
load('时变冲激响应475_0.1_3.mat')

seg_length = (size(h_matrix,2) + 1)/2;

h = h_matrix(:, seg_length:end);     % 去除负时延

imagesc(h);xlim([0,500])

plot(h(1,1:1000));
