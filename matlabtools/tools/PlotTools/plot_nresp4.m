function plot_nresp4(data, zz)

nresp = 4

figure('Position', [100 300 1200 300]);

[ha, pos] = tight_subplot(1, 4, 0.01, 0.1,  0.05);
%set(ha(2:4), 'ycolor', 'none');

for i = 1:nresp
    axes(ha(i));
    imagesc(sqrt(abs(data(:,:,zz,i))));
    hold on
    colormap('gray');
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
end

sgtitle('Esophagus in respiratory phases (gridding solution, 256x256)');
