function plot_nresp4_coronal(data)%, zz)

nresp = 4

figure('Position', [100 300 1200 300]);

[ha, pos] = tight_subplot(1, 4, 0.01, 0.1,  0.05);
%set(ha(2:4), 'ycolor', 'none');

for i = 1:nresp
    axes(ha(i));
    %imagesc(sqrt(abs(data(:,:,zz,i))));
    imagesc(squeeze(data(256/2, :, :, i))');
    hold on
    colormap('gray');
    x1 = 1; x2 = 256;
    y1 = 52; y2 = 52;
    plot([x1,x2],[y1,y2], 'w--', 'Linewidth', 1);
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
end

sgtitle('Esophagus in respiratory phases (256x80)');
