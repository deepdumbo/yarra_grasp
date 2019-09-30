function plot_nresp4_necho3(data1, data2, data3, zz)

nresp = 4;
necho = 3;

% data12 = data1.data_echo1;
% data22 = data2.data_echo2;
% data32 = data3.data_echo3;
% 
% echo1 = data12(:, :, zz, :);
% echo2 = data22(:, :, zz, :);
% echo3 = data32(:, :, zz, :);

echo1 = data1(:, :, zz, :);
echo2 = data2(:, :, zz, :);
echo3 = data3(:, :, zz, :);

nresp = 4

figure('Position', [100 300 1200 900]);

[ha, pos] = tight_subplot(3, 4, 0.01, 0.1,  0.05);
%set(ha(2:4), 'ycolor', 'none');
disp(ha);

for i = 1:nresp
    axes(ha(i));
    imagesc(sqrt(abs(echo1(:,:,1,i))));
    colormap('gray');
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
end


for i = 1:nresp
    axes(ha(i+nresp));
    imagesc(sqrt(abs(echo2(:,:,1,i))));
    colormap('gray');
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
end


for i = 1:nresp
    axes(ha(i+2*nresp));
    imagesc(sqrt(abs(echo3(:,:,1,i))));
    colormap('gray');
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
end


end