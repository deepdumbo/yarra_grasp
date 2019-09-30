function plot_radial_kspace(kdata)

% Function for the density plot of radial k-space data. Input a single
% partition

if size(kdata) > 2
    error('Please input k-space data in 2 dimensions.')
end

[nx, ntviews] = size(kdata);
traj = Trajectory_GoldenAngle(ntviews, nx, 0, [0 0], 1);

kdata = reshape(kdata, [nx*ntviews, 1]);
traj = i*reshape(traj, [nx*ntviews, 1]);

kdata = abs(kdata); 
X = real(traj);
Y = imag(traj);

figure()
kmin=min(kdata);
kmax=max(kdata);
map=colormap;
color_steps=size(map,1);

hold on
for i=1:color_steps
    ind=find(kdata<kmin+i*(kmax-kmin)/color_steps & kdata>=kmin+(i-1)*(kmax-kmin)/color_steps);
    plot(X(ind),Y(ind),'.','Color',map(i,:));
end

end