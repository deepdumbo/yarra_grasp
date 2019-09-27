%script that creates DC-centered spectra of measurements of 
%breathing pattern as sampled with GRASP.

R_os = 1.083; %oversampling of 8.3 per cent
R_pf = 5.0/8; %partial fourier filling of 3/8 of total slices
n_par = 48; %number of measured partitions in experiment w/o oversampling

%input repetition times
TRs = 1.e-3*[4.3 9.3 15 22.5 30]; 
n_exp = numel(TRs);

rs_tr4 = load('rs_tr4');
rs_tr4 = rs_tr4.Res_Signal; 
rs_tr9 = load('rs_tr9');
rs_tr9 = rs_tr9.Res_Signal; 
rs_tr15 = load('rs_tr15');
rs_tr15 = rs_tr15.Res_Signal; 
rs_tr22_5 = load('rs_tr22_5');
rs_tr22_5 = rs_tr22_5.Res_Signal; 
rs_tr30 = load('rs_tr30');
rs_tr30 = rs_tr30.Res_Signal; 

%save time series to matrix for loop
RSs = zeros(512, n_exp);
RSs(:,1) = rs_tr4;
RSs(:,2) = rs_tr9;
RSs(:,3) = rs_tr15;
RSs(:,4) = rs_tr22_5;
RSs(:,5) = rs_tr30;

%L contains all time sequences (for every TR)
L = zeros(512, n_exp);
for i = 1:n_exp
    L(:,i) = (0:511)*TRs(i)*n_par*R_os*R_pf;
end

%for i = 1:n_exp
%    figure();
%    txt = ['TR (ms) = %d', (TRs(i)*1.e3)];
%    disp(txt)
%    plot(L(:,i), RSs(:,i));
%end


M = zeros(n_exp);
M(1,:) = TRs; % Repetition times in secs
M(2,:) = TRs*n_par*R_os*R_pf; % Dwell times in secs
M(3,:) = 0.5./M(2,:); % Bandwidths in Hz

%frequency ranges
N = zeros(512, n_exp);
for i = 1:n_exp
    N(:,i) = (-512/2:512/2-1).*(M(3,i)./512);
end
disp(N);

%FFTs of time series
Q = zeros(512, n_exp);
for i = 1:n_exp
    Q(:,i) = abs(fftshift(fft(RSs(:,i))));
end

%plotting
figure();
hold on;
for i = 1:n_exp
    plot(N(:,i),Q(:,i));
end
axis([0 2.5 0 100]);
xlabel('$f$ (Hz)', 'Interpreter', 'latex');
title('Fourier Transforms of respiratory signals', 'Interpreter', 'latex');
legend('TR = 4.3 ms', 'TR = 9.3 ms', 'TR = 15 ms', ...
    'TR = 22.5 ms', 'TR = 30 ms');
