%Motion detection: step 3

%Li Feng, NYU, 12/18/2017

function [Res_Signal]=MC_Step3(Res_Signal,ZIP1);

close all

warning off

ntviews=length(Res_Signal);
ft = fittype( 'smoothingspline' );
opts = fitoptions( ft );
opts.SmoothingParam = 0.015;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate the peak positions
t=6;
[Peak,Peak_Index]=findpeaks(double(Res_Signal),'MinPeakDistance',t);
figure,plot(Res_Signal);
hold on; plot(Peak_Index,Res_Signal(Peak_Index),'ro');

% check whether the peak positions are reasonable or not.
II=6;
kk=1;clear Peak_Index_new
for ii=1:length(Peak_Index)
    t1=Res_Signal(Peak_Index(ii))>Res_Signal(max(Peak_Index(ii)-II,1):Peak_Index(ii)-1);
    t2=Res_Signal(Peak_Index(ii))>Res_Signal(Peak_Index(ii)+1:min(Peak_Index(ii)+II,length(Res_Signal)));
    if isempty(find([t1;t2]==0))
        Peak_Index_new(kk)=Peak_Index(ii);
        kk=kk+1;
    end
end
Peak_Index=Peak_Index_new;
figure,plot(Res_Signal);
hold on; plot(Peak_Index,Res_Signal(Peak_Index),'ro');

%Do a fitting and substract the fitted signal
[xData, yData] = prepareCurveData(Peak_Index,Res_Signal(Peak_Index));
[fitresult, gof] = fit( xData, yData, ft, opts );

cfval = coeffvalues(fitresult);
ftmax = feval(ft,cfval(1),(1:ntviews)');
Res_Signal=Res_Signal-ftmax;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
figure,imagesc(abs(ZIP1(:,:,1))),axis image, axis off, colormap(gray),title('Respiratory Motion')
hold on
plot(-Res_Signal*150+70,'r')
Res_Signal=single(Res_Signal);