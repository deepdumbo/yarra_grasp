function [x] = CSL1NlCg_4DDCE(x0,param)
% motion-weithed nonlinear CG slover

% starting point
x=x0;

% line search parameters
maxlsiter = 6;
gradToll = 1e-8 ;
param.l1Smooth = 1e-15;	
alpha = 0.01;  
beta = 0.6;
t0 = 1 ; 
k = 0;
% clear test
% compute g0  = grad(f(x))
g0 = grad(x,param);
dx = -g0;

% iterations
while(1)

    % backtracking line-search
	f0 = objective(x,dx,0,param);
	t = t0;
    f1 = objective(x,dx,t,param);
	lsiter = 0;
	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:))).^2 & (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		f1 = objective(x,dx,t,param);
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2, t0 = t0 * beta;end 
	if lsiter<1, t0 = t0 / beta; end
	x = (x + t*dx);

    % print some numbers for debug purposes	
    disp(sprintf('%d   , obj: %f, L-S: %d', k,f1,lsiter));
    k = k + 1;
    
    % stopping criteria (to be improved)
	if (k > param.nite) || (norm(dx(:)) < gradToll), break;end

    %conjugate gradient calculation
	g1 = grad(x,param);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	
end
return;

function res = objective(x,dx,t,param) %**********************************

% L2-norm part
w=(param.E*(x+t*dx)-param.y).*param.SG;
L2Obj=w(:)'*w(:);

% TV part along time
if param.TVWeight
    w = param.TV*(x+t*dx); 
    TVObj = sum((w(:).*conj(w(:))+param.l1Smooth).^(1/2));
else
    TVObj = 0;
end

res=L2Obj+param.TVWeight*TVObj;

function g = grad(x,param)%***********************************************

% L2-norm part
L2Grad = 2.*(param.E'*((param.E*x-param.y).*param.SG));

% TV part along time
if param.TVWeight
    w = param.TV*x;
    TVGrad = param.TV'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
else
    TVGrad=0;
end

g=L2Grad+param.TVWeight*TVGrad;
