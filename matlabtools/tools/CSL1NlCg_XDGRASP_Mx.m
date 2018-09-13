function [x] = CSL1NlCg_XDGRASP_Mx(x0,param)

% Non-linear Conjugate Gradient Algorithm adapted from the code provided by
% Miki Lustig in http://www.eecs.berkeley.edu/~mlustig/Software.html
%
% Input:
% (1) x0, starting point (gridding images)
% (2) param, reconstruction parameters
%
% Output:
% (1) x, reconstructed images
%
% (c) Li Feng, 2016, Stanford University


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
% 	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
    %MCM changed this: the ^2 operation appears wrongly placed here, namely
    %on a logical outcome. It looks like it should be over
    %alpha*t*abs(...). But the last term is itself already a (type of)
    %square operation, at least in the first iteration, where g0 = dx and
    %the multiplication of transpose by original matrix would yield a sum
    %of squares (and then why take an absolute value...)
    %Second, the square operation crashes if any of the inner components is
    %on the GPU. 
    %So we could either take out the square altogether, or move it inside
    %the parentheses.
    %Checking the workings of ConjGrad algorithms (more specifically the 
    %backtracking line search [or Armijo rule]  part of it), I find that 
    %indeed the square operation should not be there. 
    %Note that the absolute value operation is there because the inner product
    %is complex valued.
    %value doesn't belong there. Let's get rid of that too...
    while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:))) && (lsiter<maxlsiter) 
		lsiter = lsiter + 1;
		t = t * beta;
		f1 = objective(x,dx,t,param);
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2, t0 = t0 * beta;end 
	if lsiter<1, t0 = t0 / beta; end
%     test(:,:,k+1)=x(:,:,103,9);
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
w=param.E*(x+t*dx)-param.y;
L2Obj=w(:)'*w(:);

% TV part along time
if param.TVWeight_dim1
    w = param.TV_dim1*(x+t*dx); 
    TVObj_dim1 = sum((w(:).*conj(w(:))+param.l1Smooth).^(1/2));
else
    TVObj_dim1 = 0;
end

% TV part along respiration
if param.TVWeight_dim2
    w = param.TV_dim2*(x+t*dx); 
    TVObj_dim2 = sum((w(:).*conj(w(:))+param.l1Smooth).^(1/2));
else
    TVObj_dim2 = 0;
end

res=L2Obj+param.TVWeight_dim1*TVObj_dim1+param.TVWeight_dim2*TVObj_dim2;

function g = grad(x,param)%***********************************************

% L2-norm part
L2Grad = 2.*(param.E'*(param.E*x-param.y));

% TV part along time
if param.TVWeight_dim1
    w = param.TV_dim1*x;
    TVGrad_dim1 = param.TV_dim1'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
else
    TVGrad_dim1=0;
end

% TV part along respiration
if param.TVWeight_dim2
    w = param.TV_dim2*x;
    TVGrad_dim2 = param.TV_dim2'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
else
    TVGrad_dim2=0;
end

g=L2Grad+param.TVWeight_dim1*TVGrad_dim1+param.TVWeight_dim2*TVGrad_dim2;
