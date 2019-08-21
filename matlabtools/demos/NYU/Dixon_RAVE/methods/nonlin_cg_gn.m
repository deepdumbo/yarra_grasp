function [x,f1save,normdx,normdx0,breakflag] = nonlin_cg_gn(x0,par,ii,normdx0)

% Set backtracking line search parameters
alpha    = par.alpha;
beta     = par.beta;
t0       = 1;
x        = x0;

[g0,DAx] = gradient(x,par);
dx       = -g0;
gtg      = g0(:)'*g0(:);
k        = 0;

% Update normdx0 for each GN step
if ii == 1
    normdx0   = norm(dx(:));
end

f1save      = zeros(par.maxIter_cg,1);
normdx      = zeros(par.maxIter_cg,1);
breakflag   = 0;

% Iterations
while (k < par.maxIter_cg)
    
    [DAdx, FD1w, FD1dw, FD1f, FD1df] = preobjective(x,dx,par);
    f0 = objective(x, dx, DAx, DAdx, FD1w, FD1dw, FD1f, FD1df, 0, par);

    t       = t0;
    iter_ls = 0;
    f1      = objective(x, dx, DAx, DAdx, FD1w, FD1dw, FD1f, FD1df, t, par);
    
    % Backtraching line search
    c2 = -abs(g0(:).'*dx(:));
    while((f1 > f0 + alpha*t/c2)^2 && (iter_ls < par.maxIter_ls))                           %#ok
        iter_ls  = iter_ls + 1;
        t  = t * beta;
        f1 = objective(x, dx, DAx, DAdx, FD1w, FD1dw, FD1f, FD1df, t, par);
    end
    
    % Adjust initial step size for faster search
    if iter_ls > 2
        t0 = t0 * beta;
    end
    
    if iter_ls < 1
        t0 = t0 / beta;
    end
    
    % Update image
    x           = x + t*dx;
    
    % Update gradient
    [g1,DAx]    = gradient(x,par);
    
    % Update search direction
    gtg_new     = g1(:)'*g1(:);
    bk          = gtg_new/(gtg + eps);
    gtg         = gtg_new;
    g0          = g1;
    
    dx          = -g1 + bk*dx;
    k           = k + 1;
    
    f1save(k)   = f1;
    normdx(k)   = norm(dx(:));

    if (norm(dx(:))/normdx0 < par.threshold)
        breakflag = 1;
        break;
    end
    
    if par.verbose
        fprintf('\t \t \t Iteration k = %d \n',k)
        fprintf('\t \t \t iter_ls = %d \n',iter_ls)
        fprintf('\t \t \t normdx0 = %.4f \n',normdx0)
        fprintf('\t \t \t normdx = %.4f \n',normdx(k))
        fprintf('\t \t \t f1 = %.4f \n',f1)
        fprintf('\t \t \t ------------------- \n')        
    end
    
end


function [DAdx, FD1w, FD1dw, FD1f, FD1df] = preobjective(x,dx,par)

DAdx   = par.DA*dx;

if par.FD1Weight
    % water part
    FD1w  = par.FD1*(x(:,:,1) + par.x(:,:,1));
    FD1dw = par.FD1*dx(:,:,1);
    % fat part
    FD1f  = par.FD1*(x(:,:,2) + par.x(:,:,2));
    FD1df = par.FD1*dx(:,:,2);
else
    FD1w  = 0;
    FD1dw = 0;
    FD1f  = 0;
    FD1df = 0;
end

function res = objective(x, dx, DAx, DAdx, FD1w, FD1dw, FD1f, FD1df, t, par)

res = (DAx + t*DAdx - par.data);
res = res(:)'*res(:);

if par.FD1Weight
    p     = par.pnorm_im;
    temp1 = FD1w + t*FD1dw;
    temp2 = FD1f + t*FD1df;
    FD1   = (temp1.*conj(temp1) + par.mu).^(p/2) + (temp2.*conj(temp2) + par.mu).^(p/2);
else
    FD1   = 0;
end

tmp = x + t*dx;
tmp = tmp(:)'*tmp(:);

res = par.dataWeight*res + par.alpha_n*tmp +  par.FD1Weight*sum(FD1(:));


function [grad,DAx] = gradient(x,par)

grad(:,:,1) = par.FD1Weight*gFD1(x(:,:,1) + par.x(:,:,1), par);
grad(:,:,2) = par.FD1Weight*gFD1(x(:,:,2) + par.x(:,:,2), par);
grad(:,:,3) = 0*(x(:,:,3) + par.x(:,:,3));

[gradData,DAx] = gData(x,par);
grad = grad + gradData + par.alpha_n*2*x;

function [grad,DAx] = gData(x,par)
DAx = par.DA*x;
grad = par.dataWeight*2*(par.DA'*(DAx - par.data));

function grad = gFD1(x,par)
p    = par.pnorm_im;
fd   = par.FD1*x;
grad = p*fd.*(fd.*conj(fd) + par.mu).^(p/2-1);
grad = par.FD1'*grad;
