function [water,fat,fieldmap,optinfo] = wfcs(kdata, x0, par)

%% Cast to specified precision
x0      = cast(x0,par.prec);
kdata   = cast(kdata,par.prec);

%% Set initial values
x          = x0;
dx         = complex(zeros(size(x),par.prec));

obj         = zeros(par.iter_gn,1);
f1          = zeros(par.iter_gn*par.iter_cg*par.maxIter_cg,1);
normdx      = zeros(par.iter_gn*par.iter_cg*par.maxIter_cg,1);
normdx0     = 0;

%% Optimization
for io = 1:par.iter_gn
    fprintf('\t Outer iteration: %d \n',io)
    
    res = kdata - forward_operator(x, par.coilmaps, par.TE, par.fatModulation, par.FT, par.Nread, par.Nproj, par.Ncoil, par.prec);
        
    par.data      = res;
    par.x         = x;
    par.DA        = DAMP(x, par.TE, par.fatModulation, par.FT, par.coilmaps, par.imsize, par.Nread, par.Nproj, par.Ncoil, par.prec);
    
    % CG iterations
    for ii = 1:par.iter_cg
        if par.verbose
            fprintf('\t \t Inner iteration: %d \n',ii)
        end
        [dx,f1_tmp,normdx_tmp,normdx0,breakflag] = nonlin_cg_gn(dx,par,ii,normdx0);
        
        if breakflag
            break;
        end
        
        % Save convergence information
        curridx = (io-1)*par.maxIter_cg*par.iter_cg + (ii-1)*par.maxIter_cg + 1;
        f1(curridx : curridx+par.maxIter_cg-1)          = f1_tmp;
        normdx(curridx : curridx+par.maxIter_cg-1)      = normdx_tmp;
        
    end

    obj(io) = norm(res(:));
    
    % Update image
    x = x + dx;
    
    % Reduce parameters
    par.alpha_n = par.alpha_n * par.q;
    
    if par.verbose
        figure(1);
        set(gcf, 'name', sprintf('Performing optimization'));
        imshow([abs(x(:,:,1)), abs(x(:,:,2))],[],'Border','tight');
        drawnow;
    end
    
end

water       = x(:,:,1);
fat         = x(:,:,2);
fieldmap    = x(:,:,3);

optinfo.obj         = obj;
optinfo.f1          = f1;
optinfo.normdx      = normdx;