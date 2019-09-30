% ## Cost function for estimating the joint-coil image
%
% Input:
%   x     = current image estimate (in vector form)
%   param = reconstruction parameters
%
% Output:
%   f     = Cost function value
%   g     = Gradient of cost function

function [f,g]=costfunction_image(x,param)

    global iterationCounter

    % Reshape the given estimate vector into 2D image form
    est_x=vec_to_img(x,param.br);

    % Show the current image estiamte
    imagescn(real(est_x),[],[],[],3);
    iterationCounter=iterationCounter+1;
    set(gcf, 'name', sprintf('Step 2 -- Joint-Coil Iterations -- (Evalutions: %d)',iterationCounter));
    drawnow;  

    % ## Calculate the data fidelity term
    g=zeros(param.nvar,1);
    f=0;
    
    % Loop over the channels and match with data from all coils
    for ic=1:param.channels      
        % Multiply estimate with coil profile before NUFFT
        est_y=param.FT*(param.coilprofile(:,:,ic).*est_x);
        res=est_y-squeeze(param.y(:,ic,:,:));
        f=f+0.5*norm(res(:))^2;

        % Use preconditioning to accelerate convergence (needed for most
        % optimizers)
        res=param.dcf.*res;
        work=param.FT'*res;
        work=conj(param.coilprofile(:,:,ic)).*work;

        g=g+img_to_vec(work);
    end 

    % ## Calculate the penalty term    
    if (param.enablePenalty==1)
        % Use helper function to calculate value and gradient of penalty
        [penaltyVal,penaltyGrad]=TV(est_x);

        % Add to cost function
        fprintf('Data fidelity=%f  Penalty=%f\n',f, param.lambdaImage*penaltyVal);
        f=f + param.lambdaImage*penaltyVal;
        g=g + param.lambdaImage*img_to_vec(penaltyGrad);  
    end
end

