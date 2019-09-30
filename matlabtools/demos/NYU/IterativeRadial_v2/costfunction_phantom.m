% ## Cost function for phantom-reconstruction demonstration
%
% Input:
%   x     = current image estimate (in vector form)
%   param = reconstruction parameters
%
% Output:
%   f     = Cost function value
%   g     = Gradient of cost function

function [f,g]=costfunction_phantom(x,param)

    global iterationCounter

    % Reshape the given estimate vector into 2D image form
    est_x=vec_to_img(x,param.br);
        
    % Show the current image estiamte
    imagescn(real(est_x),[],[],[],3);
    iterationCounter=iterationCounter+1;
    set(gcf, 'name', sprintf('Iterations -- (Evalutions: %d)',iterationCounter));
    drawnow;  

    % Calculate the data-fidelity term
    est_y=param.FT*est_x;
    res=est_y-param.y;  
    f=0.5*norm(res(:))^2;

    % Use preconditioning to accelerate convergence (needed for most
    % optimizers)
    res=param.dcf.*res;
    work=param.FT'*res;
    g=img_to_vec(work);

    % Calculate the penalty term
    if (param.enablePenalty==1)
        % Use helper function to calculate value and gradient of penalty
        [penaltyVal,penaltyGrad]=TV(est_x);

        % Print values for information purpose
        fprintf('Data fidelity=%f  Penalty=%f\n',f,param.lambdaPhantom*penaltyVal);

        % Add to cost function
        f=f + param.lambdaPhantom*penaltyVal;
        g=g + param.lambdaPhantom*img_to_vec(penaltyGrad);  
    end    
end

