%% REFERENCE
% Consensus Equilibrium for Super-Resolution and Extreme-Scale CT SC %19, November 17–22, 2019, Denver, CO, USA

%% GD-MBIR Equation
% Algorithm 1 Parallel Mini-Batch Gradient Descent Algorithm
% INPUT: Uniformly partitioned measurements y and matrix A.
% LOCAL: G: the number of subsets. α: A fixed learning rate. ∇f i :
% A gradient of f i . vi : the ith LR reconstruction.
% OUTPUT: x: the consensus HR solution.
% 1: Initialize x, vi  
% 2: while Not Converged do
% 3: for each measurement subset i from 1 to G do in parallel
% 4: vi ←vi − α∇f i (vi )
% 5: end for
% 6: x ← Average(v1, · · · ,vG)
% 7: for each measurement subset i from 1 to G do in parallel
% 8: vi ←x
% 9: end for
% 10: end while

%%
function [x,x_stack] = GD_MBIR(A,AT,b,x0 ,G , denoiser,lambda, alpha,niter,bpos)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GD_MBIR
%       [x,x_stack] = GD_MBIR(A,AT,b,x0 ,G , denoiser,lambda, alpha,niter,bpos)
% This function implements GD_MBIR algorithm of Algorithm. 2
%
%   Input:
%       A       : system matrix operator A 
%       AT      : transpose of system matrix operator
%       b       : observation of projection data
%       x0      : initial state of image
%       G       : number of subsets
%       denoiser: denoiser operator
%       lambda  : denoiser parameter
%       alpha   : Gradent descent step size
%       niter   : number of iterations
%       bpos    : non-negetiveness of image values
%   Output:
%       x       : recontructed image
%       x_stack : stack of recontructed images across iterations
%
% YU CAO, July 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (nargin < 8)
    bpos	= true;
end

if (nargin < 7)
    niter       = 1e2;
end

ATA    = cell(1,G);
v    = cell(1,G);

for i = 1 : G

ATA{i}	= AT{i}(A{i}(ones(size(x0), 'single')));
v{i}    = x0;

end 

x_size = size(x0);

x_stack = zeros(x_size(1),x_size(2),niter);

for j = 1:niter
    
x = zeros(size(x0));

    for i = 1 : G
    
        v{i}  	= v{i} + alpha*AT{i}(b{i} - A{i}(v{i}))./ATA{i} - alpha*imgradient(v{i} - denoiser(v{i},lambda))/G;     % Update the LR reconstruction v_i
        
        if (bpos)
            v{i}(v{i} < 0) = 0;
        end
        
        x = x + v{i};
    
    end
    
   x = x./G;            % Averaging
    
    for i = 1 : G
        
        v{i} = x; 
      
    end
    
    x_stack(:,:,j) = x(:,:);
    figure(1); colormap gray;       % Instant plot
    imagesc(x);
    axis image off;
    title(num2str([j, niter], '%d / %d'));
    drawnow();
end

x   = gather(x);            % Turn back to gpuarray

end