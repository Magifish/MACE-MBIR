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
function [x,x_stack]  = GD_MBIR_3d(A,AT,b,x0 ,G , denoiser,lambda,alpha,niter,bpos)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GD_MBIR
%       [x,x_stack] = GD_MBIR_3d(A,AT,b,x0 ,G , denoiser,lambda, alpha,niter,bpos)
% This function implements GD_MBIR algorithm of Algorithm. 2 for 3D
% reconstruction.
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


if (nargin < 10)
    bpos	= true;
end

if (nargin <9)
    niter       = 1e2;
end

ATA    = cell(1,G);
v      = cell(1,G);
% u      = cell(1,G);
% w      = cell(1,G);
% ww     = cell(1,G);

for i = 1 : G

ATA{i}	= AT{i}(A{i}(ones(size(x0), 'single')));        % ATA operator
v{i}    = x0;
% u{i}    = 0.01*ones(size(x0));

end 

x = x0;

x_size = size(x0);

x_stack = zeros(x_size(1),x_size(2),x_size(3),niter);

for j = 1:niter
    
xx = zeros(size(x0)); % Temp value for next iteration of x

    
    for i = 1 : G
        
%         isa(v{i},'single')
%         class(v{i})
%         isa(single(v{i}),'single')
%         isa(b{i} - A{i}(single(v{i})),'single')
        
        v{i}  	= v{i} + alpha*AT{i}(b{i} - A{i}(v{i}))./ATA{i} - alpha*imgradient3(v{i} - denoiser(v{i},'TV',1,lambda))/G;  % Update the LR reconstruction v_i
        
        if (bpos)
            v{i}(v{i} < 0) = 0;
        end
        
        x = x + v{i};
    
    end
    
   x = x./G;
    
    for i = 1 : G
        
        v{i} = x; 
      
    end
    
    
    x_stack(:,:,:,j) = x(:,:,:);
    x_slide(:,:) = x(:,round(end/2),:);
    
    figure(1); colormap gray;
    imagesc(x_slide);
    axis image off;
    title(num2str([j, niter], '%d / %d'));
    drawnow();
end

x   = gather(x);

end