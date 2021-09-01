%% REFERENCE
% Consensus Equilibrium for Super-Resolution and Extreme-Scale CT SC %19, November 17–22, 2019, Denver, CO, USA
% Uniform sparse-view

%% CE-MBIR Equation
% Algorithm 2 Consensus Equilibrium MBIR
% INPUT: Measurements y and matrix A partitioned in any order.
% LOCAL: wi : a temporary copy for vi . ui : input change to proximal
% map function, Fi (x). Other variables are the same as in
% Algorithm 1.
% OUTPUT: x: the consensus HR solution.
% 1: Initialize x, vi , wi , and ui
% 2: while x , vi , i = 1, 2, ...,G do
% 3: for each subset i from 1 to G do in parallel
% 4: vi ← argminvi Fi (x + ui )
% 5: (w′)i ←wi
% 6: wi ← (2vi − x − ui )
% 7: wi ← ρwi +􀀀1 − ρ�(w′)i , 0 ≤ ρ ≤ 1
% 8: end for
% 9: x ← Average(w1, · · · ,wG)
% 10: for each subset i from 1 to G do in parallel
% 11: ui ←x − wi
% 12: end for
% 13: end while
%% 

%%
function [x,x_stack]  = CE_MBIR(A,AT,b,x0 ,G , denoiser,lambda,alpha,sigma2 , rho ,niter,bpos)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CE_MBIR
%       [x,x_stack] = CE_MBIR_2d(A,AT,b,x0 ,G ,denoiser,lambda, alpha,niter,bpos)
% This function implements CE_MBIR algorithm of Algorithm. 3 for 2D
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
%       sigma2  : proximal function scaling factor
%       rho     : Mann iteration parameter
%       niter   : number of iterations
%       bpos    : non-negetiveness of image values
%   Output:
%       x       : recontructed image
%       x_stack : stack of recontructed images across iterations
%
% YU CAO, July 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin < 12)
    bpos	= true;
end

if (nargin < 11)
    niter       = 1e2;
end

ATA    = cell(1,G);
v      = cell(1,G);
u      = cell(1,G);
w      = cell(1,G);
ww     = cell(1,G);

for i = 1 : G

ATA{i}	= AT{i}(A{i}(ones(size(x0), 'single')));      % ATA operator
v{i}    = x0;
u{i}    = 0.0001*ones(size(x0));                      % Initialize update difference

end 

x = x0;

x_size = size(x0);

x_stack = zeros(x_size(1),x_size(2),niter);

for j = 1:niter
    
xx = zeros(size(x0)); % Temp value for next iteration of x

    for i = 1 : G
        % 5 iterations of gradient descent
        v{i}  	= v{i} + alpha*AT{i}(b{i} - A{i}(v{i}))./ATA{i} - alpha*imgradient(v{i} - denoiser(v{i},lambda))/G - alpha*(v{i}-x-u{i})./sigma2;
        v{i}  	= v{i} + alpha*AT{i}(b{i} - A{i}(v{i}))./ATA{i} - alpha*imgradient(v{i} - denoiser(v{i},lambda))/G - alpha*(v{i}-x-u{i})./sigma2;
        v{i}  	= v{i} + alpha*AT{i}(b{i} - A{i}(v{i}))./ATA{i} - alpha*imgradient(v{i} - denoiser(v{i},lambda))/G - alpha*(v{i}-x-u{i})./sigma2;
        v{i}  	= v{i} + alpha*AT{i}(b{i} - A{i}(v{i}))./ATA{i} - alpha*imgradient(v{i} - denoiser(v{i},lambda))/G - alpha*(v{i}-x-u{i})./sigma2;
        v{i}  	= v{i} + alpha*AT{i}(b{i} - A{i}(v{i}))./ATA{i} - alpha*imgradient(v{i} - denoiser(v{i},lambda))/G - alpha*(v{i}-x-u{i})./sigma2;
        
        if (bpos)
            v{i}(v{i} < 0) = 0;
        end
        
        w{i}    = v{i};
        ww{i}   = w{i};
        w{i}    = (2 * v{i} - x - u{i});
        w{i}    = rho * w{i} + (1 - rho) * ww{i};       % Mann iteration update
        
        xx = xx + w{i};
    
    end
    
   x = xx./G;
  
    for i = 1 : G
        
        u{i} = x - w{i};                % Update difference
      
    end
    
    x_stack(:,:,j) = x(:,:);
    figure(1); colormap gray;           % Instant plot
    imagesc(x);
    axis image off;
    title(num2str([j, niter], '%d / %d'));
    drawnow();
end

x   = gather(x);

end