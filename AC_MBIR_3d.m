%% REFERENCE
% Consensus Equilibrium for Super-Resolution and Extreme-Scale CT SC %19, November 17–22, 2019, Denver, CO, USA


%% AC-MBIR Equation
% Algorithm 4 AC-MBIR
% INPUT: The partitioned CT measurements and system matrix.
% LOCAL: G: the number of subsets. vi : the ith LR reconstruction.
%        [vi ]j : the jth sub-volume of vi . [vi ]j,k : the kth SV of [vi ]j .
%        Similar notations apply to wi , ui .
% OUTPUT: x: the consensus HR solution.
% 1: Initialize x, vi , wi , and ui
% 2: while x , vi , i = 1, 2, ...,G do
% 3:    for each subset i do in parallel across clusters
% 4:        for each sub-volume [vi ]j ∈ vi , do in parallel across nodes
% 5:            for each SV [vi ]j,k ∈ [vi ]j , do in parallel across cores
% 6:                [vi ]j,k ← argmin[vi ]j,k Fi ([x]j,k + [ui ]j,k )
% 7:                [(w′)i ]j,k ← [wi ]j,k
% 8:                [wi ]j,k ← (2[vi ]j,k − [x]j,k − [ui ]j,k )
% 9:                [wi ]j,k ← ρ[wi ]j,k +􀀀1 − ρ�[(w′)i ]j,k
% 10:               [x]j,k ← Average([w1]j,k , · · · , [wG]j,k ), through MPI all reduce among nodes in a cabinet
% 11:               [ui ]j,k ← [x]j,k − [wi ]j,k
% 12:           end for
% 13:       end for
% 14:       Neighboring sub-volumes exchange boundary voxels through MPI send/receive operations.
% 15:   end for
% 16: end while

%%
function [x_c,x_stack]  = AC_MBIR_3d(A,AT,b,x0 ,geo_s,THETA,G , J ,K , denoiser,lambda,alpha,sigma2 , rho ,niter,bpos)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AC_MBIR
%       [x_c,x_stack] = AC_MBIR_3d(A,AT,b,x0 ,geo_s,THETA,G , J ,K , denoiser,lambda,alpha,sigma2 , rho ,niter,bpos)
% This function implements CE_MBIR algorithm of Algorithm. 4 for 3D
% reconstruction.
%
%   Input:
%       A       : system matrix operator A 
%       AT      : transpose of system matrix operator
%       b       : observation of projection data
%       x0      : initial state of image
%       geo_s   : geometry information for each super-voxel
%       THETA   : angle information for each super-voxel
%       G       : number of subsets
%       J       : number of sub-volumes
%       K       : number of super-voxels
%       denoiser: denoiser operator
%       lambda  : denoiser parameter
%       alpha   : Gradent descent step size
%       sigma2  : proximal function scaling factor
%       rho     : Mann iteration parameter
%       niter   : number of iterations
%       bpos    : non-negetiveness of image values
%   Output:
%       x_c     : recontructed image
%       x_stack : stack of recontructed images across iterations
%
% YU CAO, July 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin < 16)
    bpos	= true;
end

if (nargin < 15)
    niter       = 1e2;
end

x_size = size(x0);
x_c = zeros(x_size);
% b_size = size(b{1});

x0_     = cell(J,K);
x       = cell(J,K);

xx      = cell(J,K);
x_slide = cell(J,K);
% x0_{1,1}     = x0(end/J,end/K,:);

v      = cell(G,J,K);
u      = cell(G,J,K);
w      = cell(G,J,K);
ww     = cell(G,J,K);
% b_     = cell(G,J,K);
ATA    = cell(G,J,K); 

for j = 1 : J
        for k = 1 : K
%             x0_{j,k} = zeros(x_size);
            x0_{j,k}= x0((x_size(1)/J*(j-1)+1:x_size(1)/J*j),(x_size(2)/K*(k-1)+1:x_size(2)/K*k),:);        % Reallocate the initial state
        end
end

x = x0;

x_size = size(x0);

x_stack = zeros(x_size(1),x_size(2),x_size(3),niter);

for i = 1 : G
    


    for j = 1 : J
        for k = 1 : K
            ATA{i,j,k}	= AT(A(ones(size(x0_{1,1}), 'single'),geo_s{j,k},THETA(i,:)),geo_s{j,k},THETA(i,:));        % ATA operator
            v{i,j,k}    = x0_{j,k};
            u{i,j,k}    = 0.001*ones(size(x0_{1,1}));           % Initialize update difference
%             b_{i,j,k} = zeros(b_size);
%             b_{i,j,k}(b_size(1)/J*(j-1)+1:b_size(1)/J*j,b_size(2)/K*(k-1)+1:b_size(2)/K*k,:) = b{i}(b_size(1)/J*(j-1)+1:b_size(1)/J*j,b_size(2)/K*(k-1)+1:b_size(2)/K*k,:);
        end
    end
    



end 

    for j = 1 : J
        for k = 1 : K
            x{j,k} = x0_{j,k};          % Reallocate the initial state
        end
    end

%     size(v{1,1,1})
%     size(b_{1,1,1})
%     size(A(v{1,1,1},geo_s{1,1},THETA(1,:)))
%     size(ATA{1,1,1})

for iter = 1:niter

     for j = 1 : J
        for k = 1 : K
            xx{j,k} = zeros(size(x0_{1,1})); % temp value for next iteration of x
        end
     end 
     
    figure(1);colormap gray;
    title(num2str([iter, niter], '%d / %d'));
    
    for j = 1 : J
        for k = 1 : K
    
    x_slide{j,k}(:,:) = x{j,k}(:,round(end/2),:);
    
    subplot(J,K,2*(j-1)+k);
    imagesc(x_slide{j,k});
    axis image off;
    drawnow();      
        end
    end

    
 
        
   for j = 1 : J
        for k = 1 : K
%         isa(v{i},'single')
%         class(v{i})
%         isa(single(v{i}),'single')
%         isa(b{i} - A{i}(single(v{i})),'single')
        
                for i = 1 : G
%         volumeViewer(v{i,j,k});
        % 5 iterations of gradient descent
        v{i,j,k} = v{i,j,k} + alpha*AT((b{i,j,k} - A(v{i,j,k},geo_s{j,k},THETA(i,:))),geo_s{j,k},THETA(i,:))./ATA{i,j,k} - alpha*imgradient3(v{i,j,k} - denoiser(v{i,j,k},'TV',1,lambda))/G - alpha*(v{i,j,k}-x{j,k}-u{i,j,k})./sigma2;
        v{i,j,k} = v{i,j,k} + alpha*AT((b{i,j,k} - A(v{i,j,k},geo_s{j,k},THETA(i,:))),geo_s{j,k},THETA(i,:))./ATA{i,j,k} - alpha*imgradient3(v{i,j,k} - denoiser(v{i,j,k},'TV',1,lambda))/G - alpha*(v{i,j,k}-x{j,k}-u{i,j,k})./sigma2;
        v{i,j,k} = v{i,j,k} + alpha*AT((b{i,j,k} - A(v{i,j,k},geo_s{j,k},THETA(i,:))),geo_s{j,k},THETA(i,:))./ATA{i,j,k} - alpha*imgradient3(v{i,j,k} - denoiser(v{i,j,k},'TV',1,lambda))/G - alpha*(v{i,j,k}-x{j,k}-u{i,j,k})./sigma2;
        v{i,j,k} = v{i,j,k} + alpha*AT((b{i,j,k} - A(v{i,j,k},geo_s{j,k},THETA(i,:))),geo_s{j,k},THETA(i,:))./ATA{i,j,k} - alpha*imgradient3(v{i,j,k} - denoiser(v{i,j,k},'TV',1,lambda))/G - alpha*(v{i,j,k}-x{j,k}-u{i,j,k})./sigma2;
        v{i,j,k} = v{i,j,k} + alpha*AT((b{i,j,k} - A(v{i,j,k},geo_s{j,k},THETA(i,:))),geo_s{j,k},THETA(i,:))./ATA{i,j,k} - alpha*imgradient3(v{i,j,k} - denoiser(v{i,j,k},'TV',1,lambda))/G - alpha*(v{i,j,k}-x{j,k}-u{i,j,k})./sigma2;
        %         volumeViewer(v{i,j,k});
        
        
        if (bpos)
            v{i,j,k}(v{i,j,k} < 0) = 0;
        end
        
%         volumeViewer(v{i,j,k});
        
        w{i,j,k}    = v{i,j,k};
        ww{i,j,k}   = w{i,j,k};
        w{i,j,k}    = (2 * v{i,j,k} - x{j,k} - u{i,j,k});
        w{i,j,k}    = rho * w{i,j,k} + (1 - rho) * ww{i,j,k};       % Mann iteration update
        
        xx{j,k} = xx{j,k} + w{i,j,k};
                end
        x{j,k} = xx{j,k}./G;
        
%         volumeViewer(x{j,k});
        
               for i = 1 : G    
        u{i,j,k} = x{j,k} - w{i,j,k};       % Update difference
      
               end
         end
    end
    

    
%     for j = 1 : J
%         for k = 1 : K
%     
%     x_slide{j,k}(:,:) = x{j,k}(:,:,round(end/2));
%     
%     subplot(J,K,2*(j-1)+k);
%     imagesc(x_slide{j,k});
%     axis image off;
%     drawnow();      
%         end
%     end
    % Instant plot of super-voxels
    for j = 1 : J
        for k = 1 : K
            
             x_c(x_size(1)/J*(j-1)+1:x_size(1)/J*j,x_size(2)/K*(k-1)+1:x_size(2)/K*k,:)=x{j,k}; %(x_size(1)/J*(j-1)+1:x_size(1)/J*j,x_size(2)/K*(k-1)+1:x_size(2)/K*k,:);
             
        end
    end
    x_stack(:,:,:,j) = x_c(:,:,:);
    x_c_slide(:,:) = x_c(:,round(end/2),:);
    
    figure(2);  colormap gray;
    imagesc(x_c_slide);         
%    axis image off;
%    title(num2str([iter, niter], '%d / %d'));
    drawnow();
end


    
x_c   = gather(x_c);

end