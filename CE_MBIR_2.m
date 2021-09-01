%% REFERENCE
% https://en.wikipedia.org/wiki/Algebraic_reconstruction_technique

%% ART Equation
% x^(k+1) = x^k + lambda * AT(b - A(x))/ATA

%%
function x  = CE_MBIR_2(A,AT,b,x0 ,G , denoiser,lambda,alpha,sigma2 , rho ,tol,niter,bpos)



if (nargin < 13)
    bpos	= true;
end

if (nargin < 12)
    niter       = 1e2;
end

ATA    = cell(1,G);
v      = cell(1,G);
u      = cell(1,G);
w      = cell(1,G);
ww     = cell(1,G);
u_max  = ones(1,G);

umax   = 1;

parfor i = 1 : G

ATA{i}	= AT{i}(A{i}(ones(size(x0), 'single')));
v{i}    = x0;
u{i}    = (1e-2) * ones(size(x0));

u_max(i) = max(u{i},[],'all');

end 
iter = 0;
x = x0;
while umax > tol && iter <niter
    iter = iter +1;
    
    xx = zeros(size(x0)); % temp value for next iteration of x

    parfor i = 1 : G
    
        v{i}  	= v{i} + alpha*AT{i}(b{i} - A{i}(v{i}))./ATA{i} - alpha*imgradient(v{i} - denoiser(v{i},lambda))/G - alpha*(v{i}-x-u{i})./sigma2;
        
        if (bpos)
            v{i}(v{i} < 0) = 0;
        end
        
        w{i}    = v{i};
        ww{i}   = w{i};
        w{i}    = (2 * v{i} - x - u{i});
        w{i}    = rho * w{i} + (1 - rho) * ww{i}; 
        
        xx = xx + w{i};
    
    end
    
   x = xx./G;
  
    parfor i = 1 : G
        
        u{i} = x - w{i}; 
        u_max(i) = max(u{i},[],'all');
            
    end
    
    umax = max(u_max,[],'all');
    
    figure(1); colormap gray;
    imagesc(x);
    axis image off;
    title(num2str([umax, iter], 'x - v(i) = %d, iteration: %d'));
    drawnow();
end

x   = gather(x);

end