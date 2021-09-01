%% REFERENCE
% https://en.wikipedia.org/wiki/Algebraic_reconstruction_technique

%% ART Equation
% x^(k+1) = x^k + lambda * AT(b - A(x))/ATA

%%
function x  = GD_MBIR_2(A,AT,b,x0 ,G , denoiser,lambda, alpha,niter,bpos)



if (nargin < 8)
    bpos	= true;
end

if (nargin < 7)
    niter       = 1e2;
end

ATA    = cell(1,G);
% v    = cell(1,G);

v = repmat(x0,[1 1 G]);

for i = 1 : G

ATA{i}	= AT{i}(A{i}(ones(size(x0), 'single')));

end 



for j = 1:niter
    
x = zeros(size(x0));

    for i = 1 : G
    
        v(:,:,i) 	= v(:,:,i) + alpha*AT{i}(b(:,:,i) - A{i}(v(:,:,i)))./ATA{i} - alpha*imgradient(v(:,:,i) - denoiser(v(:,:,i),lambda))/G;
        
        if (bpos)
            v(v(:,:,i) < 0) = 0;
        end
        
        x = x + v(:,:,i);
    
    end
    
   x = x./G;
    
    for i = 1 : G
        
        v(:,:,i) = x; 
      
    end
    
    figure(1); colormap gray;
    imagesc(x);
    axis image off;
    title(num2str([j, niter], '%d / %d'));
    drawnow();
end

x   = gather(x);

end