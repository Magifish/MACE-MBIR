function out = nrmse(x,y)
MSE = mean( (x-y).^2, 'all' );
out = sqrt(MSE)/(max(x,[],'all')-min(x,[],'all'));
end