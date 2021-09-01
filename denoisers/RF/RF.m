function xbaseline = RF(y,sigma,ref)
if ~exist('ref','var')
    sigma_s     = 3;
    sigma_r     = sigma;
    noise_sigma = sigma;
    max_itr     = 1;
    sigma_g     = 2;
    ref  = RF_1st_mex_opt(y, sigma_s, sigma_r, noise_sigma, sigma_g, max_itr);
    
    sigma_s     = 6;
    sigma_r     = 1*sigma;
    noise_sigma = 0.02*sigma;
    max_itr     = 1;
    sigma_g     = 0.5;
    xbaseline  = RF_3rd_mex_opt(y, sigma_s, sigma_r, noise_sigma, sigma_g, max_itr, ref);    
else
    sigma_s     = 6;
    sigma_r     = 1*sigma;
    noise_sigma = 0.0*sigma;
    max_itr     = 1;
    sigma_g     = 0.1;
    xbaseline  = RF_3rd_mex_opt(y, sigma_s, sigma_r, noise_sigma, sigma_g, max_itr, ref);
end
end
