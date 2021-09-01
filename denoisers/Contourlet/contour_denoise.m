function out = contour_denoise(y,sigma,ref)
if ~exist('ref','var')
    ref = y;
    ref_flag = 0;
else
    ref_flag = 1;
end

y     = y*255;
sigma = sigma*255;
ref   = ref*255;

Pyr_mode = 1; 
eval(['load SDmm_' num2str(size(y,1)) '_' num2str(Pyr_mode)]);

dfilt = 'pkva';
nlev_SD = [2 3 4];
smooth_func = @rcos;
Pyr_mode = 1; 

% Take the ContourletSD transform
Y = ContourletSDDec(y, nlev_SD, Pyr_mode, smooth_func, dfilt);
R = ContourletSDDec(ref, nlev_SD, Pyr_mode, smooth_func, dfilt);

% The redundancy of the transform can be verified as follows.
% dstr1 = whos('Y');
% dstr2 = whos('Xn');
% dstr1.bytes / dstr2.bytes

% Apply hard thresholding on coefficients
Yth = Y;
for m = 2:length(Y)
  thresh = 3*sigma + sigma * (m == length(Y));
  for k = 1:length(Y{m})
      if ref_flag==0
        Yth{m}{k} = Y{m}{k}.* (abs(Y{m}{k}) > thresh*E{m}{k});
      else
        Yth{m}{k} = Y{m}{k}.* ((R{m}{k}.^2)./(R{m}{k}.^2 + sigma^2));
      end
  end
end

% ContourletSD reconstruction
Xd = ContourletSDRec(Yth, Pyr_mode, smooth_func, dfilt);
out = Xd/255;