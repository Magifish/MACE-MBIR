function out = wrapper_BM3D(in,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out = wrapper_BM3D(in,sigma)
% performs BM3D denoising
% 
% Require BM3D package
%
% Download:
% http://www.cs.tut.fi/~foi/GCF-BM3D/
%
% Xiran Wang and Stanley Chan
% Copyright 2016
% Purdue University, West Lafayette, In, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = in;
[~,zhat] = BM3D(1,y,sigma);
out = zhat;

end