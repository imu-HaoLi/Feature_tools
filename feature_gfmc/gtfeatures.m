function [gfcc, g] = gtfeatures(input_signal)
% Generate gammatone features (GF) and gammatone frequency cepstral
% coefficients (GFCC).
% Written by Yang Shao, and adapted by Xiaojia Zhao in Oct'11

sampFreq = 16e3;

numChannel = 64;

gt = gen_gammaton(sampFreq, numChannel);  % get gammatone filterbank
    
sig = input_signal;

sig = reshape(sig, 1, length(sig));  % gammatone filtering function requires input to be row vector

g=fgammaton(sig, gt, sampFreq, numChannel);     % gammatone filter pass and decimation to get GF features

gfcc = gtf2gtfcc(g, 1, 13);  % apply dct to get GFCC features with 0th coefficient removed

%writeHTK(strcat('entry','.gtf'), g);     % write GF features
%writeHTK(strcat('entry','.gtfcc'), gfcc);     % write GFCC features
end
