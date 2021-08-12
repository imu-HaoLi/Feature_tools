function cepDpDD = gfmc(input_signal)
% gfmc, a gfcc variation with modulation information

sampFreq = 16e3;
numChannel = 64;

gt = gen_gammaton(sampFreq, numChannel);  % get gammatone filterbank
    
sig = input_signal;

sig = reshape(sig, 1, length(sig));  % gammatone filtering function requires input to be row vector

g=fgammaton(sig, gt, sampFreq, numChannel);     % gammatone filter pass and decimation to get GF features
        
cep = gtf2gtfcc(g,1,31);  % apply dct to get GFCC features without 0th coefficient removed
%cep = g;

%% compute modulation coefficients
winLength = 16;      % default window length in sample points which is 20 ms for 16 KHz sampling frequency
winShift = 1;  

numChan = size(cep,1);
sigLength = length(sig);     % number of channels and input signal length          
increment = winLength/winShift;    % special treatment for first increment-1 frames
M = floor(sigLength/160);     % number of time frames

cep = cep(:,1:M);

% calculate energy for each frame in each channel
modulation = zeros(numChan,M);

% 2-16 HZ energy
start_fft = ceil(2/(100/256));
end_fft = floor(16/(100/256));

for m = 1:M
   for i = 1:numChan
       if m < increment      % shorter frame lengths for beginning frames
           fft_tmp = fft(cep(i,1:m*winShift),256);  
           fft_tmp = fft_tmp(start_fft:end_fft);
           modulation(i,m) = fft_tmp*fft_tmp';
       else
           startpoint = (m-increment)*winShift;
           fft_tmp = fft(cep(i,startpoint+1:startpoint+winLength) , 256);
           fft_tmp = fft_tmp(start_fft:end_fft);
           modulation(i,m) = fft_tmp*fft_tmp';
       end
   end
end

modulation = gtf2gtfcc(modulation,1,31);

cep = [cep; modulation];

del = deltas(cep);

ddel = deltas(deltas(cep,5),5);

cepDpDD = [cep;del;ddel];

end
