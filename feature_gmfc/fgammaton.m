function [f, tmp]=fgammaton(sig, gamma, sampFreq, numChannel , window_size)
% Filter input signal with a gammatone filterbank and generate GF features
% sig: input signal
% gamma: gammtone filter impulse responses
% sampFreq: sampling frequency
% numChannel: number of channels
% Written by Yang Shao, and adapted by Xiaojia Zhao in Oct'11
% adapted by masood in 2016

if nargin < 5
       window_size = 10;
end

[ignore,stepBound]=size(sig);

d = 2^ceil(log2(length(sig)));

tmp =ifft((fft((ones(numChannel,1)*sig)',d).*fft(gamma(:,1:stepBound)',d)));

if window_size == 10
    f=(abs(resample(abs(tmp(1:stepBound,:)),100,sampFreq)')).^(1/3);
elseif window_size == 50
    in10ms = resample(abs(tmp(1:stepBound,:)),100,sampFreq);
    in10ms = [ ones(2,1)*in10ms(1 , :); in10ms ; ones(2,1)*in10ms(end , :) ];
    
    tmp = [];
    for i=1:(size(in10ms,1)-4)
      tmp = [tmp ; sum(in10ms(i:i+4 , :))];
    end
    
    dLambda = 0.4;

    d_c_0 = 0.01;

    tmp2 = zeros(size(tmp));
    for l = 1:size(tmp,1) %time frame index
       for m = 1:64 % channel
           ad_mu = [];
            if m == 1
              ad_M = tmp(l , m);
            else
              % Equation (3) in the IS2010 paper
              ad_M  = dLambda * ad_M + (1 - dLambda) * tmp(l , m);
              ad_P2 = max(tmp(l , m) - ad_M,  d_c_0 * ad_M);
              d_w =  ad_P2 / max(tmp(l , m), realmin('double')); 
              tmp2(l , m) = d_w * tmp(l , m);
            end
       end
    end
    
    f = abs(tmp2').^(1/3);
else
    disp('not supported');
    assert(1==2);
end

%%%%%%%%%%%%%%%%%

%r = tmp(1:stepBound,:).';
%
%winLength = 320;      % default window length in sample points which is 20 ms for 16 KHz sampling frequency
%
%[numChan,sigLength] = size(r);     % number of channels and input signal length
%
%winShift = winLength/2;            % frame shift (default is half frame)
%increment = winLength/winShift;    % special treatment for first increment-1 frames
%M = floor(sigLength/winShift);     % number of time frames
%
%% calculate energy for each frame in each channel
%f = zeros(numChan,M);
%for m = 1:M
%    for i = 1:numChan
%        if m < increment        % shorter frame lengths for beginning frames
%            f(i,m) = r(i,1:m*winShift)*r(i,1:m*winShift)';
%        else
%            startpoint = (m-increment)*winShift;
%            f(i,m) = r(i,startpoint+1:startpoint+winLength)*r(i,startpoint+1:startpoint+winLength)';
%        end
%    end
%end
%
%f = f.^(1/3);
