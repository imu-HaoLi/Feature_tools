function [cep, auto] = acmfcc( time_seq )
    %return matrix auto, each line column is the autocorrelation of a frame
    win = 320;
    offset = win/2;
    n_frame = floor(length(time_seq)/offset);   
    % turin input into horizonal vector
    time_seq = reshape(time_seq,1,length(time_seq));
    auto = zeros(win,n_frame);
    
    preemph = 0.97;
    if preemph ~= 0
        time_seq = filter([1 -preemph], 1, time_seq);
    end
    
    %special processing of the first frame
    x = [ zeros(1,offset) time_seq( 1:offset)]; 
    x = x.*transpose(hamming(win)); % hamming windowing
    temp = xcorr(x,x);
    auto_slice = temp(length(x):end);        
    auto(:,1) = auto_slice';
    
    for i=2:n_frame
        x = time_seq(offset*(i-2)+1:offset*(i-2)+win);    
        x = x.*transpose(hamming(win)); % hamming windowing
        temp = xcorr(x,x);
        auto_slice = temp(length(x):end);        
        auto(:,i) = auto_slice';
    end  
    
    auto = auto./win;
    
    %discard the lower lags (the lags < 2ms) remove first 32 lags
    low_lags = 32;
    
    % create double dynamic range    
    temp_double_ham = hamming((win - low_lags)/2);
    double_hamming = [xcorr(temp_double_ham,temp_double_ham) ; 0];
    
    high_lag_auto = [auto(low_lags+1:win,:).*repmat(double_hamming,1,size(auto,2)) ; zeros(low_lags,size(auto,2))];
    
    cep = zeros(39, size(high_lag_auto,2));

    for i=1:size(high_lag_auto,2)
        cep(1:31,i)=mag_melfcc(high_lag_auto(:,i), 16000, 'maxfreq', 8000, 'numcep', 31, 'nbands', 20, 'fbtype', 'htkmel', 'dcttype', 3, 'wintime', 0.020, 'hoptime', 0.010, 'sumpower',0, 'lifterexp', -22);
    end
    
    %cep(32:62,:) = deltas(cep(1:31,:));
    %cep(63:93,:) = deltas(deltas(cep(1:31,:),5),5);
    
    cep = cep(: , 2:end);
end
