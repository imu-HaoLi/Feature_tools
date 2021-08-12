function [cep, output1] = rasmfcc( time_seq )
    %return matrix auto, each line column is the autocorrelation of a frame
    L = 2; % filter length
    
    win = 320;
    offset = win/2;
    n_frame = floor(length(time_seq)/offset);   
    % turin input into horizonal vector
    time_seq = reshape(time_seq,1,length(time_seq));
    auto = zeros(win,n_frame);
    output1 = zeros(win,n_frame);
          
    %special processing of the first frame
    x = [ zeros(1,offset) time_seq( 1:offset)];    
    temp = xcorr(x,x);
    auto_slice = temp(length(x):end);        
    auto(:,1) = auto_slice';
    
    for i=2:n_frame
        x = time_seq(offset*(i-2)+1:offset*(i-2)+win);    
        temp = xcorr(x,x);
        auto_slice = temp(length(x):end);        
        auto(:,i) = auto_slice';
    end  
    

    for i=1:size(auto,1)
        auto(i,:) = auto(i,:)*1/(win-(i-1));
    end
    
    T_L = sum((-L:L).^2);    
    pad_auto = [zeros(win, L) auto zeros(win, L)];
    
    for i=(L+1):n_frame+L
        sum1 = 0;
        for j=-L:L
            sum1 = sum1 + j*pad_auto(:,i+j);
        end
        output1(:,i-L) = sum1/T_L;
    end
    
    %cep = zeros(93, size(output1,2));
    
    for i=1:size(output1,2)
        
        cep(1:31,i)=melfcc(output1(:,i), 16000, 'maxfreq', 8000, 'numcep', 31, 'nbands', 20, 'fbtype', 'htkmel', 'dcttype', 3, 'wintime', 0.020, 'hoptime', 0.010, 'sumpower',0, 'lifterexp', -22);
  
    end
    
%    cep(32:62,:) = deltas(cep(1:31,:));
 %   
  %  cep(63:93,:) = deltas(deltas(cep(1:31,:),5),5);
   
   
	cep = cep(: , 2:end);
end
