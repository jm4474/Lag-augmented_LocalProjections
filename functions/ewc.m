function Omega = ewc(Y, bw)

    % Equal-Weighted Cosine HAR estimator
    % See Lazarus, Lewis, Stock & Watson (JBES 2018), p. 545
    
    % Here calculated equivalently using Fast Fourier Transform for speed
    
    % Inputs:
    % Y         T x n   data matrix
    % bw        1 x 1   integer bandwidth
    
    % Outputs:
    % Omega     n x n   HAR var-cov matrix 
    
    
    [T,n] = size(Y);
    rffts = zeros(bw,n);
    
    for j=1:n
        the_Yt = reshape([zeros(T,1) Y(:,j)]',2*T,1); % Insert 0s before every obs.
        the_rfft = real(fft(the_Yt,4*T)); % Compute real part of FFT, padding with T trailing 0s
        rffts(:,j) = the_rfft(2:bw+1,:);
    end
    
    Omega = (2/bw)*(rffts'*rffts);
    
end