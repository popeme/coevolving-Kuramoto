function Ybold_reg = convert_to_BOLD(ths,lts,step,N)

load('/home/spornslab/Documents/Maria/Kura_comm/Data/BOLDHRF.mat')

hrf = BOLDHRF(1:20000);                     % define length of hrf (20 sec)
lhrf = length(hrf);

% Convert to BOLD (convolution method), lowpass filter, and downsample
Ybold_ds = zeros(N,lts);
for n=1:N
    Yboldtemp = conv(sin(ths(n,:)),hrf,'valid'); % 'valid'
    Yboldtemp = lowpass(Yboldtemp,0.25,1000,'Steepness',0.95); % lowpass, 0.25 Hz
    for t=1:lts
        Ybold_ds(n,t) = mean(Yboldtemp((t-1)*step+1:t*step));
    end
end
ths = ths(:,lhrf+1:end);  % discard first 'lhrf' to align with 'valid'

% Regress global mean from downsampled Ybold
Yboldglob = mean(Ybold_ds);
Ybold_reg = zeros(lts,N);
for n=1:N
    [~,~,Ybold_reg(:,n)] = regress(Ybold_ds(n,:)',[ones(1,lts)', Yboldglob']);
end
end