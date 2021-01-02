function gk = get_channelVector(GSantennapos,dronepos,pathloss,freq)
c = 3e8;
lambda = c/freq;
M = length(GSantennapos);
K = length(dronepos);
for k=1:K
    for i=1:M
        x = GSantennapos(i)-dronepos(k);
        d_kl(k,i) = sqrt(sum(x.^2));
        gk(k,i) = sqrt(pathloss(k,i))*exp(-1i*2*pi*d_kl(i)/lambda);
    end
end
end