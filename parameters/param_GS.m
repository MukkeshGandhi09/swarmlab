p_GS.is_GS_required=true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for the ground station
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if p_GS.is_GS_required
    p_GS.Mx = 4;
    p_GS.My = 4;
    p_GS.N0 = 2e-20;         % Noise power
    p_GS.pu = 10;            % SNR
    p_GS.bw = 1e7;          % Bandwidth
    p_GS.freq = 2.4e9;       % operating frequency
    c = physconst('LightSpeed');
    p_GS.tau_dl=0;
    p_GS.tau_ctrl=0;
    p_GS.Bc=3e6;
    
    lambda = c/p_GS.freq;
    p_GS.delx = lambda/2;
    p_GS.dely = lambda/2;
    
    p_GS.X = zeros(1,p_GS.Mx*p_GS.My);
    p_GS.Y = zeros(1,p_GS.Mx*p_GS.My);
    p_GS.Z = zeros(1,p_GS.Mx*p_GS.My);
    
    for p = 1:p_GS.Mx
        for q = 1:p_GS.My
            i = (q-1)*p_GS.Mx+p;
            p_GS.X(i) = (p-1)*p_GS.delx;
            p_GS.Y(i) = (q-1)*p_GS.dely;
            p_GS.Z(i) = 0;
        end
    end
    p_GS.Az = [10 20 30 10 40 0 60 23 88 78 10 90 120 160 60 79];
    p_GS.El = [23 88 78 10 90 12 60 60 79 10 20 30 10 40 0 60];
end
clear c p q i lambda