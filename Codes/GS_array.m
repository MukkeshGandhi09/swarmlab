Mx = 16;
My = 16;
c = physconst('LightSpeed');
freq = 10^7;

lambda = c/freq;
delx = lambda/2;
dely = lambda/2;

X = zeros(1,Mx*My);
Y = zeros(1,Mx*My);
Z = zeros(1,Mx*My);

for p = 1:Mx
    for q = 1:My
        i = (q-1)*Mx+p;
        X(i) = (p-1)*delx;
        Y(i) = (q-1)*dely;
        Z(i) = 0;        
    end
end

