function h_kl=get_polarization_loss(antennaorient,droneorient,dronepos,antennapos,f0)
dlen=1;
rotmp=eul2rotm(droneorient,'ZYX');
rotmx=eul2rotm(antennaorient,'ZYX');
c=physconst('LightSpeed');
p=dronepos-antennapos;
x1=mtimes(rotmx,p);
dkl=norm(p);
thetax_kl=acos(x1(3)/dkl);
psix_kl=acos(x1(2)/dkl);
x1=mtimes(rotmp,-p);
Rxp=(rotmx.')*rotmp;
thetap_kl=acos(x1(3)/dkl);
psip_kl=acos(x1(2)/dkl);
constmat=[0,0;0,1;1,0];
x=p(1);
y=p(2);
z=p(3);
magxy=norm([x,y]);
magxz=norm([x,z]);
firstmat=1/dkl*[ -x*z/magxy, -y*z/magxy,    magxy; ...
                 -x*y/magxz,      magxz, -z*y/magxz ];
T=firstmat*Rxp*constmat;

    function f=F(theta,f)
        if mod(theta,pi)==0
            f=0;
        else
            f=(cos((pi*dlen/(c/f))*cos(theta))-cos(pi*dlen/(c/f)))/sin(theta);
        end
    end

E_ltheta=1;
E_lpsi=1;
E_ktheta=1;
E_kpsi=1;
El=[(E_ltheta*F(thetax_kl,f0));(E_lpsi*F(psix_kl,f0))];
Ek=[(E_ktheta*F(thetap_kl,f0));(E_kpsi*F(psip_kl,f0))];

h_kl=El'*T*Ek;


end