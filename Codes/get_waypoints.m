function waypoints=get_waypoints(pos,Drones_x,Drones_y,x1,x2,y1,y2,GSD,FL,PS)
H=GSD*FL/PS;
rpy=2664;
rpx=1496;
dx=rpx*GSD;
dy=rpy*GSD;
Tx=abs((x2-x1)/Drones_x);
Ty=abs((y2-y1)/Drones_y);

if dy<Ty
    Y=[dy/2 Ty-(dy/2)];
    j=1;
    X=dx/2;
    waypoints=zeros(length(X)*length(Y),3);
    m=1;
    n=2;
    while X<Tx
        waypoints(j , :)=[X,Y(m),H];
        waypoints(j+1,:)=[X,Y(n),H];
        m=m+n;
        n=m-n;
        m=m-n;
        X=X+dx;
        j=j+2;
    end
else
    waypoints(: , 1)=dx/2:dx:Tx;
    waypoints(: , 2)=dy/2;
    waypoints(: , 3)=H;
end
if isempty(waypoints)
    waypoints=[Tx/2,Ty/2,H];
end
waypoints=pos+waypoints;
end