function angle = getangle(dronepos,antennapos,droneorient)

rotm=eul2rotm(droneorient,'ZYX');
x=antennapos-dronepos;
x1=mtimes(rotm,x);
angle_arg=[x1(1)/norm(x1),x1(2)/norm(x1),x1(3)/norm(x1)];
angle=[acos(angle_arg(1)),acos(angle_arg(2)),acos(angle_arg(3))];

end
