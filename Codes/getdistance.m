function d = getdistance(antennapos,dronepos)
    d=sqrt(sum((antennapos-dronepos ).^2));
end