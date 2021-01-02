function gain = get_HW_Dipole_AntennaGain(ang)
gain = cos(pi*cos(ang)/2)/sin(ang);
end