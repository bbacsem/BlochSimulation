function angle_list = phase_shift_angle(i,angle)

for k = 1:i
phi = (k - 1) * k / 2 * angle / 360;
phi = phi - round(phi);
angle_list(k) =  2.0 * pi * phi;
end