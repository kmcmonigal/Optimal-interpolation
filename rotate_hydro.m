% This will rotate u and v components to the angle of the hydrographic line

function [across,along,dist,angle] = rotate_hydro(lat,lon,u,v)

[dist,angle] = sw_dist(lat,lon,'km');

along = u.*cos(angle)+v.*sin(angle);
across = -u.*sin(angle)+v.*cos(angle);

