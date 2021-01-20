% Obtain the projected points with plane norm.
% hui.chen@kaust.edu.sa
function output = project_nd(plane_norm, point)
    f = plane_norm;
    p = point;
    t0 = -(p'*f)/(norm(f)^2);
    output = p + f*t0;
end
