function P = points_on_line(p1, p2, t_vals)
    
    if nargin < 3
        t_vals = 2/3;
    end

    P = zeros(length(t_vals), length(p1));
    for i = 1:length(t_vals)
        t = t_vals(i);
        P(i,:) = t*p1 + (1 - t)*p2;
    end
end