% Author: Kaivalya Bakshi; Date: 18 May 2017; v1 Date 18 May 2017: v mean for 2 particles plotted

function comm_rate = cr(beta, x_j, x_i)
% Cucker Smale communication rate psi(a,b) = 1/(1 + |a - b|^2)^beta

comm_rate = 1/(1 + norm(x_j - x_i)^2)^beta;

end