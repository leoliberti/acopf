% Y = Ymatrix(r,x,b,tau,theta)
%
%    Construct the complex 2x2 Y matrix for AC Ohm's law I = YV
%      given MatPower's branch parameters r,x,b,tau,theta
function Y = Ymatrix(r,x,b,tau,theta)
    Y = zeros(2);
    ys = 1/(r+i*x);
    Y(1,1) = (ys + i*b/2) / tau^2;
    Y(1,2) = -ys*tau*exp(i*theta);
    Y(2,1) = -ys*tau*exp(-i*theta);
    Y(2,2) = ys + i*b/2;
end
