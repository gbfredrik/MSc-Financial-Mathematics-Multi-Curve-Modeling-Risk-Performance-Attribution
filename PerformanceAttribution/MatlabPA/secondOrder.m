function [H] = secondOrder(fx, fxph, fxmh, h)

H = (fxph - 2*fx + fxmh) / h^2;

end