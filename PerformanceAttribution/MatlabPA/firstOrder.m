function [g] = firstOrder(fxph, fxmh, h)

g = (fxph - fxmh) / (2*h);

end