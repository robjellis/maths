function [sl_vals] = sl(vals)

% signed log values

sl_vals = sign(vals) .* log(abs(vals)+1);