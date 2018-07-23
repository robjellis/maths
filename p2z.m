function [zval] = p2z(pval)

% convert a P-value to the corresponding Z-value (for ease of display;
% also, makes more sense than a -log10 transform)

zval = -1 * norminv(pval);

