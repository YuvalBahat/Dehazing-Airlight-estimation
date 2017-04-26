function modulus = one2Nmod(x,y)
% modulus = one2Nmod(x,y)
% Same as mod, but the result is in the range [1,N] instead of [0,N-1].
% This means it is the same as mod except for those values for which mod
% returns 0 - Those will now return N.

modulus = mod(x-1,y)+1;