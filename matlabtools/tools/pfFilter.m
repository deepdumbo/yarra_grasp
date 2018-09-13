function filter = pfFilter(filterLen, slopeLen)

if nargin<2
    lSlopeLen = filterLen / 8;
    rSlopeLen = lSlopeLen;
else
    [lSlopeLen, rSlopeLen] = slopeLen;
end

filter = ones(1,filterLen);

filter(1:lSlopeLen)         = sin( (1:lSlopeLen) * pi / (2*lSlopeLen) );
filter(end-rSlopeLen+1:end) = sin( (rSlopeLen:-1:1) * pi / (2*rSlopeLen) );


