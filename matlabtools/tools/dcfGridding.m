function dcf = dcfGridding(spokes, nx)
% Calculate the density compensation matrix for gridding reconstruction

dcfRow=ones(nx,1);
for i=1:nx
   dcfRow(i)=abs(nx/2-(i-0.5));
end
% dcfRow=pi/spokes*dcfRow;                % Not sure what's the reason for this normalization
dcfRow = dcfRow/max(dcfRow);            % ...just normalize to 1 instead
dcf=repmat(dcfRow,1,spokes);

return