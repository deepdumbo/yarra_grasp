function kdata = forward_operator(x, coilmaps, TE, fatModulation, FT, Nread, Nproj, Ncoil, prec)
% Forward operator for non-Cartesian fat-water separation for multiple coils

kdata  = complex(zeros(Nread, Nproj, Ncoil, length(TE), prec));

water       = x(:,:,1);
fat         = x(:,:,2);
fieldmap    = x(:,:,3);

for ic = 1:Ncoil
    for ie = 1:length(TE)
                
        % Signal model with chemical shift in k-space
        kdata(:,:,ic,ie) = ...
            FT{ie} * (coilmaps(:,:,ic).*exp(1i*2*pi*fieldmap*TE(ie)).*water) + ...
            bsxfun(@times, fatModulation(:,1,ie), FT{ie} * (coilmaps(:,:,ic).*exp(1i*2*pi*fieldmap*TE(ie)).*fat));
        
    end
end