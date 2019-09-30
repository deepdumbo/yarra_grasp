function d_kdata = d_forward_operator(dx, x, TE, fatModulation, FT, coilmaps, Nread, Nproj, Ncoil, prec)

d_kdata     = complex(zeros(Nread, Nproj, Ncoil, length(TE), prec));

water       = x(:,:,1);
fat         = x(:,:,2);
fieldmap    = x(:,:,3);

d_water     = dx(:,:,1);
d_fat       = dx(:,:,2);
d_fieldmap  = dx(:,:,3);

for ic = 1:Ncoil
    for ie = 1:length(TE)
        
        % Signal model with chemical shift in k-space
        coilmaps_expfieldmap = coilmaps(:,:,ic).*exp(1i*2*pi*fieldmap*TE(ie));
        
        d_kdata_tmp1 = 1i*2*pi*TE(ie) * water .* coilmaps_expfieldmap .* d_fieldmap + ...
            coilmaps_expfieldmap .* d_water;
        
        d_kdata_tmp2 = 1i*2*pi*TE(ie) * fat .* coilmaps_expfieldmap .* d_fieldmap + ...
            coilmaps_expfieldmap .* d_fat;
        
        d_kdata(:,:,ic,ie) = (FT{ie} * d_kdata_tmp1) +  bsxfun(@times, fatModulation(:,1,ie), FT{ie} * d_kdata_tmp2);
        
    end
end
