function dx = d_adjoint_operator(d_kdata, x, TE, fatModulation, FT, coilmaps, Nread, Nproj, Ncoil, prec)

water       = x(:,:,1);
fat         = x(:,:,2);
fieldmap    = x(:,:,3);

d_water     = complex(zeros(size(water), prec));
d_fat       = complex(zeros(size(fat), prec));
d_fieldmap  = complex(zeros(size(fieldmap), prec));

dx          = complex(zeros(size(x), prec));

for ic = 1:Ncoil
    for ie = 1:length(TE)
        
        % Signal model with chemical shift in k-space
        coilmaps_expfieldmap = coilmaps(:,:,ic).*exp(1i*2*pi*fieldmap*TE(ie));
        
        temp1 = FT{ie}' * d_kdata(:,:,ic,ie);
        temp2 = FT{ie}' * bsxfun(@times, d_kdata(:,:,ic,ie), conj(fatModulation(:,1,ie)));
        
        d_water     = d_water   + conj(coilmaps_expfieldmap) .* temp1;
        d_fat       = d_fat     + conj(coilmaps_expfieldmap) .* temp2;
        d_fieldmap  = d_fieldmap ...
            + conj(1i*2*pi*TE(ie) * coilmaps_expfieldmap .* water) .* temp1 ...
            + conj(1i*2*pi*TE(ie) * coilmaps_expfieldmap .* fat) .* temp2;
        
    end
end

dx(:,:,1) = d_water;
dx(:,:,2) = d_fat;
dx(:,:,3) = d_fieldmap;

end
