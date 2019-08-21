function res = mtimes(a,b)

if a.adjoint
    res = d_adjoint_operator(b, a.x, a.TE, a.fatModulation, a.FT, a.coilmaps, a.Nread, a.Nproj, a.Ncoil, a.prec);    
else
    res = d_forward_operator(b, a.x, a.TE, a.fatModulation, a.FT, a.coilmaps, a.Nread, a.Nproj, a.Ncoil, a.prec);
end
