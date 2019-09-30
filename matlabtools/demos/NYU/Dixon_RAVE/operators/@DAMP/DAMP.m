function res = DAMP(x, TE, fatModulation, FT, coilmaps, imsize, Nread, Nproj, Ncoil, prec)

res.adjoint         = 0;
res.x               = x;
res.TE              = TE;
res.fatModulation   = fatModulation;
res.imsize          = imsize;
res.FT              = FT;
res.coilmaps        = coilmaps;
res.Nread           = Nread;
res.Nproj           = Nproj;
res.Ncoil           = Ncoil;
res.prec            = prec;

res                 = class(res,'DAMP');

