function fatModulation = calculateFatModulation(par)

par.rel_amp(1) = [];
par.f_wf(1) = [];

%% Modelling of TEs in k-space

fatModulation = zeros(par.Nread,1,par.Neco); % Nread (x Nproj) x Neco

for ie = 1:par.Neco
    %% Calculate timemap
    if par.timemap.usetimemap
        timemap = par.TE(ie) + (- par.Nread/2 : par.Nread/2 - 1) * par.samptime;
        
        % Mirror acquisition times for every even echo due to bipolar readout
        if mod(ie+1,2)
            timemap = timemap(end:-1:1);
        end
        
    else
        timemap = par.TE(ie);
    end
    
    %% Calculate demodulation map
    timemap = timemap.';
    for is = 1:length(par.rel_amp)
        fatModulation(:,1,ie) = fatModulation(:,1,ie) +  par.rel_amp(is) .* exp(1i*2*pi*par.f_wf(is)*timemap);
    end
    
end