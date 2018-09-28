%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to perform Coil Compression using singular value decomposition.
%
% Marnix Maas, RadboudUMC, 27-9-2018
% Adapted from Li Feng's Demo1_Unstreaking.m, provided with the paper:
% Benkert T, Tian Y, Huang C, DiBella EVR, Chandarana H, Feng L.
% Magn Reson Med. 2017 Nov 28. doi: 10.1002/mrm.27030.
% 
% Input parameters:
% - kdata:        k-space data. Shoud have same order of dimensions as is
%                 read in by mapVBVD, ie [nx,nc,ntviews,nz,(ne)]
% - ncc:          desired number of 'combined' coil elements
% - doErr:        Flag to select whether the error introduced by coil
%                 compression should be calculated. This is done by
%                 calculating the Frobenius norm of the difference between
%                 compressed and original data matrices as a function of 
%                 the number of combined coil elements, and is output in
%                 the parameter err. Default is 0.
% 
% Output parameters:
% - kdata_out     Compressed k-space data
% - S             Singular values found during compression
% - err           Compression error as a function of number of combined
%                 coil elements
%
% USE: [kdata_out, S, err] = coilCompress(kdata, ncc, doErr)

function [kdata_out, S, err] = coilCompress(kdata, ncc, doErr)
% Hoe dit werkt:
% - Beschrijf elk datapunt als een vector in coil-ruimte (elke coil
% vertegenwoordigt een as in deze ruimte)
% - Schrijf de data daartoe om naar een Nd (datapunten) x Nc (coils)
% matrix, die de datapunten dus als rijvectoren bevat.
% - Bereken de SVD van deze ruimte.
% - De hieruit volgende matrix V bevat de z.g. 'right singular vectors',
% waarvan de kolommen de eenheidsvectoren bevatten van een nieuwe basis
% voor de datavectoren. Deze basisvectoren zijn lineaire combinaties van de
% oorspronkelijke eenheidsvectoren, die zo zijn gekozen dat de variatie in
% de data langs de betreffende as steeds zo groot mogelijk is, afgezien van
% die langs eerdere assen. (Dus variatie langs as 1 is grootst, langs as 2
% kleiner, etc, maar altijd maximum bereikbare van de overgebleven
% dimensies). In dit geval kunnen deze basisvectoren worden gezien als
% 'gecombineerde spoelelementen'.
% Het product van een data- (=rij-)vector met een kolom van V geeft de
% coordinaat van die vector langs de betreffende V-as. Product van de
% rijvector met matrix V geeft dus dezelfde vector uitgedrukt in de basis
% V, en product D*V geeft alle datavectoren uitgedrukt in de basis V.
% - Waar V de nieuwe basisvectoren bevat op volgorde van data-variatie,
% bevat het SVD-resultaat S de grootte van die variaties. Latere elementen
% van S bevatten lagere waarden, die dus tot uitdrukking brengen dat de
% data langs de bijbehorende V-vector weinig variatie (=informatie) bevat.
% Dit betekent weer dat weglaten van die componenten weinig gevolgen heeft
% voor de daadwerkelijke inhoud van de data
% - Bekijk dus de 'singular values' (uit S), en besluit welke wel en niet
% zullen worden meegenomen
% - Vermenigvuldig D met de kolommen van V die de meeste informatie
% bevatten: data is gecomprimeerd tot een aantal 'gecombineerde
% spoelelementen', dat dus lager is dan het totaal.

% Parse inputs
if nargin<3
    doErr = 0;
end

% Permute dimensions:
% from [nx,nc,ntviews,nz,...]
% to   [nx,ntviews,nz,...,nc]
ndk     = length(size(kdata));              % Number of dimensions
order   = [1 3:ndk 2];                      % Permutation order
kdata   = permute(kdata,order);

% Get data dimensions
sk      = size(kdata);
nc      = sk(end);                          % Other elements contain nx, 
                                            % ntviews etc, but this ensures 
                                            % independence from exact data 
                                            % dimensionality

% Reshape data, do SVD
D=reshape(kdata,prod(sk(1:end-1)),nc);
[U,S,V]=svd(D,'econ');

% Optional, calculate Frobenius norm with original data to assess effect
% of compression on data quality as function of number of 'combined' coil
% elements used
if doErr
    for j=1:nc
        s1 = diag(S);
        s1(j:end) = 0;
        S1 = diag(s1);
        err(j) = norm(U*S1*V' - D, 'fro');
    end
end

% Project data onto new axes & keep only those we want
D = D*V(:,1:ncc);

% Free up some memory
clear U V kdata
% Reshape back
kdata_out = reshape(D, [sk(1:end-1), ncc]);
% Permute back
order = [1 ndk 2:ndk-1];
kdata_out = permute(kdata_out, order);
S = diag(S);                                % Provide singular values in vector form

    
    
    
    