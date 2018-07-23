function C = permn_rje(V,N)
nv = length(V) ;
C = zeros(nv^N,N) ; % declaration
for ii=1:N,
    cc = 1 ;
    for jj=1:(nv^(ii-1)),
        for kk=1:nv,
            for mm=1:(nv^(N-ii)),
                C(cc,ii) = V(kk) ;
                cc = cc + 1 ;
            end
        end
    end
end