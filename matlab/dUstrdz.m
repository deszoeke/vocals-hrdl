function dUdz = dUstrdz(u, v, dz)
% uvec . duvec/dz / ||uvec||
dUdz = ( midz(u) .* diff(u, 1, 2) + midz(v) .* diff(v, 1, 2) ) ./ ( sqrt(midz(u).^2 + midz(v).^2) * dz );
end
