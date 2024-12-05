function S = SProd2D(w, Urel, dUstrdz, dt, dz)
% S = SProd2D(w, Urel, dUstrdz, dt, dz)
% Shear Production assuming nondivergent
% 2D z-x streamwise flow.
%
% Urel relative to lidar is for Gallilean transform to get uprime
% dUstrdz is the streamwise shear component

% quantities defined on mid-z levels
dwdz = diff(w, 1, 2) / dz;
uprime = dt * midz(Urel) .* anom( cumsum(dwdz) );
wprime = anom( midz(w) );

S = dUstrdz .* nanmean( uprime .* wprime );
end

% internal functions - don't shadow external ones
% function xmid = midz(x)
% xmid = (x(:, 1:end-1) + x(:, 2:end)) * 0.5
% end

function m = nanmean(x, varargin)
m = mean(x, varargin{:}, "omitnan");
end

% function a = anom(x, varargin)
% a = x - nanmean(x, varargin{:})
% end

% use this for argument dUstrdz
% function dUdz = dUstrdz(u, v, dz)
% % uvec . duvec/dz / ||uvec||
% dUdz = ( midz(u) .* diff(u, 1, 2) + midz(v) .* diff(v, 1, 2) ) / ( sqrt(midz(u).^2 .+ midz(v).^2) * dz );
%end
