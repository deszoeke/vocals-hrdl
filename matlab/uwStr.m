function uw = uwStr(w, Urel, dt, dz)
% uw = uwStr(w, Urel, dt, dz)
% streamwise stress assuming nondivergent
% 2D z-x streamwise flow.
%
% Urel relative to lidar is for Gallilean transform to get uprime

% quantities defined on mid-z levels
dwdz = diff(w, 1, 2) / dz;
uprime = dt * midz(Urel) .* anom( cumsum(dwdz) );
wprime = anom( midz(w) );

uw = nanmean( uprime .* wprime );
end

% internal functions - don't shadow external ones
% function xmid = midz(x)
% xmid = (x(:, 1:end-1) + x(:, 2:end)) * 0.5
% end

function m = nanmean(x, varargin)
m = mean(x, varargin{:}, "omitnan");
end

% function a = anom(x, varargin)
% a = x - nanmean(x, varargin{:});
% end
