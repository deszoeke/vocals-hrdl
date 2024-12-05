function a = anom(x, varargin)
a = x - nanmean(x, varargin{:});
end
