function rgb = randomRGB(numChannels, maxVal, minVal, RNG)
%RANDOMRGB generates random colors as hex strings
%   s = RANDOMRGB() will produce a random three-channel color like '07EDBA'.
%
%   RANDOMRGB(NUMCHANNELS) to specify the number of channels, e.g., 4 for RGBA.
%
%   RANDOMRGB(NUMCHANNELS, MAXVAL, MINVAL) allows one to specify the max
%   (default 255) and min (default 0) values allowed. if MAXVAL < MINVAL, they
%   are swapped.
%
%   RANDOMRGB(NUMCHANNELS, MAXVAL, MINVAL, RNG) specifies the random number
%   generator stream. If omitted, the default MATLAB-wide global stream is used.
%
%   If any of these four arguments is omitted or is empty, the default values
%   will be used.

if ~exist('numChannels', 'var') || isempty(numChannels), numChannels = 3; end
if ~exist('maxVal', 'var') || isempty(maxVal), maxVal = 255; end
if ~exist('minVal', 'var') || isempty(minVal), minVal = 0; end
if ~exist('RNG', 'var') || isempty(RNG), RNG = RandStream.getGlobalStream(); end

color = RNG.randi(sort([minVal maxVal]), 1, numChannels);
rgbCell = arrayfun(@(i) sprintf('%02s', dec2hex(i)), color, 'un', 0);
rgb = strcat(rgbCell{:});
