function s = rgbToHtmlHex(rgb)
%rgbToHtmlHex Convert RGB triple to HTML hex string

%   Steve Eddins
%   Copyright 2011 The MathWorks, Inc.

rgb = rgb(:)';
rgb = im2uint8(rgb);
s = dec2hex(rgb, 2)';
s = ['#' s(:)'];
