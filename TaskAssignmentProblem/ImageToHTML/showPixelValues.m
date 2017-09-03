function showPixelValues(varargin)
%   showPixelValues Displays pixels and pixel values like impixelregion
%   showPixelValues displays the pixels of an image, together with a pixel
%   grid and pixel values, like impixelregion does.  Unlike impixelregion,
%   however, showPixelValues produces an ordinary figure that be captured
%   when an M-file script is published.
%
%   Syntaxes:
%       showPixelValues(I)
%       showPixelValues(BW)
%       showPixelValues(RGB)
%       showPixelValues(X,map)
%       showPixelValues(...,center)
%
%   Variable names I, BW, RGB, and X,map represent grayscale, binary, RGB, 
%   and indexed images.
%
%   center is a two-element vector, [X Y], indicating where to center the
%   displayed pixel region.
%
%   Example
%   =======
%       rgb = imread('peppers.png');
%       showPixelValues(rgb, [355 80])

%   Steve Eddins
%   Copyright 2007-2014 The MathWorks, Inc.

error(nargchk(1,3,nargin))

A = varargin{1};
map = [];
center = [];

if (nargin > 1) 
    if numel(varargin{2}) == 3
        map = varargin{2};
    else
        center = varargin{2};
    end
end

if (nargin > 2)
    center = varargin{3};
end

% Display the image.
if isempty(map)
    hIm = imshow(A);
else
    hIm = imshow(A,map);
end

hImageFigure = ancestor(hIm, 'figure');

hPixelRegionFig = impixelregion(hIm);
set(hPixelRegionFig, 'Position', get(0, 'DefaultFigurePosition'));

if ~isempty(center)
    % Adjust the center of the view.  To do this we have to find the
    % imscrollpanel and get its API.
    hPixelRegionIm = findobj(hPixelRegionFig, 'type', 'image');
    if verLessThan('matlab','8.2')
        hScrollable = ancestor(hPixelRegionIm, 'uipanel');
        hImscrollPanel = get(hScrollable, 'Parent');
    else
        hImscrollPanel = ancestor(hPixelRegionIm, 'uipanel');
    end
    spApi = iptgetapi(hImscrollPanel);
    currentMag = spApi.getMagnification();
    spApi.setMagnificationAndCenter(currentMag, center(1), center(2));
end

% Execute the "print to figure" operation.
menu = findobj(hPixelRegionFig, 'Label', '&Print to Figure');
fcn = get(menu, 'Callback');
fcn();

delete(hPixelRegionFig);
delete(hImageFigure);

