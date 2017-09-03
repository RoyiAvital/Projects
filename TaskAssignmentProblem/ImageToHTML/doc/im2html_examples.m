%% im2html Examples

%%
% Display a table of values for a gray-scale image.

I = imread('pout.tif');
disp(im2html(I(125:134, 104:114)))

%%
% Display a table of values for a gray-scale image using auto-scale syntax.

I = magic(10);
disp(im2html(I,[]))

%%
% Display a table of values for an RGB image.

rgb = imread('peppers.png');
disp(im2html(rgb(88:97,200:209,:)))

%%
% Display a table of values for an indexed image.

[X,map] = imread('trees.tif');
disp(im2html(X(1:4,1:4),map))

%%
% Suppress the display of pixel values.

I = imread('pout.tif');
disp(im2html(I(125:134, 104:114),'ShowPixelValues','off'))

%%
% Steve Eddins
% Copyright 2011 The MathWorks, Inc.