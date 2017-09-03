function s = im2html(varargin)
%im2html Convert image pixels to HTML table showing pixel colors and values
%   SYNTAX
%
%       im2html(I)
%       im2html(I,[])
%       im2html(I,[LOW HIGH])
%       im2html(BW)
%       im2html(RGB)
%       im2html(X,MAP)
%
%       im2html(___, Name, Value)
%
%       s = im2html(___)
%
%   DESCRIPTION
%
%   im2html produces HTML text containing a table that displays the image
%   pixel colors and values. When called with no output argument, im2html
%   displays the HTML in the MATLAB Web Browser. When called with an output
%   argument, im2html returns the HTML text as a string. This form can be
%   useful when publishing a MATLAB file.
%
%   im2html takes the following optional name/value pairs:
%
%       ShowPixelValues    'on' (the default) or 'off'
%
%       StyleSheet         '' (the default) or string containing URL
%                          When '', the output HTML includes a <style>
%                          block in the <head> section that includes styles
%                          for the <table> and <td> classes. When a URL,
%                          the output HTML contains a link to the specified
%                          style sheet.
%
%       OutputFile         string
%                          When specified, im2html writes the HTML text
%                          into the desired file, and it does not
%                          automatically display the HTML file in the
%                          system browser.
%
%       MaximumDimension   positive integer (default 20)
%                          The largest number of rows or columns that
%                          im2html will write out. Rows or columns above
%                          this value will be ignored.
%
%   EXAMPLES
%
%   Display a table of values for a gray-scale image.
%
%       I = imread('pout.tif');
%       im2html(I(125:134, 104:114))
%
%   Display a table of values for a gray-scale image using auto-scale
%   syntax.
%
%       I = magic(10);
%       im2html(I,[])
%
%   Display a table of values for an RGB image.
%
%       rgb = imread('peppers.png');
%       im2html(rgb(88:97,200:209,:))
%
%   Save generated HTML text to a file.
%
%       rgb = imread('peppers.png');
%       im2html(rgb(88:97,200:209,:),'OutputFile','sample.html')
%
%   Display generated HTML text to the Command Window. This can be used to
%   embed an image pixel table into the output of a published MATLAB
%   file.
%
%       [X,map] = imread('trees.tif');
%       disp(im2html(X(1:4,1:4),map))
%
%   Use a custom style sheet.
%
%       I = imread('pout.tif');
%       im2html(I(125:134, 104:114), ...
%           'StyleSheet', ...
%           'http://blogs.mathworks.com/images/steve/2011/im2html_sample_stylesheet.css')    

%   Steve Eddins
%   Copyright 2011 The MathWorks, Inc.

[imshow_args, options] = parseInputs(varargin{:});

%
% Display image into invisible figure using imshow. Create a clean-up
% object so that the figure gets deleted no matter what errors happen
% later. Get an imagemodel object, which is used later to get screen RGB
% pixel values and formatted pixel-value text.
%
fig = figure('Visible', 'off', 'HandleVisibility', 'off');
cleaner = onCleanup(@() close(fig));
ax = axes('Parent', fig);
model = imagemodel(imshow(imshow_args{:}, 'Parent', ax));

%
% Initialize object to accumulate HTML text.
%
html = se.utils.HtmlText;

openTag(html, 'html');

%
% Put style information into the <head>.
% 
openTag(html, 'head');
if ~isempty(options.StyleSheet)
    % Include a link to the specified stylesheet.
    addLine(html, sprintf('<link rel="stylesheet" type="text/css" href="%s" />', ...
        options.StyleSheet));
else
    % Add default styles into the generated HTML.
    addStyles(html);
end
closeTag(html, 'head');

openTag(html, 'body');

openTag(html, 'table', tableAttributes());

M = getImageHeight(model);
N = getImageWidth(model);

if ((M > options.MaximumDimension) || (N > options.MaximumDimension))
    warning('se:im2html:ImageTooBig', ...
        'Using only the first %d rows and columns of the input image.', ...
        options.MaximumDimension);
    M = min(M, options.MaximumDimension);
    N = min(N, options.MaximumDimension);
end

if strcmp(options.ShowPixelValues, 'on')
    [pixel_text, max_text_width, max_text_height] = getPixelValueText(model);
    cell_width_em = max(max_text_width, max_text_height);
else
    cell_width_em = 3;
end

for p = 1:M
    openTag(html, 'tr', 'class="image-table-row"');
    
    for q = 1:N
        comment(html, sprintf('ROW: %d  COL: %d', p, q));
        openTag(html, 'td', ...
            tableCellAttributes(getScreenPixelRGBValue(model,p,q), ...
            cell_width_em));
        if strcmp(options.ShowPixelValues, 'on')
            addLine(html, pixel_text{p,q});
        end
        closeTag(html, 'td');
        addBlankLine(html);
    end
    
    closeTag(html, 'tr');
    addBlankLine(html);
end
closeTag(html, 'table');
closeTag(html, 'body');
closeTag(html, 'html');

ss = toString(html);

processOutput(ss, nargout, options);

if nargout > 0
    s = ss;
end

%==========================================================================
function [pixel_text, max_text_width, max_text_height] = getPixelValueText(model)
%   Takes an imagemodel object and returns a p-by-q cell array of pixel
%   text labels, suitably filtered and processed for direct insertion into
%   an HTML table <td> tag. Also returns the maximum number of characters
%   on any text line, and maximum number of text lines for any pixel.

cellLabelFun = getPixelRegionFormatFcn(model);
M = getImageHeight(model);
N = getImageWidth(model);

pixel_text = cell(M, N);
for p = 1:M
    for q = 1:N
        text_pq = cellLabelFun(p, q);
        pixel_text{p,q} = text_pq{1};
    end
end

[max_text_width, max_text_height] = findBiggestText(pixel_text);

for p = 1:M
    for q = 1:N
        label = pixel_text{p,q};
        label = strrep(label, ' ', '&nbsp;');
        label = strrep(label, '<', '&lt;');
        label = regexprep(label, '\n', '<br \\>');
        pixel_text{p,q} = label;
    end
end
    
%==========================================================================
function [imshow_args, options] = parseInputs(varargin)
%   Returns a cell array of arguments to be used to pass to imshow. Returns
%   an options struct containing name-value pair parsing results.

first_string_idx = find(cellfun(@ischar, varargin, 'UniformOutput', true), 1, 'first');
if isempty(first_string_idx)
    first_string_idx = numel(varargin) + 1;
end

imshow_args = varargin(1:(first_string_idx-1));
name_value_args = varargin(first_string_idx:end);

parser = getNameValueParser;
parse(parser, name_value_args{:});

options = parser.Results;

%==========================================================================
function parser = getNameValueParser()
%   Returns the (cached) name-value pair parser for im2html.

persistent cached_parser

if isempty(cached_parser)
    cached_parser = createNameValueParser();
end

parser = cached_parser;

%==========================================================================
function parser = createNameValueParser()
%    Creates a new name-value pair parser for im2html.

parser = inputParser;
addParamValue(parser, 'ShowPixelValues', 'on');
addParamValue(parser, 'StyleSheet', '');
addParamValue(parser, 'MaximumDimension', 20);
addParamValue(parser, 'OutputFile', '');

%==========================================================================
function [max_width, max_height] = findBiggestText(c)
%   Given a cell array of pixel-value label text, return the maximum number
%   of characters on any line of the labels for all the pixels. Return the
%   maximum number of lines of labels for all the pixels.

max_width = 0;
max_height = 0;
[M,N] = size(c);
for p = 1:M
    for q = 1:N
        lines = textscan(c{p,q}, '%s', 'Delimiter', '\n');
        lines = lines{1};
        max_height = max(max_height, numel(lines));
        max_width = max(max_width, max(cellfun('prodofsize', lines)));
    end
end

%==========================================================================
function addStyles(html)
%   Write CSS style info into the HTML text object.

addBlankLine(html);

openTag(html, 'style');
addLine(html, '.image-table {');
addLine(html, '  border: thin solid black;');
addLine(html, '  table-layout: fixed;');
addLine(html, '}');
addBlankLine(html);

addLine(html, '.light-pixel-cell {');
addLine(html, '  font-size: small;');
addLine(html, '  font-family: Consolas, ''Courier New'', monospace;');
addLine(html, '  border-style: solid;');
addLine(html, '  border-width: thin;');
addLine(html, '  border-color: black;');
addLine(html, '  padding: 1;');
addLine(html, '  height: 4.5em;');
addLine(html, '  width: 4.5em;');
addLine(html, '  text-align: center;');
addLine(html, '  color: black;');
addLine(html, '}');
addBlankLine(html);

addLine(html, '.dark-pixel-cell {');
addLine(html, '  font-size: small;');
addLine(html, '  font-family: Consolas, ''Courier New'', monospace;');
addLine(html, '  border-style: solid;');
addLine(html, '  border-width: thin;');
addLine(html, '  border-color: black;');
addLine(html, '  padding: 1;');
addLine(html, '  height: 4.5em;');
addLine(html, '  width: 4.5em;');
addLine(html, '  text-align: center;');
addLine(html, '  color: white;');
addLine(html, '}');
addBlankLine(html);

closeTag(html, 'style');
addBlankLine(html);

%==========================================================================
function s = tableCellAttributes(rgb, cell_width_em)
%   Given a single triple of RGB values and the desired cell width in em
%   units, return a <td> attributes string. The attributes string sets the
%   background color of the table cell, sets the cell width and height, and
%   sets the <td> class to be either light-pixel-cell or dark-pixel-cell
%   depending on the gray-scale equivalent value of the RGB color.

gray = rgb2gray(im2double(rgb));
if gray(1) > 0.5
    cell_class = 'light-pixel-cell';
else
    cell_class = 'dark-pixel-cell';
end
cell_width_string = sprintf('%.1fem', cell_width_em);
s = sprintf('class="%s" style="background-color:%s; width:%s; height:%s"', ...
    cell_class, se.utils.rgbToHtmlHex(rgb), ...
    cell_width_string, cell_width_string);

%==========================================================================
function s = tableAttributes()
%   Return a <table> attributes string. The attributes string sets the
%   <table> class to image-table and sets the cellspacing to "1". The
%   attribute cellspacing is used instead of CSS styles because of uneven
%   browser support.

s = 'class="image-table"';

%==========================================================================
function processOutput(s, num_outputs, options)
%   If the user has asked for an output file, then save the HTML string to
%   a file. Otherwise, if the user did not request an output argument,
%   display the HTML file using the system browser.

if ~isempty(options.OutputFile)
    filename = options.OutputFile;
    fid = fopen(filename, 'w');
    fprintf(fid, '%s', s);
    fclose(fid);
    return;
end

if num_outputs == 0
    web(['text://', s]);
end
