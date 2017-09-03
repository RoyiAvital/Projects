function test_suite = im2html_test
initTestSuite

function grayscaleScalar_test
s = im2html(1.0);
assertHtmlTableSize(s, 1, 1);

function grayscaleRectangular_test
s = im2html([0.5 1; 1 0.5; 0.25 0.75]);
assertHtmlTableSize(s, 3, 2);

function binaryImage_test
s = im2html(logical([1 0 1]));
assertEqual(getTableCellContents(s, 1, 2), '0');
assertEqual(getTableCellBackgroundColor(s, 1, 3), '#FFFFFF');

function rgbImage_test
rgb = cat(3,[1 1 1],[0.5 0.5 0.5],[0 0 0]);
s = im2html(rgb);
assertHtmlTableSize(s, 1, 3);
assertTrue(~isempty(regexp(getTableCellContents(s, 1, 1), '1\.0')));
assertEqual(getTableCellBackgroundColor(s, 1, 2), '#FF8000');

function indexedImage_test
X = [1 2 3];
map = [1 0 0; 0 1 0; 0 0 1];
s = im2html(X,map);
assertEqual(getTableCellBackgroundColor(s, 1, 1), '#FF0000');
assertEqual(getTableCellBackgroundColor(s, 1, 3), '#0000FF');

function fixedScaleGrayscale_test
I = [0 2; 3 4];
s = im2html(I,[0 2]);
assertEqual(getTableCellBackgroundColor(s, 2, 1), '#FFFFFF');

function autoScaleGrayscale_test
I = [0 2; 3 4];
s = im2html(I,[]);
assertEqual(getTableCellBackgroundColor(s, 1, 2), '#808080');
assertEqual(getTableCellBackgroundColor(s, 2, 2), '#FFFFFF');

function defaultShowPixelValues_test
s = im2html(1);
assertEqual(getTableCellContents(s, 1, 1), '1');

function showPixelValuesOn_test
s = im2html(1, 'ShowPixelValues', 'on');
assertEqual(getTableCellContents(s, 1, 1), '1');

function showPixelValuesOff_test
s = im2html(1, 'ShowPixelValues', 'off');
assertTrue(isempty(getTableCellContents(s, 1, 1)));

function stylesDefault_test
s = im2html(1);
assertTrue(hasStyleBlock(s));

function styleSheetEmpty_test
s = im2html(1,'StyleSheet','');
assertTrue(hasStyleBlock(s));

function customStyleSheet_test
% Style sheet doesn't need to actually exist in order to verify correct
% functionality here.
s = im2html(1,'StyleSheet','bogus.css');
assertFalse(hasStyleLinkBlock(s));

function outputFile_test
filename = [tempname '.html'];
cleaner = onCleanup(@() delete(filename));
im2html([1 1 1], 'OutputFile', filename);
s = fileread(filename);
assertHtmlTableSize(s, 1, 3);

function maximumDimensionDefault_test
s = warning('off', 'se:im2html:ImageTooBig');
cleaner = onCleanup(@() warning(s));
s = im2html(rand(30,1));
assertHtmlTableSize(s, 20, 1);

s = im2html(rand(1,30));
assertHtmlTableSize(s, 1, 20);

function pixelCellClassGrayscale_test
s = im2html([0.1 0.9]);
assertEqual(getTableCellClass(s, 1, 1), 'dark-pixel-cell');
assertEqual(getTableCellClass(s, 1, 2), 'light-pixel-cell');

function pixelCellClassRgb_test
s = im2html(repmat([0.1 0.9], 1, 3));
assertEqual(getTableCellClass(s, 1, 1), 'dark-pixel-cell');
assertEqual(getTableCellClass(s, 1, 2), 'light-pixel-cell');

function assertHtmlTableSize(s, M, N)
[MM, NN] = getHtmlTableSize(s);
assertEqual(M, MM);
assertEqual(N, NN);

function [M,N] = getHtmlTableSize(s)
num_tr_tags = numel(regexp(s, '<tr'));
num_td_tags = numel(regexp(s, '<td'));
M = num_tr_tags;
N = num_td_tags / M;

function text = getTableCellContents(s, r, c)
row_matches = regexp(s, '<tr.*?>(.*?)</tr>', 'match');
rth_row = row_matches{r};
td_contents_matches = regexp(rth_row, '<td.*?>(.*?)</td>', 'tokens');
text = td_contents_matches{c};
text = text{1};
text = regexprep(text,'^\s*','');
text = regexprep(text,'\s*$','');

function text = getTableCellAttributes(s, r, c)
row_matches = regexp(s, '<tr.*?>(.*?)</tr>', 'match');
rth_row = row_matches{r};
td_attributes_matches = regexp(rth_row, '<td(.*?)>', 'tokens');
text = td_attributes_matches{c};
text = text{1};
text = regexprep(text, '^s\*', '');
text = regexprep(text, '\s*$', '');

function text = getTableCellClass(s, r, c)
td_attributes = getTableCellAttributes(s, r, c);
td_class = regexp(td_attributes, 'class="(.*?)"', 'tokens');
td_class = td_class{1};
if isempty(td_class)
    text = '';
else
    text = td_class{1};
end

function text = getTableCellBackgroundColor(s, r, c)
td_attributes = getTableCellAttributes(s, r, c);
td_bgcolor = regexp(td_attributes, 'background-color\s*:\s*(\#[\dA-Fa-f]{6})', 'tokens');
td_bgcolor = td_bgcolor{1};
if isempty(td_bgcolor)
    text = '';
else
    text = td_bgcolor{1};
end

function tf = hasStyleBlock(s)
tf = ~isempty(regexp(s, '<head>.*?<style>.*?</style>.*?</head>', 'once'));

function tf = hasStyleLinkBlock(s)
tf = ~isempty(regexp(s, '<head>.*?<link>.*?</link>.*?</head>', 'once'));
