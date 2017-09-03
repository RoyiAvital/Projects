%HtmlText HTML text container

%   Steve Eddins
%   Copyright 2011 The MathWorks, Inc.

classdef HtmlText < handle
    
    properties (Hidden)
        Queue
        OpenTags
    end
    
    methods
        function self = HtmlText
            self.Queue = se.utils.Deque;
            self.OpenTags = se.utils.Deque;
        end
        
        function addLine(self, line)
            self.Queue.pushBack(sprintf('%s%s', ...
                se.utils.spaces(self.indentationLevel), ...
                line));
        end
        
        function openTag(self,tag_type,tag_attributes)
            tag_str = tag_type;
            if nargin >= 3
                tag_str = sprintf('%s %s',tag_str,tag_attributes);
            end
            line = sprintf('%s<%s>',...
                se.utils.spaces(self.indentationLevel()),...
                tag_str);
            pushBack(self.Queue, line);
            
            pushBack(self.OpenTags, tag_type);
        end
        
        function closeTag(self,tag_type)
            last_tag = popBack(self.OpenTags);
            if ~strcmp(last_tag, tag_type)
                warning('se:HtmlText:MismatchedTags',...
                    'Tag type "%s" does not match currently open tag type "%s".', ...
                    tag_type, last_tag);
            end
            pushBack(self.Queue, ...
                sprintf('%s</%s>', se.utils.spaces(self.indentationLevel()), tag_type));
        end
        
        function comment(self,comment_text)
            addLine(self, sprintf('<!-- %s -->', comment_text));
        end
        
        function addBlankLine(self)
            self.Queue.pushBack('');
        end
        
        function level = indentationLevel(self)
            level = 2*self.OpenTags.NumItems;
        end
        
        function s = toString(self)
            lines = self.Queue.Items;
            s = [sprintf('%s\n', lines{1:end-1}), lines{end}];
        end
    end
    
end