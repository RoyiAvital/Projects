%Deque Double-ended queue

%   Steve Eddins
%   Copyright 2011 The MathWorks, Inc.

classdef Deque < handle
    
    properties (Dependent = true)
        Items
    end
    
    properties (Access = private)
       array = cell(1, 10); 
       first = 1;
       last = 0;
       GrowthFactor = 2;
       PackFraction = 0.25;
    end
    
    properties (SetAccess = private)
        NumItems = 0
    end
    
    methods
        function value = get.Items(self)
            value = self.array(self.itemIndices());
        end
    end
    
    methods (Access = private)
    
        function p = arrayIsFull(self)
            p = (self.NumItems == numel(self.array));
        end
        
        function growArray(self)
            K = numel(self.array);
            self.array{max(ceil(self.GrowthFactor * K), 2)} = [];
            if self.first > self.last
                num_items_to_move = K - self.first + 1;
                self.array(end - num_items_to_move + 1:end) = ...
                    self.array(self.first:K);
            end
        end
        
        function packArray(self)
            self.array = self.array(self.itemIndices);
            self.first = 1;
            self.last = self.NumItems;
        end
        
        function incrementFirst(self)
            self.first = moduloIncrementOneBased(self.first, numel(self.array));
        end
        
        function incrementLast(self)
            self.last = moduloIncrementOneBased(self.last, numel(self.array));
        end
        
        function decrementFirst(self)
            self.first = moduloDecrementOneBased(self.first, numel(self.array));
        end
        
        function decrementLast(self)
            self.last = moduloDecrementOneBased(self.last, numel(self.array));
        end
        
        function indices = itemIndices(self)
            indices = mod((1:self.NumItems) + self.first - 2, numel(self.array)) + 1;
        end
    end
      
    methods
        
        function p = isEmpty(self)
           p = self.NumItems == 0; 
        end
        
        function pushBack(self, item)
            if self.arrayIsFull()
                self.growArray();
            end
            self.incrementLast();
            self.array{self.last} = item;            
            self.NumItems = self.NumItems + 1;
        end
        
        function pushFront(self, item)
           if self.arrayIsFull()
               self.growArray();
           end
           self.decrementFirst();
           self.array{self.first} = item;
           self.NumItems = self.NumItems + 1;
        end
        
        function val = popBack(self)
            if self.isEmpty()
               error('Deque:PopBackOnEmpty', ....
                   'popBack method called on empty deque.');
            end
            val = self.array{self.last};
            self.decrementLast();
            self.NumItems = self.NumItems - 1;
            
            if self.NumItems < self.PackFraction * numel(self.array);
                self.packArray();
            end
        end
        
        function val = popFront(self)
            if self.isEmpty()
               error('Deque:PopFrontOnEmpty', ...
                   'popFront method called on empty deque.');
            end
            val = self.array{self.first};
            self.incrementFirst();
            self.NumItems = self.NumItems - 1;
            
            if (self.NumItems / numel(self.array)) <= self.PackFraction
                self.packArray();
            end
        end
        
        function val = back(self)
            if self.isEmpty()
                error('Deque:BackOnEmpty', ...
                    'back method called on empty deque.');
            end
            val = self.array{self.last};
        end
        
        function val = front(self)
            if self.isEmpty()
                error('Deque:FrontOnEmpty', ...
                    'front method called on empty deque.');
            end
            val = self.array{self.first};
        end
    end
    
end

function np = moduloIncrementOneBased(n, N)
np = mod(n, N) + 1;
end

function np = moduloDecrementOneBased(n, N)
np = mod(n-2, N) + 1;
end