function [y] = zzmice_str( x, empty )
   if( isequal(nargin,1) )
      empty = false;
   elseif( ~isequal(nargin,2) )
      error( 'MICE(USAGE): _`y`_ = zzmice_str( _`x`_, [empty] )' )
   end
   if( iscellstr(x) )
      x = char(x);
   end
   if( ischar(x) )
      y = x;
   else
      error( [ 'MICE(BADARG): Improper type of input ' ...
               'argument passed to function. Value ' ...
               'or values expected as string cell or ' ...
               'character array.' ] )
   end
   if( ~empty && numel(x) == 0  )
      error( 'MICE(BADARG): Attempt to use null string as input.' )
   end