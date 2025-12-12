function [y] = zzmice_dp( x, nanok )
   if( isequal(nargin,1) )
      nanok = false;
   elseif( ~isequal(nargin,2) )
      error( 'MICE(USAGE): [_y_] = zzmice_dp( _x_, [nanok] )' )
   end
   if( ~isnumeric(x) )
      error( [ 'MICE(BADARG): Improper type of input argument passed to '  ...
               'function. Value or values expected as double precision '   ...
               'or integer.' ] )
   end
   if( ~nanok && ~all( isfinite( x(:) ) ) )
      error( [ 'MICE(NOTFINITE): Improper type of input argument passed '  ...
               'to function. Value or values expected as finite double '   ...
               'precision or integer.' ] )
   end
   y = double(x);