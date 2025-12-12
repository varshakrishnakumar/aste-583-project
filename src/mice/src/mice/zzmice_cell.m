function [y] = zzmice_cell(x, cell_type)
   if( ~isequal(nargin,2) )
      error( 'MICE(USAGE): _y_ = zzmice_cell( _x_, `cell_type`)' )
   end
   if( isequal( x, [] ) )
      y = cast( zeros(0,1), cell_type);
      return;
   end
   x_dim     = size(x);
   [ N, col] = size(x);
   if( isnumeric(x) && (col == 1) && isequal(size(x_dim), [1,2]) )
     y = cast( x, cell_type);
     return;
   else
      error( ['MICE(BADARG): Improper type of input '        ...
              'argument passed to function. Input expected ' ...
              'as numeric Nx1 array.' ] )
   end