function y = zzmice_logical (x)
   if( ~isequal(nargin,1) )
      error( 'MICE(USAGE): [_y_] = zzmice_logical( _x_ )' )
   end
   if (islogical(x) && ~any(isnan(x)) )
      y = x;
   elseif( isnumeric(x) && all(isfinite(x))  )
      y =  [ double(x) | zeros(size(x)) ];
   else
      error( ['MICE(BADARG): Improper type of input ' ...
              'argument passed to function. Value ' ...
              'or values expected as finite numeric or logical.'] )
   end