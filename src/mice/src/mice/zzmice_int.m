function [y] = zzmice_int(x, range )
   is_valid = (isnumeric(x) || islogical(x)) && all( isfinite( x(:) ) );
   if( is_valid )
      switch nargin
         case 0
            error( 'MICE(USAGE): _y_ = zzmice_int( _x_, [range])' )
         case 1
            y = int32(x);
         case 2
            y = int32(x);
            if( ~isequal( size(range), [1,2] ) )
               error( ['MICE(BADVAL): The range input requires ' ...
                      'dimension 1x2.' ] )
            end
            if( range(1) >= range(2) )
               error( [ 'MICE(BADVAL): The range input requires ' ...
                       'range(1) < range(2).' ] )
            end
            if ( any( y(:) < range(1) ) || any( y(:) > range(2) ) )
               txt = sprintf( ['MICE(BADVAL): Integer input value not ' ...
                               'within required range [ %ld, %ld ]' ],  ...
                               range(1), range(2) );
               error(txt)
            end
         otherwise
            error( 'MICE(USAGE): [_y_] = zzmice_int( _x_, [range])' )
      end
   else
      error( ['MICE(BADARG): Improper type of input ' ...
              'argument passed to function. Value ' ...
              'or values expected as a finite numeric or logical.'] )
   end