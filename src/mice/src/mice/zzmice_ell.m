function [y] = zzmice_ell(x)
   if( ~isequal(nargin,1) )
      error( 'MICE(USAGE): y = zzmice_ell(x)' )
   end
   if( isstruct(x) )
      if( ~isequal( ndims( x ), 2) )
         error( ['MICE(BADARG): Improper type of input '    ...
                'argument passed to zzmice_ell. Incorrect ' ...
                'dimensions, not [1,N].'] )
      end
      if( ~isequal( size( x, 1 ), 1)  )
         error( ['MICE(BADARG): Improper type of input '    ...
                'argument passed to zzmice_ell. Incorrect ' ...
                'shape, not [1,N].'] )
      end
      fields = { 'center'; 'semiMajor'; 'semiMinor' };
      names  = fieldnames( x );
   else
      error( ['MICE(BADARG): Improper type of input ' ...
              'argument passed to zzmice_ell. Not a structure.'] )
   end
   if( isequal(fields,names) )
      for i=1:numel(x)
         if( ~isequal( size( x(i).center ), [3,1] ) )
            error( ['MICE(BADFIELD): Improper dimensions for ' ...
                    'center field, not [3,1].'] )
         end
         if( ~isequal( size( x(i).semiMajor ), [3,1] ) )
            error( ['MICE(BADFIELD): Improper dimensions for ' ...
                    'semiMajor field, not [3,1].'] )
         end
         if( ~isequal( size( x(i).semiMinor ), [3,1] ) )
            error( ['MICE(BADFIELD): Improper dimensions for ' ...
                    'semiMinor field, not [3,1].'] )
         end
         x(i).center    = double(x(i).center);
         x(i).semimajor = double(x(i).semiMajor);
         x(i).semiminor = double(x(i).semiMinor);
      end
      y = x;
   else
      error( ['MICE(BADARG): Improper type of input '    ...
              'argument passed to zzmice_ell. Improper ' ...
              'form for ellipse structure.'] )
   end