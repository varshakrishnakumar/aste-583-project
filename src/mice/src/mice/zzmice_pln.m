function [y] = zzmice_pln(x)
   if( ~isequal(nargin,1) )
      error( 'MICE(USAGE): y = zzmice_pln( x)' )
   end
   if( isstruct(x) )
      if( ~isequal( size( x ), [1,1] )  )
         error( ['MICE(BADARG): Improper type of input '    ...
                'argument passed to zzmice_pln. Incorrect ' ...
                'shape, not [1,1].'] )
      end
      fields = { 'normal'; 'constant' };
      names  = fieldnames( x );
   else
      error( ['MICE(BADARG): Improper type of input ' ...
              'argument passed to zzmice_pln. Not a structure.'] )
   end
   if( isequal(fields,names) )
      if( ~isequal( size( x.normal ), [3,1] ) )
         error( ['MICE(BADFIELD): Improper dimensions for ' ...
                    'normal field, not [3,1].'] )
      end
      if( ~isequal( size( x.constant ), [1,1] ) )
         error( ['MICE(BADFIELD): Improper dimensions for ' ...
                    'constant field, not [1,1].'] )
      end
      y.normal   = zzmice_dp(x.normal);
      y.constant = zzmice_dp(x.constant);
   else
      error( ['MICE(BADARG): Improper type of input '    ...
              'argument passed to zzmice_pln. Improper ' ...
              'form for plane structure.'] )
   end