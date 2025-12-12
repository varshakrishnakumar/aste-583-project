function [ids] = cspice_pckfrm( pckfnm, room, ids_i )
   switch nargin
      case 2
         pckfnm = zzmice_str(pckfnm);
         room   = zzmice_int(room);
      case 3
         pckfnm = zzmice_str(pckfnm);
         room   = zzmice_int(room);
         ids_i  = zzmice_int(ids_i);
         is_valid =  (  (ndims(ids_i) == 2) && (size(ids_i, 2) == 1) )     ...
                        ||                                                 ...
                        isequal( size(ids_i), [0,0] );
         if (~is_valid )
            error( 'MICE(BADARG): Argument ''ids_i'' must have size Nx1.' )
         end
      otherwise
         error( 'Usage: [ids] = cspice_pckfrm( _`pckfnm`_, room, [ids_i] )' )
   end
if ( nargin == 2 )
   try
      ids = mice('pckfrm_c', pckfnm, room );
   catch spiceerr
      rethrow(spiceerr)
   end
else
   try
      ids = unique( [ [ids_i]; mice('pckfrm_c', pckfnm, room ) ]  );
   catch spiceerr
      rethrow(spiceerr)
   end
end