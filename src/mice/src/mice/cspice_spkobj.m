function [ids] = cspice_spkobj( spkfnm, room, ids_i )
   switch nargin
      case 2
         spkfnm = zzmice_str(spkfnm);
         room   = zzmice_int(room);
      case 3
         spkfnm = zzmice_str(spkfnm);
         room   = zzmice_int(room);
         ids_i  = zzmice_int(ids_i);
         is_valid =  (  (ndims(ids_i) == 2) && (size(ids_i, 2) == 1) )     ...
                        ||                                                 ...
                        isequal( size(ids_i), [0,0] );
         if (~is_valid )
            error( 'MICE(BADARG): Argument ''ids_i'' must have size Nx1.' )
         end
      otherwise
         error( 'Usage: [ids] = cspice_spkobj( _`spkfnm`_, room, [ids_i])' )
   end
if ( nargin == 2 )
   try
      ids = mice('spkobj_c', spkfnm, room );
   catch spiceerr
      rethrow(spiceerr)
   end
else
   try
      ids = unique( [ [ids_i]; mice('spkobj_c', spkfnm, room ) ]  );
   catch spiceerr
      rethrow(spiceerr)
   end
end