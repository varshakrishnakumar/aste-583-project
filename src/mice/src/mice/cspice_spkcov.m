function [cover] = cspice_spkcov( spkfnm, idcode, room, cover_i )
   switch nargin
      case 3
         spkfnm = zzmice_str(spkfnm);
         idcode = zzmice_int(idcode);
         room   = zzmice_int(room, [1, int32(inf)/2] );
      case 4
         spkfnm  = zzmice_str(spkfnm);
         idcode  = zzmice_int(idcode);
         room    = zzmice_int(room, [1, int32(inf)/2] );
         cover_i = zzmice_win(cover_i);
      otherwise
         error ( [ 'Usage: [cover] = cspice_spkcov( _`spkfnm`_, '          ...
                                     'idcode, room, [cover_i] )' ] )
   end
if ( nargin == 3 )
   try
      [cover] = mice('spkcov_c', spkfnm, idcode, room );
   catch spiceerr
      rethrow(spiceerr)
   end
else
   try
      cover = mice('spkcov_c', spkfnm, idcode, room );
   catch spiceerr
      rethrow(spiceerr)
   end
   cover = cspice_wnunid( cover, cover_i );
end