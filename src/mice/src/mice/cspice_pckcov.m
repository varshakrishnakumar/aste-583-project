function [cover] = cspice_pckcov( pckfnm, idcode, room, cover_i )
   switch nargin
      case 3
         pckfnm  = zzmice_str(pckfnm);
         idcode  = zzmice_int(idcode);
         room    = zzmice_int(room, [1, int32(inf)/2] );
      case 4
         pckfnm  = zzmice_str(pckfnm);
         idcode  = zzmice_int(idcode);
         room    = zzmice_int(room, [1, int32(inf)/2] );
         cover_i = zzmice_win(cover_i);
      otherwise
         error ( [ 'Usage: [cover] = cspice_pckcov( _`pckfnm`_, '          ...
                                     'idcode, room, [cover_i] )' ] )
   end
   if ( nargin == 3 )
      try
         [cover] = mice('pckcov_c', pckfnm, idcode, room );
      catch spiceerr
         rethrow(spiceerr)
      end
   else
      try
         cover = mice('pckcov_c', pckfnm, idcode, room );
      catch spiceerr
         rethrow(spiceerr)
      end
      cover = cspice_wnunid( cover, cover_i );
   end