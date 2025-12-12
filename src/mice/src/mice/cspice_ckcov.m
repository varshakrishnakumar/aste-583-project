function [cover] = cspice_ckcov( ckfnm, idcode, needav, level,             ...
                                 tol,   timsys, room,   cover_i )
   switch nargin
      case 7
         ckfnm  = zzmice_str(ckfnm);
         idcode = zzmice_int(idcode);
         needav = zzmice_int(needav);
         level  = zzmice_str(level);
         tol    = zzmice_dp(tol);
         timsys = zzmice_str(timsys);
         room   = zzmice_int(room, [1, int32(inf)/2] );
      case 8
         ckfnm  = zzmice_str(ckfnm);
         idcode = zzmice_int(idcode);
         needav = zzmice_int(needav);
         level  = zzmice_str(level);
         tol    = zzmice_dp(tol);
         timsys = zzmice_str(timsys);
         room   = zzmice_int(room, [1, int32(inf)/2] );
         cover_i= zzmice_win(cover_i);
      otherwise
         error ( ['Usage: [cover] = cspice_ckcov( _`ckfnm`_, idcode, '     ...
                          'needav, level, tol, timsys, room, [cover_i] )'] )
   end
if ( nargin == 7 )
   try
      cover = mice('ckcov_c', ckfnm, idcode, needav,                       ...
                   level,     tol,   timsys, room    );
   catch spiceerr
      rethrow(spiceerr)
   end
else
   try
      cover = mice('ckcov_c', ckfnm, idcode, needav,                       ...
                   level,     tol,   timsys, room );
   catch spiceerr
      rethrow(spiceerr)
   end
   cover = cspice_wnunid( cover, cover_i );
end