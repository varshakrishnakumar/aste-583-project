function [idset] = cspice_bltfrm( frmcls, room )
   switch nargin
      case 2
         frmcls = zzmice_int(frmcls);
         room = zzmice_int(room, [1, int32(inf)] );
      otherwise
         error ( 'Usage: [idset] = cspice_bltfrm( frmcls, room )' )
   end
   try
      [idset] = mice('bltfrm_c', frmcls, room);
   catch spiceerr
      rethrow(spiceerr)
   end