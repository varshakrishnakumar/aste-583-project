function [bodids] = cspice_dskobj( dskfnm, room, bodids_i )
   switch nargin
      case 2
         dskfnm = zzmice_str(dskfnm);
         room   = zzmice_int(room, [1, int32(inf)/2] );
      case 3
         dskfnm   = zzmice_str(dskfnm);
         room     = zzmice_int(room, [1, int32(inf)/2] );
         bodids_i = zzmice_cell( bodids_i, 'int32' );
      otherwise
         error( ['Usage: [bodids] = '                                      ...
                        'cspice_dskobj( _`dskfnm`_, room, [bodids_i])'] )
   end
if ( nargin == 2 )
   try
      bodids = mice('dskobj_c', dskfnm, room );
   catch spiceerr
      rethrow(spiceerr)
   end
else
   try
      bodids = unique( [ [bodids_i]; mice('dskobj_c', dskfnm, room ) ] );
   catch spiceerr
      rethrow(spiceerr)
   end
end