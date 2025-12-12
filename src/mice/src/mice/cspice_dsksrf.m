function [srfids] = cspice_dsksrf( dskfnm, bodyid, room, srfids_i )
   switch nargin
      case 3
         dskfnm   = zzmice_str(dskfnm);
         bodyid   = zzmice_int(bodyid);
         room     = zzmice_int(room, [1, int32(inf)/2] );
      case 4
         dskfnm   = zzmice_str(dskfnm);
         bodyid   = zzmice_int(bodyid);
         room     = zzmice_int(room, [1, int32(inf)/2] );
         srfids_i = zzmice_cell(srfids_i, 'int32');
      otherwise
         error ( [ 'Usage: [srfids] = cspice_dsksrf( _`dskfnm`_, '         ...
                                     'bodyid, room, [srfids_i])' ] )
   end
if ( nargin == 3 )
   try
      [srfids] = mice('dsksrf_c', dskfnm, bodyid, room );
   catch spiceerr
      rethrow(spiceerr)
   end
else
   try
      [srfids] = unique( [ [srfids_i];                                     ...
                         mice('dsksrf_c', dskfnm, bodyid, room) ]);
   catch spiceerr
      rethrow(spiceerr)
   end
end