function [ code, found] = cspice_srfscc( srfstr, bodyid )
   switch nargin
      case 2
         srfstr = zzmice_str(srfstr);
         bodyid = zzmice_int(bodyid);
      otherwise
         error ( 'Usage: [ code, found] = cspice_srfscc( `srfstr`, bodyid )' )
   end
   try
      [srfscc] = mice('srfscc_s', srfstr, bodyid);
      code   = reshape( [srfscc.code],  1, [] );
      found    = reshape( [srfscc.found], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end