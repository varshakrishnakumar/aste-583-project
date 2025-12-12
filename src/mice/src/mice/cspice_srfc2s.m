function [srfstr, isname] = cspice_srfc2s(code, bodyid)
   switch nargin
      case 2
         code   = zzmice_int(code);
         bodyid = zzmice_int(bodyid);
      otherwise
         error ( 'Usage: [`srfstr`, found] = cspice_srfc2s(code, bodyid)' )
   end
   try
      [srfc2s] = mice('srfc2s_s', code, bodyid);
      srfstr   = char( srfc2s.name );
      isname   = reshape( [srfc2s.found], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end