function [utcstr] = cspice_et2utc(et, format, prec )
   switch nargin
      case 3
         et     = zzmice_dp(et);
         format = zzmice_str(format);
         prec   = zzmice_int(prec);
      otherwise
         error ( 'Usage: [_`utcstr`_] = cspice_et2utc(_et_, `format`, prec)' )
   end
   try
      [utcstr] = mice('et2utc_c',et,format,prec);
   catch spiceerr
      rethrow(spiceerr)
   end