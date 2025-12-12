function [srfstr, isname] = cspice_srfcss(code, bodstr)
   switch nargin
      case 2
         code   = zzmice_int(code);
         bodstr = zzmice_str(bodstr);
      otherwise
         error ( 'Usage: [`srfstr`, isname] = cspice_srfcss(code, `bodstr`)' )
   end
   try
      [srfcss] = mice('srfcss_s', code, bodstr);
      srfstr   = char( srfcss.name );
      isname   = reshape( [srfcss.found], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end