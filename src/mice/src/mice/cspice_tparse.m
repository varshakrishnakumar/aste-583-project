function [sp2000, errmsg] = cspice_tparse( string )
   switch nargin
      case 1
         string = zzmice_str(string);
      otherwise
         error ( 'Usage: [sp2000, `errmsg`] = cspice_tparse( `string` )' )
   end
   try
      [sp2000, errmsg] = mice('tparse_c', string);
   catch spiceerr
      rethrow(spiceerr)
   end