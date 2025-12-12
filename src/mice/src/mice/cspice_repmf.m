function [out] = cspice_repmf( in, marker, value, sigdig, format )
   switch nargin
      case 5
         in     = zzmice_str(in, true);
         marker = zzmice_str(marker, true);
         value  = zzmice_dp(value);
         sigdig = zzmice_int(sigdig);
         format = zzmice_str(format);
      otherwise
         error ( [ 'Usage: [`out`] = '                                      ...
                   'cspice_repmf( `in`, `marker`, value, sigdig, `format` ' ...
                   ')' ] )
   end
   try
      [out] = mice('repmf_c', in, marker, value, sigdig, format);
   catch spiceerr
      rethrow(spiceerr)
   end