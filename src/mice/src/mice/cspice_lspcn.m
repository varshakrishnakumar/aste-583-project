function [lspcn] = cspice_lspcn( body, et, abcorr )
   switch nargin
      case 3
         body   = zzmice_str(body);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);
      otherwise
         error( 'Usage: [_lspcn_] = cspice_lspcn( `body`, _et_, `abcorr`)' )
   end
   try
      [lspcn] = mice('lspcn_c', body, et, abcorr );
   catch spiceerr
      rethrow(spiceerr)
   end