function [rectan] = cspice_radrec( range, ra, dec )
   switch nargin
      case 3
         range = zzmice_dp(range);
         ra    = zzmice_dp(ra);
         dec   = zzmice_dp(dec);
      otherwise
         error ( ['Usage: [_rectan(3)_] = cspice_radrec( _range_, ',  ...
                                                    '_ra_, _dec_ )'] )
   end
   try
      [rectan] = mice('radrec_c', range, ra, dec);
   catch spiceerr
      rethrow(spiceerr)
   end