function [dvnorm] = cspice_dvnorm(state)
   switch nargin
      case 1
         state = zzmice_dp(state);
      otherwise
         error ( 'Usage: [_dvnorm_] = cspice_dvnorm(_state(6)_)' )
   end
   try
      [dvnorm] = mice('dvnorm_c', state);
   catch spiceerr
      rethrow(spiceerr)
   end