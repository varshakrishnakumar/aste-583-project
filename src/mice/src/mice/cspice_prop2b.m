function [pvprop] = cspice_prop2b( gm, pvinit, dt )
   switch nargin
      case 3
         gm = zzmice_dp(gm);
         pvinit = zzmice_dp(pvinit);
         dt = zzmice_dp(dt);
      otherwise
         error ( 'Usage: [pvprop(6)] = cspice_prop2b( gm, pvinit(6), dt )' )
   end
   try
      [pvprop] = mice('prop2b_c', gm, pvinit, dt);
   catch spiceerr
      rethrow(spiceerr)
   end