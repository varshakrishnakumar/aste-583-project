function [cmat, av, clkout, found] = cspice_ckgpav(inst, sclkdp, tol, ref)
   switch nargin
      case 4
         inst   = zzmice_int(inst);
         sclkdp = zzmice_dp(sclkdp);
         tol    = zzmice_dp(tol);
         ref    = zzmice_str(ref);
      otherwise
         error ( [ 'Usage: [_cmat(3,3)_, _av(3)_, _clkout_, _found_] = ' ...
                   'cspice_ckgpav(inst, _sclkdp_, tol, `ref`)' ] )
   end
   try
      [cmat, av, clkout, found] = mice('ckgpav_c', inst, sclkdp, tol, ref);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end