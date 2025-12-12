function [cmat, clkout, found] = cspice_ckgp(inst, sclkdp, tol, ref)
   switch nargin
      case 4
         inst   = zzmice_int(inst);
         sclkdp = zzmice_dp(sclkdp);
         tol    = zzmice_dp(tol);
         ref    = zzmice_str(ref);
      otherwise
         error ( [ 'Usage: [_cmat(3,3)_, _clkout_, _found_] = ' ...
                   'cspice_ckgp(inst, _sclkdp_, tol, `ref`)' ] )
   end
   try
      [cmat, clkout, found] = mice('ckgp_c', inst, sclkdp, tol, ref);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end