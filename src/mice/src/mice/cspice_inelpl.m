function [ nxpts, xpt1, xpt2] = cspice_inelpl( ellips, plane )
   switch nargin
      case 2
         ellips = zzmice_ell(ellips);
         plane   = zzmice_pln(plane);
      otherwise
         error ( ['Usage: [ nxpts, xpt1, xpt2] = ' ...
                         ' cspice_inelpl( ellips, plane )'] )
   end
   try
      [ nxpts, xpt1, xpt2] = mice( 'inelpl_c', ellips, plane );
   catch spiceerr
      rethrow(spiceerr)
   end