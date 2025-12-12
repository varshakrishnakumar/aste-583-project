function [center, smajor, sminor] = cspice_el2cgv( ellips )
   switch nargin
      case 1
         ellips = zzmice_ell(ellips);
      otherwise
         error ( ['Usage: [center(3), smajor(3), sminor(3)] = ' ...
                  'cspice_el2cgv( ellips )'] )
   end
   try
      [center, smajor, sminor] = mice('el2cgv_c', ellips );
   catch spiceerr
      rethrow(spiceerr)
   end