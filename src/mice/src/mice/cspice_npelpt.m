function [ pnear, dist ] = cspice_npelpt( point, ellips )
   switch nargin
      case 2
         point  = zzmice_dp(point);
         ellips = zzmice_ell(ellips);
      otherwise
         error ( 'Usage: [pnear, dist] = cspice_npelpt( point, ellips )' )
   end
   try
      [npelpt] = mice( 'npelpt_s', point, ellips );
      pnear    = reshape( [npelpt.pos], 3, [] );
      dist     = reshape( [npelpt.alt], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end