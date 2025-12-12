function [ettarg, elapsd] = cspice_ltime(etobs, obs, dir, targ)
   switch nargin
      case 4
         etobs = zzmice_dp(etobs);
         obs   = zzmice_int(obs);
         targ  = zzmice_int(targ);
         dir   = zzmice_str(dir);
      otherwise
         error ( ['Usage: [_ettarg_, _elapsd_] = ' ...
                  'cspice_ltime( _etobs_, obs, `dir`, targ)'] )
   end
   try
      [ettarg, elapsd] = mice('ltime_c',etobs, obs, dir, targ);
   catch spiceerr
      rethrow(spiceerr)
   end