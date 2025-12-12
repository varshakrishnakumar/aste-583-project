function [state, lt] = cspice_spkcpo(target, et,     outref, ...
                                     refloc, abcorr, obspos, ...
                                     obsctr, obsref )
   switch nargin
      case 8
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         outref = zzmice_str(outref);
         refloc = zzmice_str(refloc);
         abcorr = zzmice_str(abcorr);
         obspos = zzmice_dp(obspos);
         obsctr = zzmice_str(obsctr);
         obsref = zzmice_str(obsref);
      otherwise
         error ( ['Usage: [ state(6), lt] = cspice_spkcpo( ' ...
                          '`target`, et, `outref`, '         ...
                          '`refloc`, `abcorr`, obspos(3), '  ...
                          '`obsctr`, `obsref` )'] )
   end
   try
      [starg] = mice('spkcpo_s',target, et, outref, refloc, ...
                                abcorr, obspos, obsctr, obsref);
      state   = reshape( [starg.state], 6, [] );
      lt      = reshape( [starg.lt   ], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end