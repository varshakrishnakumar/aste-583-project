function [state, lt] = cspice_spkcvo(target, et,     outref, ...
                                     refloc, abcorr, obssta, ...
                                     obsepc, obsctr, obsref )
   switch nargin
      case 9
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         outref = zzmice_str(outref);
         refloc = zzmice_str(refloc);
         abcorr = zzmice_str(abcorr);
         obssta = zzmice_dp(obssta);
         obsepc = zzmice_dp(obsepc);
         obsctr = zzmice_str(obsctr);
         obsref = zzmice_str(obsref);
      otherwise
         error ( ['Usage: [ state(6), lt] = cspice_spkcvo( ' ...
                          '`target`, et, `outref`, '         ...
                          '`refloc`, `abcorr`, obssta(6), '  ...
                          'obsepc, `obsctr`, `obsref` )'] )
   end
   try
      [starg] = mice('spkcvo_s',target, et,     outref, ...
                                refloc, abcorr, obssta, ...
                                obsepc, obsctr, obsref );
      state   = reshape( [starg.state], 6, [] );
      lt      = reshape( [starg.lt   ], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end