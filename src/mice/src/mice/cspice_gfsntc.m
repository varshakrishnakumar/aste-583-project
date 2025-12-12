function [result] = cspice_gfsntc( target, fixref,  method, abcorr,   ...
                                   obsrvr, dref,    dvec,   crdsys,   ...
                                   coord,  relate,  refval, adjust,   ...
                                   step,   nintvls, cnfine        )
   switch nargin
      case 15
         target  = zzmice_str(target);
         fixref  = zzmice_str(fixref);
         method  = zzmice_str(method);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         dref    = zzmice_str(dref);
         dvec    = zzmice_dp(dvec);
         crdsys  = zzmice_str(crdsys);
         coord   = zzmice_str(coord);
         relate  = zzmice_str(relate);
         refval  = zzmice_dp(refval);
         adjust  = zzmice_dp(adjust);
         step    = zzmice_dp(step);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );
         cnfine  = zzmice_win(cnfine);
      otherwise
         error ( [ 'Usage: [result] = cspice_gfsntc( `target`, `fixref`, ' ...
                                             '`method`, `abcorr`, '        ...
                                             '`obsrvr`, `dref`, dvec[3], ' ...
                                             '`crdsys`, `coord`, '         ...
                                             '`relate`, refval, adjust, '  ...
                                             'step, nintvls, cnfine )' ])
   end
   try
      [result] = mice('gfsntc_c', target, fixref,  method, abcorr, ...
                                  obsrvr, dref,    dvec,   crdsys, ...
                                  coord,  relate,  refval, adjust, ...
                                  step,   nintvls, [zeros(6,1); cnfine] );
   catch spiceerr
      rethrow(spiceerr)
   end