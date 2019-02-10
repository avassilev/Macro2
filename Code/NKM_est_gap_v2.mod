var y trend_y cycle_y Pi betta;

varexo v1 v2 v3 v4 ;

parameters c1 c2 c3;
c1 = 0.7;
c2 = 0.6;
c3 = 1;



model(linear);

y = trend_y + cycle_y;
trend_y = trend_y(-1) + betta(-1) + v1;
betta = betta(-1) + v2;
cycle_y = c1*cycle_y(-1) + v3;

Pi = (1-c2)*Pi(+1) + c2*Pi(-1) + c3*cycle_y + v4;


end;

initval;
y = 8.001577;
trend_y = y;
betta = 0;
cycle_y  = 0;
Pi = 0.01;
end;


steady(nocheck);

estimated_params;
stderr v1     , inv_gamma_pdf,  0.005  , inf;
stderr v2     , inv_gamma_pdf,  0.005  , inf;
stderr v3     , inv_gamma_pdf,  0.005  , inf;
stderr v4     , inv_gamma_pdf,  0.005  , inf;

c1        , normal_pdf   ,  0.7    , 0.03;
c2        , beta_pdf   ,  0.6    , 0.2;
c3        , normal_pdf   ,  1   , 1;
end;

varobs y Pi; 


estimation(datafile=estdata,mh_replic=2000, mode_compute=4, mh_nblocks=2, mh_jscale=0.3, filtered_vars, smoother, diffuse_filter) trend_y cycle_y Pi; 
