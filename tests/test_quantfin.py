from quantfin import npv, irr
assert abs(npv(0.1,[-100,60,60])-4.13)<0.05
assert abs(irr([-100,60,60])-0.1306)<1e-3
