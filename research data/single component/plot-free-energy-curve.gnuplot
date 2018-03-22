reset
set xrange [-0.3:1.7]
unset key

tc = 0.4
theta = tc*0.833323529411598
rho1 = 0.290682751988810
rho2 = 1.851178911767298

rho(x) = rho1 + x*(rho2-rho1)
		
pc = 3*tc/8
a = 27*tc*tc/(64*pc)
b = tc/(8*pc)

psi(x) = rho(x)*theta*log(rho(x)/(1-b*rho(x))) - a*rho(x)*rho(x) - theta*rho(x)

plot psi(x)
