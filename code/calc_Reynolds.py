#!/usr/bin/env python3

#at standard temperature and pressure, converted to SI units (i.e. m-kg-sec)
D = 1e-6 #micrometer
rho = 0.9982 * 1e-3 * 1e6 #g/cm^3  1g=10^-3kg, 1/cm^3*(100cm/1m)^3
v = 1e-6 #micrometer / sec

#"kinematic viscosity" - not the right units
#eta = 0.01 #cm^2 / sec

#"dynamics viscosity" - same units as Taylor
eta = 1.0016e-3 #mPa *sec = mN/m^2 *sec


#unit check:  m*m/sec*kg/(m^3)*1/(N/m^2*sec)
#blech! 
#unit check:  kg/m * m^2 / (N*sec^2) 
#unit check:  m*kg / (N*sec^2) = m*kg / (kg*m s^2 /s^2)
#kg cancels!  m cancels! s^2 cancels!
#unitless Reynold's number


#1 micrometer = 10^-6 meters
#1 centimeter = 10^-2 meters
#1 micrometer = 10^-4 centimeters

R = D*v*rho/eta
print(R)
