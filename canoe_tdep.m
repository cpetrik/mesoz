Eap = 37.4*1e3; %37.4 kJ/mol = 37.4*1e3 J/mol
R = 8.314; %J/mol/K gas constant

T1 = 298.15; %Tref
T2 = 298.15 + 10;
Tref = 298.15;
T3 = 298.15 - 10;
T4 = 298.15 + 20;
T5 = 298.15 - 20;

R1 = exp(-1 * (Eap/R) * ((1/T1) - (1/Tref)) );
R2 = exp(-1 * (Eap/R) * ((1/T2) - (1/Tref)) );
R3 = exp(-1 * (Eap/R) * ((1/T3) - (1/Tref)) );
R4 = exp(-1 * (Eap/R) * ((1/T4) - (1/Tref)) );
R5 = exp(-1 * (Eap/R) * ((1/T5) - (1/Tref)) );

Q10a = (R2/R1) ^ (10/(T2-T1))
Q10b = (R2/R3) ^ (10/(T2-T3))
Q10c = (R3/R1) ^ (10/(T3-T1))
Q10d = (R4/R1) ^ (10/(T4-T1))
Q10e = (R5/R1) ^ (10/(T5-T1))