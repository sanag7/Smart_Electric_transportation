syms d2 d3 V3

Y = [-1i*20 1i*10 1i*10; 1i*10 -1i*20 1i*10; 1i*10 1i*10 -1i*20];


IY = inv(Y);

P2 = 1.05*10*sin(d2)+1.05*V3*10*sin(d2-d3);

P3 = V3*10*sin(d3)+V3*1.05*10*sin(d3-d2);

Q3 = -(V3*10*cos(d3)+V3*1.05*10*cos(d3-d2)+V3^2*(-20));


S2 = 0.6661;

R = [P2; P3; Q3];

JR = jacobian(R, [d2,d3,V3]);

X = [d2 d3 V3];

Dx = [1; 1; 1];

x = [0; 0; 1];

inter = 0;

while max(abs(Dx))>=.0000001 && inter<10
    
    d2 = x(1);
    d3 = x(2);
    V3 = x(3);
    PL3 = 2.8653;
    QL3 = 1.2244;
    inter = inter + 1;
    
    F = [10.5*sin(d2)+10.5*V3*sin(d2-d3)-0.6661 
        V3*10*sin(d3)+V3*10.5*sin(d3-d2)+PL3 
        -V3*10*cos(d3)-V3*10.5*cos(d2-d3)+V3*V3*20+QL3];
    
    DC = -F;
    
    J = [10.5*cos(d2)+10.5*V3*cos(d2-d3) -10.5*V3*cos(d2-d3) 10.5*sin(d2-d3);
        -10.5*V3*cos(d3-d2) V3*10*cos(d2)+V3*10.5*cos(d3-d2) 10*sin(d3)+10.5*cos(d3-d2);
        -V3*10.5*sin(d3-d2) V3*10*sin(d3)+V3*10.5*sin(d3-d2) -10*cos(d3)-10.5*cos(d3-d2)+2*(20)*V3];
    
    Dx = J\DC;
    
    x=x+Dx;
end

PG1 = 10.5*sin(-d2)+ V3*10*sin(-d3);

QG1 = -(-20+10.5*cos(d2)+V3*10*cos(d3));

QG2 = -(10.5*cos(d2)+1.05^2*(-20)+1.05*V3*10*cos(d2-d3));

I12 = abs((1-1.05*exp(1j*d2))*(-1j*10));

I13 = abs((1-V3*exp(1j*d3))*(-1j*10));

I23 = abs((1.05*exp(1j*d2)-V3*exp(1j*d3))*(-1j*10));

fprintf("when PL3=%g and QL3=%g,\n\nPG1 =%g, \nQG1 =%g,\nd2 =%g,\nQG2 =%g,\nd3 =%g, \n|V3| =%g,\n|i12|=%g, \n|i13|=%g,\n|i23|=%g\n\n",PL3,QL3,PG1,QG1,(180/pi)*d2,QG2,(180/pi)*d3,V3,I12,I13,I23)
        
Dx = [1; 1; 1];

x = [0; 0; 1];

inter = 0;

while max(abs(Dx))>=.00000001 && inter<10
    
    d2c = x(1);
    d3c = x(2);
    V3c = x(3);
    PL3change = 0;
    QL3change = 0.3;
    PL3c = 2.8653+PL3change;
    QL3c = 1.2244+QL3change;
    inter = inter + 1;
    
    F = [10.5*sin(d2c)+10.5*V3c*sin(d2c-d3c)-0.6661 
        V3c*10*sin(d3c)+V3c*10.5*sin(d3c-d2c)+PL3c 
        -V3c*10*cos(d3c)-V3c*10.5*cos(d2c-d3c)+V3c*V3c*20+QL3c];
    
    DC = -F;
    
    J = [10.5*cos(d2c)+10.5*V3c*cos(d2c-d3c) -10.5*V3c*cos(d2c-d3c) 10.5*sin(d2c-d3c);
        -10.5*V3c*cos(d3c-d2c) V3c*10*cos(d2c)+V3c*10.5*cos(d3c-d2c) 10*sin(d3c)+10.5*cos(d3c-d2c);
        -V3c*10.5*sin(d3c-d2c) V3c*10*sin(d3c)+V3c*10.5*sin(d3c-d2c) -10*cos(d3c)-10.5*cos(d3c-d2c)+2*(20)*V3];
    
    Dx = J\DC;
    
    x=x+Dx;
end

PG1c = 10.5*sin(-d2c)+ V3c*10*sin(-d3c);

QG1c = -(-20+10.5*cos(d2c)+V3c*10*cos(d3c));

QG2c = -(10.5*cos(d2c)+1.05^2*(-20)+1.05*V3c*10*cos(d2c-d3c));

I12c = abs((1-1.05*exp(1j*d2c))*(-1j*10));

I13c = abs((1-V3c*exp(1j*d3c))*(-1j*10));

I23c = abs((1.05*exp(1j*d2c)-V3c*exp(1j*d3))*(-1j*10));

crPG1=100*(PG1c-PG1)/PG1;

crQG1=100*(QG1c-QG1)/QG1;

crd2=100*(d2c-d2)/d2;
crQG2=100*(QG2c-QG2)/QG2;
crV3 = 100*(V3c-V3)/V3;
crd3 = 100*(d3c-d3)/d3;
crI12=100*(I12c-I12)/I12;
crI13=100*(I13c-I13)/I13;
crI23=100*(I23c-I23)/I23;


fprintf("when PL3=%g and QL3=%g,\nwhich means PL3change=%g and QL3change=%g \n\nPG1 =%g, %g percent\nQG1 =%g, %g percent\nd2 =%g,%g percent\nQG2  =%g, %g percent\nd3  =%g, %g percent\n|V3| =%g,%g percent\n|i12|=%g, %g percent\n|i13|=%g, %g percent\n|i23|=%g, %gpercent\n\n",PL3c,QL3c,PL3change,QL3change,PG1c,crPG1,QG1c,crQG1,(180/pi)*d2c,crd2,QG2c,crQG2,(180/pi)*d3c,crd3,V3c,crV3,I12c,crI12,I13c,crI13,I23c,crI23)
