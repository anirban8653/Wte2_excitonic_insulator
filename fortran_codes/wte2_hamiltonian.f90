Module wte2_hamiltonian

Implicit None

Contains

Function hamiltonian(kx,ky)

Complex*16, Dimension(8,8):: hamiltonian
complex*16, parameter :: iota = (0.D0, 1.D0)
real*8, parameter :: mup = -1.65 !this is for semimetallic case. Put -1.75 for QSHI case.
real*8, parameter :: mud = 0.74
real*8, parameter :: mu = 0.08
Real*8,Intent(In)::kx,ky



hamiltonian(1,1)=((mud,0.)-0.41*Cos(3.477*kx)-0.41*Cos(3.477*kx)-  &
0.008*Sin(3.477*kx)-0.008*Sin(3.477*kx))-mu
hamiltonian(1,2)=(0.,0.)-(0.,0.012)*exp(iota*(0.+2.43711*ky))-  &
(0.,0.012)*exp((0.,-6.249)*ky+(0,1)*(0.+2.43711*ky))-  &
(0.,0.28)*exp(iota*(0.+2.43711*ky))*Sin(3.477*kx)
hamiltonian(1,3)=(0.,0.)+0.51*exp((0.,-6.249)*ky+(0,1)*(-1.7385*kx+  &
3.99936*ky))*(1.+exp((0.,3.477)*kx))
hamiltonian(1,4)=(0.,0.)+0.39*exp(iota*(-1.7385*kx+1.56225*ky))*(1.-  &
exp((0.,3.477)*kx))+0.29*exp(iota*(-1.7385*kx+  &
1.56225*ky))*(exp((0.,-3.477)*kx)-exp((0.,6.954)*kx))
hamiltonian(1,5)=(0.,0.)+(0.,0.031)*Sin(3.477*kx)+  &
(0.,0.031)*Sin(3.477*kx)
hamiltonian(1,6)=(0.,0.)-(0.051,0.)*exp(iota*(0.+2.43711*ky))-  &
(0.05,0.)*exp((0.,-6.249)*ky+(0,1)*(0.+2.43711*ky))
hamiltonian(1,7)=(0.,0.)
hamiltonian(1,8)=(0.,0.)-(0.011,0.)*exp(iota*(-1.7385*kx+1.56225*ky))*(1+  &
exp((0.,3.477)*kx))

hamiltonian(2,1)=(0.,0.)+(0.,0.012)/exp(iota*(0.+2.43711*ky))+  &
(0.,0.012)*exp((0.,6.249)*ky-(0,1)*(0.+  &
2.43711*ky))+  &
((0.,0.28)*Sin(3.477*kx))/exp(iota*(0.+  &
2.43711*ky))
hamiltonian(2,2)=((mup,0.)+1.13*Cos(3.477*kx)+0.13*Cos(6.249*ky)+  &
1.13*Cos(3.477*kx)+0.13*Cos(6.249*ky)-  &
0.01*Sin(3.477*kx)-0.01*Sin(3.477*kx))-mu
hamiltonian(2,3)=(0.,0.)-0.39*exp(iota*(-1.7385*kx+1.56225*ky))*(1.-  &
exp((0.,3.477)*kx))-0.29*exp(iota*(-1.7385*kx+  &
1.56225*ky))*(exp((0.,-3.477)*kx)-exp((0.,6.954)*kx))
hamiltonian(2,4)=(0.,0.)+0.4*exp(iota*(-1.7385*kx-0.87486*ky))*(1.+  &
exp((0.,3.477)*kx))
hamiltonian(2,5)=(0.,0.)+(0.051,0.)/exp(iota*(0.+2.43711*ky))+  &
(0.05,0.)*exp((0.,6.249)*ky-(0,1)*(0.+  &
2.43711*ky))
hamiltonian(2,6)=(0.,0.)+(0.,0.04)*Sin(3.477*kx)+  &
(0.,0.04)*Sin(3.477*kx)
hamiltonian(2,7)=(0.,0.)-(0.011,0.)*exp(iota*(-1.7385*kx+1.56225*ky))*(1+  &
exp((0.,3.477)*kx))
hamiltonian(2,8)=(0.,0.)

hamiltonian(3,1)=(0.,0.)+0.51*exp((0.,6.249)*ky-  &
(0,1)*(-1.7385*kx+3.99936*ky))*(1.+  &
exp((0.,-3.477)*kx))
hamiltonian(3,2)=(0.,0.)+conjg(-0.39*exp(iota*(-1.7385*kx+1.56225*ky))*(1.-  &
exp((0.,3.477)*kx))-0.29*exp(iota*(-1.7385*kx+  &
1.56225*ky))*(exp((0.,-3.477)*kx)-exp((0.,6.954)*kx)))
hamiltonian(3,3)=((mud,0.)-0.41*Cos(3.477*kx)-0.41*Cos(3.477*kx)+  &
0.008*Sin(3.477*kx)+0.008*Sin(3.477*kx))-mu
hamiltonian(3,4)=(0.,0.)+(0.,0.012)/exp(iota*(0.+2.43711*ky))+  &
(0.,0.012)*exp((0.,6.249)*ky-(0,1)*(0.+2.43711*ky))-  &
((0.,0.28)*Sin(3.477*kx))/exp(iota*(0.+2.43711*ky))
hamiltonian(3,5)=(0.,0.)
hamiltonian(3,6)=(0.,0.)+((0.011,0.)*(1+  &
exp((0.,-3.477)*kx)))/exp(iota*(-1.7385*kx  &
+1.56225*ky))
hamiltonian(3,7)=(0.,0.)-(0.,0.031)*Sin(3.477*kx)-  &
(0.,0.031)*Sin(3.477*kx)
hamiltonian(3,8)=(0.,0.)+(0.051,0.)/exp(iota*(0.+2.43711*ky))+  &
(0.05,0.)*exp((0.,6.249)*ky-(0,1)*(0.+2.43711*ky))

hamiltonian(4,1)=(0.,0.)+conjg(0.39*exp(iota*(-1.7385*kx+1.56225*ky))*(1.-  &
exp((0.,3.477)*kx))+0.29*exp(iota*(-1.7385*kx+  &
1.56225*ky))*(exp((0.,-3.477)*kx)-exp((0.,6.954)*kx)))
hamiltonian(4,2)=(0.,0.)+(0.4*(1.+  &
exp((0.,-3.477)*kx)))/exp(iota*(-1.7385*kx  &
-0.87486*ky))
hamiltonian(4,3)=(0.,0.)-(0.,0.012)*exp(iota*(0.+2.43711*ky))-  &
(0.,0.012)*exp((0.,-6.249)*ky+(0,1)*(0.+  &
2.43711*ky))+(0.,0.28)*exp(iota*(0.+  &
2.43711*ky))*Sin(3.477*kx)
hamiltonian(4,4)=((mup,0.)+1.13*Cos(3.477*kx)+0.13*Cos(6.249*ky)+  &
1.13*Cos(3.477*kx)+0.13*Cos(6.249*ky)+  &
0.01*Sin(3.477*kx)+0.01*Sin(3.477*kx))-mu
hamiltonian(4,5)=(0.,0.)+((0.011,0.)*(1+  &
exp((0.,-3.477)*kx)))/exp(iota*(-1.7385*kx  &
+1.56225*ky))
hamiltonian(4,6)=(0.,0.)
hamiltonian(4,7)=(0.,0.)-(0.051,0.)*exp(iota*(0.+2.43711*ky))-  &
(0.05,0.)*exp((0.,-6.249)*ky+(0,1)*(0.+  &
2.43711*ky))
hamiltonian(4,8)=(0.,0.)-(0.,0.04)*Sin(3.477*kx)-  &
(0.,0.04)*Sin(3.477*kx)

hamiltonian(5,1)=(0.,0.)-(0.,0.031)*Sin(3.477*kx)-  &
(0.,0.031)*Sin(3.477*kx)
hamiltonian(5,2)=(0.,0.)+(0.051,0.)*exp(iota*(0.+2.43711*ky))+  &
(0.05,0.)*exp((0.,-6.249)*ky+(0,1)*(0.+2.43711*ky))
hamiltonian(5,3)=(0.,0.)
hamiltonian(5,4)=(0.,0.)+(0.011,0.)*exp(iota*(-1.7385*kx+1.56225*ky))*(1+  &
exp((0.,3.477)*kx))
hamiltonian(5,5)=((mud,0.)-0.41*Cos(3.477*kx)-0.41*Cos(3.477*kx)+  &
0.008*Sin(3.477*kx)+0.008*Sin(3.477*kx))-mu
hamiltonian(5,6)=(0.,0.)+(0.,0.012)*exp(iota*(0.+2.43711*ky))+  &
(0.,0.012)*exp((0.,-6.249)*ky+(0,1)*(0.+2.43711*ky))-  &
(0.,0.28)*exp(iota*(0.+2.43711*ky))*Sin(3.477*kx)
hamiltonian(5,7)=(0.,0.)+0.51*exp((0.,-6.249)*ky+(0,1)*(-1.7385*kx+  &
3.99936*ky))*(1.+exp((0.,3.477)*kx))
hamiltonian(5,8)=(0.,0.)+0.39*exp(iota*(-1.7385*kx+1.56225*ky))*(1.-  &
exp((0.,3.477)*kx))+0.29*exp(iota*(-1.7385*kx+  &
1.56225*ky))*(exp((0.,-3.477)*kx)-exp((0.,6.954)*kx))

hamiltonian(6,1)=(0.,0.)-(0.051,0.)/exp(iota*(0.+2.43711*ky))-  &
(0.05,0.)*exp((0.,6.249)*ky-(0,1)*(0.+  &
2.43711*ky))
hamiltonian(6,2)=(0.,0.)-(0.,0.04)*Sin(3.477*kx)-  &
(0.,0.04)*Sin(3.477*kx)
hamiltonian(6,3)=(0.,0.)+(0.011,0.)*exp(iota*(-1.7385*kx+1.56225*ky))*(1+  &
exp((0.,3.477)*kx))
hamiltonian(6,4)=(0.,0.)
hamiltonian(6,5)=(0.,0.)-(0.,0.012)/exp(iota*(0.+2.43711*ky))-  &
(0.,0.012)*exp((0.,6.249)*ky-(0,1)*(0.+  &
2.43711*ky))+  &
((0.,0.28)*Sin(3.477*kx))/exp(iota*(0.+  &
2.43711*ky))
hamiltonian(6,6)=((mup,0.)+1.13*Cos(3.477*kx)+0.13*Cos(6.249*ky)+  &
1.13*Cos(3.477*kx)+0.13*Cos(6.249*ky)+  &
0.01*Sin(3.477*kx)+0.01*Sin(3.477*kx))-mu
hamiltonian(6,7)=(0.,0.)-0.39*exp(iota*(-1.7385*kx+1.56225*ky))*(1.-  &
exp((0.,3.477)*kx))-0.29*exp(iota*(-1.7385*kx+  &
1.56225*ky))*(exp((0.,-3.477)*kx)-exp((0.,6.954)*kx))
hamiltonian(6,8)=(0.,0.)+0.4*exp(iota*(-1.7385*kx-0.87486*ky))*(1.+  &
exp((0.,3.477)*kx))

hamiltonian(7,1)=(0.,0.)
hamiltonian(7,2)=(0.,0.)-((0.011,0.)*(1+  &
exp((0.,-3.477)*kx)))/exp(iota*(-1.7385*kx  &
+1.56225*ky))
hamiltonian(7,3)=(0.,0.)+(0.,0.031)*Sin(3.477*kx)+  &
(0.,0.031)*Sin(3.477*kx)
hamiltonian(7,4)=(0.,0.)-(0.051,0.)/exp(iota*(0.+2.43711*ky))-  &
(0.05,0.)*exp((0.,6.249)*ky-(0,1)*(0.+2.43711*ky))
hamiltonian(7,5)=(0.,0.)+0.51*exp((0.,6.249)*ky-  &
(0,1)*(-1.7385*kx+3.99936*ky))*(1.+  &
exp((0.,-3.477)*kx))
hamiltonian(7,6)=(0.,0.)+conjg(-0.39*exp(iota*(-1.7385*kx+1.56225*ky))*(1.-  &
exp((0.,3.477)*kx))-0.29*exp(iota*(-1.7385*kx+  &
1.56225*ky))*(exp((0.,-3.477)*kx)-exp((0.,6.954)*kx)))
hamiltonian(7,7)=((mud,0.)-0.41*Cos(3.477*kx)-0.41*Cos(3.477*kx)-  &
0.008*Sin(3.477*kx)-0.008*Sin(3.477*kx))-mu
hamiltonian(7,8)=(0.,0.)-(0.,0.012)/exp(iota*(0.+2.43711*ky))-  &
(0.,0.012)*exp((0.,6.249)*ky-(0,1)*(0.+2.43711*ky))-  &
((0.,0.28)*Sin(3.477*kx))/exp(iota*(0.+2.43711*ky))

hamiltonian(8,1)=(0.,0.)-((0.011,0.)*(1+  &
exp((0.,-3.477)*kx)))/exp(iota*(-1.7385*kx  &
+1.56225*ky))
hamiltonian(8,2)=(0.,0.)
hamiltonian(8,3)=(0.,0.)+(0.051,0.)*exp(iota*(0.+2.43711*ky))+  &
(0.05,0.)*exp((0.,-6.249)*ky+(0,1)*(0.+  &
2.43711*ky))
hamiltonian(8,4)=(0.,0.)+(0.,0.04)*Sin(3.477*kx)+  &
(0.,0.04)*Sin(3.477*kx)
hamiltonian(8,5)=(0.,0.)+conjg(0.39*exp(iota*(-1.7385*kx+1.56225*ky))*(1.-  &
exp((0.,3.477)*kx))+0.29*exp(iota*(-1.7385*kx+  &
1.56225*ky))*(exp((0.,-3.477)*kx)-exp((0.,6.954)*kx)))
hamiltonian(8,6)=(0.,0.)+(0.4*(1.+  &
exp((0.,-3.477)*kx)))/exp(iota*(-1.7385*kx  &
-0.87486*ky))
hamiltonian(8,7)=(0.,0.)+(0.,0.012)*exp(iota*(0.+2.43711*ky))+  &
(0.,0.012)*exp((0.,-6.249)*ky+(0,1)*(0.+  &
2.43711*ky))+(0.,0.28)*exp(iota*(0.+  &
2.43711*ky))*Sin(3.477*kx)
hamiltonian(8,8)=((mup,0.)+1.13*Cos(3.477*kx)+0.13*Cos(6.249*ky)+  &
1.13*Cos(3.477*kx)+0.13*Cos(6.249*ky)-  &
0.01*Sin(3.477*kx)-0.01*Sin(3.477*kx))-mu

Return

End Function

End Module wte2_hamiltonian