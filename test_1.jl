using Printf
using Plots
using Statistics
using LinearAlgebra

#machine parameters
p=8 #order of matrix power calculation
m=14 #discretization number
mx=0.01986
my=0.02008
kx=1.60312
ky=1.155697
sx=408866
sy=413445

#cutting parameters
kt=644*10^6
kr=0.368
d=0.008
ae=0.0004
z=1

#calculating cutting teeth
aepu=ae/d
aepd=2-ae/d
phi_inu=0; phi_outu=acos(1-2*aepu);
phi_ind=acos(1-2(aepd-1)); phi_outd=pi;

function phi(j,t,Omega)
    mod(t*Omega+j*2*pi/z,2*pi)
end

function gu(phi::Float64)
    gu_val::Int8=0;
    if phi < phi_inu
        gu_val=0
    elseif phi >= phi_inu && phi <= phi_outu
        gu_val=1
    elseif phi > phi_outu
        gu_val=0
    else
end
return(gu_val)
end

function gd(phi::Float64)
    gd_val::Int8=0;
    if phi < phi_ind
        gd_val=0
    elseif phi >= phi_ind && phi <= phi_outd
        gd_val=1
    elseif phi > phi_outd
        gd_val=0
    else
end
return (gd_val)
end

#delay magnitude
function tau(Omega)
    2*pi/(Omega*z)
end
#timestep
function deltat(Omega)
    tau(Omega)/m
end

#how to plot
x=[n for n in 0:0.2:2*pi]
y=[gu(n) for n in 0:0.2:2*pi]
println(y)
pplot=plot(x,y,seriestype=:scatter)

#tool coefficient matrices
Mm=[mx 0; 0 my];
Km=[kx 0; 0 ky];
Sm=[sx 0; 0 sy];
println(Mm," ",Km ," ",Sm)


function W(t,Omega,w)
    W_f=zeros(Float64,2,2)
    for j = 0:(z-1)
        W_f=W_f+gu(phi(j,t,Omega))*[-sin(2*phi(j,t,Omega))-kr+kr*cos(2*phi(j,t,Omega)) -1-cos(2*phi(j,t,Omega))-kr*sin(2*phi(j,t,Omega)); 1-cos(2*phi(j,t,Omega))-kr*sin(2*phi(j,t,Omega)) sin(2*phi(j,t,Omega))-kr-kr*cos(2*phi(j,t,Omega))]
    end
    return ((w*kt/2)*W_f)
end

println(W(0.001,200.0,100.0))

function A(t,Omega,w)
    vcat(hcat(-inv(Mm)*Km,-inv(Mm)*(Sm-W(t,Omega,w))),hcat(Matrix{Float64}(I, 2, 2),zeros(2,2)))
end

function B(t,Omega,w)
    vcat(hcat(zeros(2,2),-inv(Mm)*W(t,Omega,w)),hcat(zeros(2,2),zeros(2,2)))
end

function expM(M)
n=size(M)[1]
Mval=zeros(n,n)
    for j = 0:p
        Mval=Mval+(M^j)/(factorial(j))
    end
return(Mval)
end

function P(i,Omega,w)
    expM(A(i*deltat(Omega),Omega,w)*deltat(Omega))
end

function R(i,Omega,w)
    0.5*((expM(A(i*deltat(Omega),Omega,w)*deltat(Omega))-Matrix{Float64}(I, 4, 4))*inv(A(i*deltat(Omega),Omega,w))*B(i*deltat(Omega),Omega,w))
end

function Zw(i,Omega,w)
    vcat(hcat(P(i,Omega,w),zeros(4,2*(m+2)-8),R(i,Omega,w)[:,3:4],R(i,Omega,w)[:,3:4]),hcat(zeros(2*(m+2)-4,2),Matrix{Float64}(I, 2*(m+2)-4, 2*(m+2)-4),zeros(2*(m+2)-4,2)))
end

function Zws(Omega,w)
    Z=Zw(0,Omega,w)
    for j = 1:(m-1)
        Z=Zw(j,Omega,w)*Z
    end
    return(Z)
end

eigv=eigvals(Zws(1700,0.001))

eigRe=[real(eigv[n]) for n in 1:length(eigv)]
eigIm=[imag(eigv[n]) for n in 1:length(eigv)]
ploteig=plot(eigRe,eigIm,seriestype=:scatter)
Re_circle=[cos(n) for n in 0:0.05:2*pi]
Im_circle=[sin(n) for n in 0:0.05:2*pi]
plot!(ploteig,Re_circle,Im_circle)
display(eigv)
display(ploteig)
