using Printf
using Plots
using Statistics
using LinearAlgebra

dim=2
mult=3
n=30

i=5

tau=2*pi
omega=1

#Delayed Mathieu parameters
kappa=0.1
delta=3.2
epsilon=1
b=0.5

#Cauchy matrices
function A(t)
 [-kappa -(delta+epsilon*cos(omega*t)); 1 0]
end

B=[0 b; 0 0]

deltat=tau/n
T=2*pi/omega

function checkT(tau,T)
    if tau <= T
        println("OK")
    else
        println("Error")
    end
end

checkT(tau,T)

#initial function
function v_ini(t)
    transpose([1*sin(omega*(t+tau/2)), 2*sin(omega*(t+tau/3)), 3*sin(omega*(t+tau/4)), 4*sin(omega*(t+tau/5)), 5*sin(omega*(t+tau/6)), 6*sin(omega*(t+tau/7))])
end

vect=zeros(n+1,mult*dim)
for j = 0:n
        vect[j+1,:]=v_ini(-tau+j*deltat)
end
vect
function simulation(init)
    k=Int16(floor(i*T/deltat))
    solmatrix=zeros(k,mult*dim)
    sol=vcat(init,solmatrix)

    for m = 1:dim:(dim*mult-(dim-1))
        for j = 0:(k-1)
            Y1=sol[n+1+j,m:m+(dim-1)]
            Y2=Y1+deltat/2*(A(j*deltat)*Y1+B*sol[j+1,m:m+(dim-1)])
            Y3=Y1+deltat/2*(A((j+1/2)*deltat)*Y2+B*(sol[j+1,m:m+(dim-1)]+sol[j+2,m:m+(dim-1)])/2)
            Y4=Y1+deltat*(A((j+1/2)*deltat)*Y3+B*(sol[j+1,m:m+(dim-1)]+sol[j+2,m:m+(dim-1)])/2)
            Y=Y1+deltat/6*((A(j*deltat)*Y1+B*sol[j+1,m:m+(dim-1)])+2*(A((j+1/2)*deltat)*Y2+B*(sol[j+1,m:m+(dim-1)]+sol[j+2,m:m+(dim-1)])/2)+2*(A((j+1/2)*deltat)*Y3+B*(sol[j+1,m:m+(dim-1)]+sol[j+2,m:m+(dim-1)])/2)+(A((j+1)*deltat)*Y4+B*sol[j+2,m:m+(dim-1)]))
            sol[n+2+j,m:m+(dim-1)]=transpose(Y)
        end
    end
    tvect=transpose(zeros(1,n+1+k))
    for j = 1:(n+1+k)
        tvect[j]=-n*deltat+j*deltat
    end
    solt=hcat(tvect,sol)
    return(solt)
end

solt1=simulation(vect)
@btime simulation(vect)
@time simulation(vect)

tplot=[solt1[j,1] for j in 1:size(solt1,1)]
v1plot=[solt1[j,2] for j in 1:size(solt1,1)]
v2plot=[solt1[j,3] for j in 1:size(solt1,1)]
v3plot=[solt1[j,4] for j in 1:size(solt1,1)]
v4plot=[solt1[j,5] for j in 1:size(solt1,1)]
v5plot=[solt1[j,6] for j in 1:size(solt1,1)]
v6plot=[solt1[j,7] for j in 1:size(solt1,1)]
plot1=plot(tplot,hcat(v1plot,v2plot,v3plot,v4plot,v5plot,v6plot))

function Ha(solt)
    rat=Int16(floor(T/tau))
    ha=zeros(n,(i+1)*mult*dim)
    for j = 1:i+1
        ha[:,(j-1)*mult*dim+1:j*mult*dim]=solt[(j-1)*n+1:j*n,2:mult*dim+1]
    end
    return(ha)
end

ha=Ha(solt1)

function eigenv(haa)
    eigena=zeros(mult*dim,1)
    Hglob=zeros(mult*dim,i*mult*dim)
    for j = 1:i
        Hglob[:,(j-1)*mult*dim+1:j*mult*dim]=inv(transpose(haa[:,(j-1)*mult*dim+1:j*mult*dim])*haa[:,(j-1)*mult*dim+1:j*mult*dim])*(transpose(haa[:,(j-1)*mult*dim+1:j*mult*dim])*haa[:,j*mult*dim+1:(j+1)*mult*dim])
    end
    for j = 1:i
        eigena=hcat(eigena,eigvals(Hglob[:,(j-1)*mult*dim+1:j*mult*dim]))
    end
    eigena=eigena[:,2:i+1]
    return(eigena)
end

U=eigenv(ha)

eigRe1=[real(U[n,1]) for n in 1:size(U,1)]
eigIm1=[imag(U[n,1]) for n in 1:size(U,1)]
plot1=plot(eigRe1,eigIm1,seriestype=:scatter)
Re_circle=[cos(n) for n in 0:0.05:2*pi]
Im_circle=[sin(n) for n in 0:0.05:2*pi]
plot!(plot1,Re_circle,Im_circle)


eigRe2=[real(U[n,2]) for n in 1:size(U,1)]
eigIm2=[imag(U[n,2]) for n in 1:size(U,1)]
plot2=plot(eigRe2,eigIm2,seriestype=:scatter)
plot!(plot2,Re_circle,Im_circle)

eigRe3=[real(U[n,3]) for n in 1:size(U,1)]
eigIm3=[imag(U[n,3]) for n in 1:size(U,1)]
plot3=plot(eigRe3,eigIm3,seriestype=:scatter)
plot!(plot3,Re_circle,Im_circle)

eigRe4=[real(U[n,4]) for n in 1:size(U,1)]
eigIm4=[imag(U[n,4]) for n in 1:size(U,1)]
plot4=plot(eigRe4,eigIm4,seriestype=:scatter)
plot!(plot4,Re_circle,Im_circle)

eigRe5=[real(U[n,5]) for n in 1:size(U,1)]
eigIm5=[imag(U[n,5]) for n in 1:size(U,1)]
plot5=plot(eigRe5,eigIm5,seriestype=:scatter)
plot!(plot5,Re_circle,Im_circle)
