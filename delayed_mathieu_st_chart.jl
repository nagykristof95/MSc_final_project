using Printf
using Plots
using Statistics
using LinearAlgebra
using MDBM
using PyPlot
pygui(true);
using Interpolations

dim=2 #DoF of the system (after Cauchy transcription)
mult=3 #Ns=dim*mult
n=20 #timestep number
#choosing numerical simulation method
#RK4: 4th order Runge-Kutta
method="RK4"
i=50 #number of iteration

#parameters of the delayed Mathieu equation
tau=2*pi
omega=1
kappa=0.1
delta=3.2
epsilon=1
b=0.5
T=2*pi/omega

#timestep
dt=tau/(n-1)

#matrix coefficients
function A(t)
 [-kappa -(delta+epsilon*cos(omega*t)); 1 0]
end
B=[0 b; 0 0]

#checking the relation of tau and T
function checkT(tau,T)
    if tau <= T
        println("OK")
    else
        println("Error")
    end
end

checkT(tau,T)
k=Int16(floor(T/tau)) #k factor


#initial function ()
function v_ini(t)
    transpose([1.5*sin(omega*(t+tau/2)), 2.3*sin(omega*(t+tau/3)), -3*sin(omega*(t+tau/4)), 4*sin(omega*(t+tau/5)), 5*sin(omega*(t+tau/6)), 6*sin(omega*(t+tau/7))])
end
vect0=zeros(n,mult*dim)
for j = 0:n-1
        global vect0[j+1,:]=v_ini(-tau+j*dt)
end
vect0
vect1=rand(n,mult*dim)*(-10)

#interpolate
function it(A)
    matrdim=size(A,2)-1
    scaleitp=A[:,1]
    knots=(scaleitp,)
    im=[interpolate(knots, A[:,2],Gridded(Linear()))]
    for j = 3:matrdim+1
        imtemp=[interpolate(knots, A[:,j], Gridded(Linear()))]
        im=vcat(im,imtemp)
    end
    return(im)
end

#substitution in interpolating function
function sub(it,t)
    subdim=size(it,1)
    out=zeros(subdim)
    for j = 1:subdim
        out[j]=it[j](t)
    end
    return(out)
end

function simulation(init)
        sol00=init
        ev=zeros(mult*dim,mult*dim)
        eval=zeros(ComplexF64,mult*dim,1)
        for g = 1:i
        sol=zeros(n-1,mult*dim)
        tvect0=collect((-tau:dt:k*tau))
        sol0=hcat(tvect0,vcat(sol00,sol))
            for m = 2:dim:(dim*mult-(dim-2))
            sol0m=hcat(sol0[:,1],sol0[:,m:m+dim-1])
                for j = 0:k*(n-2)
                    sol0mit=it(sol0m)

                    Y1=sub(sol0mit,j*dt)
                    Y2=Y1+(dt/2)*(A(j*dt)*Y1+B*sub(sol0mit,-tau+j*dt))
                    Y3=Y1+(dt/2)*(A((j+0.5)*dt)*Y2+B*sub(sol0mit,-tau+(j+0.5)*dt))
                    Y4=Y1+dt*(A((j+0.5)*dt)*Y3+B*sub(sol0mit,-tau+(j+0.5)*dt))
                    Y=Y1+(dt/6)*((A(j*dt)*Y1+B*sub(sol0mit,-tau+j*dt))+2*(A((j+0.5)*dt)*Y2+B*sub(sol0mit,-tau+(j+0.5)*dt))+2*(A((j+0.5)*dt)*Y3+B*sub(sol0mit,-tau+(j+0.5)*dt))+(A((j+1)*dt)*Y4+B*sub(sol0mit,-tau+(j+1)*dt)))

                    sol0m[n+j+1,2:2+(dim-1)]=transpose(Y)
                end
                sol0[:,m:m+dim-1]=sol0m[:,2:2+(dim-1)]
            end
            #itt van egy eredmény, itt kéne a mátrixos számítást megcsinálni
            #sol00=sol0[k*n:n+k*n-1,2:mult*dim+1]
            #mc=vcat(mc,sol00[2:n,:])
            S=sol0[1:n,2:mult*dim+1]
            V=sol0[k*n:n+k*n-1,2:mult*dim+1]
            H=pinv(S)*V
            Am=eigvecs(H)
            sol00comp=zeros(ComplexF64,n,mult*dim)
            for o = 1:n
                sol00comp[o,:]=Am*V[o,:]
            end
            sol00=sol00comp
            ev=hcat(ev,Am)
            evm=eigvals(H)
            eval=hcat(eval,evm)
        end
        return(eval[:,2:i+1])
end

function ploteigvals(ev)
    plotev=zeros(size(ev,1),2*i)
    for j=1:i
        plotev[:,2j-1]=[real(ev[n,j]) for n in 1:size(ev,1)]
        plotev[:,2j]=[imag(ev[n,j]) for n in 1:size(ev,1)]
    end
    return(plotev)
end

aaa=simulation(vect1)[:,47:50]
size(aaa,1)

ploteigvals(aaa)
aaaa=simulation(vect0)[:,8:10]
norm(eigvecs([1 2; 3 4])[:,2])

eigRe1=[real(U[n,1]) for n in 1:size(U,1)]
eigIm1=[imag(U[n,1]) for n in 1:size(U,1)]
[1 2; 3 4]*[1; 2]

pinv([1 1; 2 2])
