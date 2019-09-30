using Printf
using Statistics
using LinearAlgebra
using MDBM
using PyPlot
pygui(true);
using Interpolations
using Plots

dim=2 #DoF of the system (after Cauchy transcription)
mult=3 #Ns=dim*mult
n=40 #timestep number
#choosing numerical simulation method
#RK4: 4th order Runge-Kutta
method="RK4"
i=20 #number of iteration
iT=8

#parameters of the delayed Mathieu equation
tau=2*pi
omega=1
kappa=0.1
delta=1.5
epsilon=1
b=0.25
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
vect1=rand(n,mult*dim)*(-10)+ones(n,mult*dim)*5

#interpolate
function it(A)
    matrdim=size(A,2)-1
    scaleitp=real(A[:,1])
    knots=(scaleitp,)
    ARe=real(A)
    AIm=imag(A)

    imRe=[interpolate(knots, ARe[:,2],Gridded(Linear()))]
    for j = 3:matrdim+1
        imRetemp=[interpolate(knots, real(ARe[:,j]), Gridded(Linear()))]
        imRe=vcat(imRe,imRetemp)
    end
    imIm=[interpolate(knots, AIm[:,2],Gridded(Linear()))]
    for j = 3:matrdim+1
        imImtemp=[interpolate(knots, real(AIm[:,j]), Gridded(Linear()))]
        imIm=vcat(imIm,imImtemp)
    end
    return(hcat(imRe,imIm))
end

#substitution in interpolating function
function sub(it,t)
    subdim=size(it,1)
    out=zeros(ComplexF64,subdim)
    for j = 1:subdim
        out[j]=it[j,1](t)+it[j,2](t)*im
    end
    return(out)
end
collect(-tau:dt:iT*tau)
function simulation(init)
        solreturn=zeros(ComplexF64,1,1+mult*dim) #matrix for time integration solution (test)
        eval=zeros(ComplexF64,mult*dim,1) #eigenvalue matrix

        sol00=init #S0 inital random array
        #iteration loop
        for g = 1:i
        sol=zeros(iT*n-iT,mult*dim) #empty solution vector
        tvect0=collect(-tau:dt:iT*tau+dt/2) #discretized time vector
        sol0=hcat(tvect0,vcat(sol00,sol)) #constructing solution array
            #time integration
            for m = 2:dim:(dim*mult-(dim-2)) #choosing "one physical system"
            sol0m=hcat(sol0[:,1],sol0[:,m:m+dim-1])
                #time integration steps
                for j = 0:iT*n-iT-1
                    sol0mit=it(sol0m)

                    Y1=sub(sol0mit,j*dt)
                    Y2=Y1+(dt/2)*(A(j*dt)*Y1+B*sub(sol0mit,-tau+j*dt))
                    Y3=Y1+(dt/2)*(A((j+0.5)*dt)*Y2+B*sub(sol0mit,-tau+(j+0.5)*dt))
                    Y4=Y1+dt*(A((j+0.5)*dt)*Y3+B*sub(sol0mit,-tau+(j+0.5)*dt))
                    Y=Y1+(dt/6)*((A(j*dt)*Y1+B*sub(sol0mit,-tau+j*dt))+2*(A((j+0.5)*dt)*Y2+B*sub(sol0mit,-tau+(j+0.5)*dt))+2*(A((j+0.5)*dt)*Y3+B*sub(sol0mit,-tau+(j+0.5)*dt))+(A((j+1)*dt)*Y4+B*sub(sol0mit,-tau+(j+1)*dt)))

                    sol0m[n+j+1,2:2+(dim-1)]=transpose(Y) #filling solution array of the physical system
                end
                sol0[:,m:m+dim-1]=sol0m[:,2:2+(dim-1)] #filling iterational solution array
            end
            solreturn=vcat(solreturn,sol0) #filling global solution array (test)
            S=sol0[(iT-1)*n-((iT-1)-1):(iT)*n-(iT-1),2:mult*dim+1] #taking array Sj (n x Ns)
            V=sol0[iT*n-(iT-1):(iT+1)*n-iT,2:mult*dim+1] #taking array Vj (n x Ns)
            #H=pinv(S)*V #calculating matrix Hj (Ns x Ns)
            H=pinv(S)*V
            Am=eigvecs(H) #calculation eigenvectors of Hj, creating matrix Am (Ns x Ns)
            sol00=zeros(ComplexF64,n,mult*dim) #calculating matrix Sj+1 (n x Ns)
            for o = 1:n
                sol00[o,:]=Am*V[o,:]
            end
            sol00=sol00
            evm=eigvals(H) #obtaining dominant eigenvalues in the iteration step
            eval=hcat(eval,evm) #filling up the eigenvalue array
        end
        return(eval[:,2:i+1]) #returning the global eigenvalue array
        #return(solreturn)
end

aa=simulation(vect1)

aa[6000:6100,:]

function evplot(aaa,k)
U=aaa[:,k]
eigRe1=[real(U[n,1]) for n in 1:size(U,1)]
eigIm1=[imag(U[n,1]) for n in 1:size(U,1)]
plot1=Plots.plot(eigRe1,eigIm1,seriestype=:scatter)
Re_circle=[cos(n) for n in 0:0.05:2*pi]
Im_circle=[sin(n) for n in 0:0.05:2*pi]
return(Plots.plot!(plot1,Re_circle,Im_circle))
end

evplot(aa,1)
evplot(aa,2)
evplot(aa,3)
evplot(aa,4)
evplot(aa,5)
evplot(aa,6)
evplot(aa,7)
evplot(aa,8)
evplot(aa,9)
evplot(aa,10)
evplot(aa,12)
evplot(aa,11)
evplot(aa,13)
evplot(aa,14)
evplot(aa,15)
evplot(aa,16)
evplot(aa,17)
evplot(aa,18)
evplot(aa,19)
evplot(aa,20)
evplot(aa,21)
evplot(aa,22)
evplot(aa,23)
evplot(aa,24)
evplot(aa,25)
evplot(aa,26)
evplot(aa,27)
evplot(aa,28)
evplot(aa,29)
evplot(aa,30)
evplot(aa,31)
evplot(aa,32)



solt1=real(aa)
v1plot=[solt1[j,2] for j in 1:size(solt1,1)]
v2plot=[solt1[j,3] for j in 1:size(solt1,1)]
v3plot=[solt1[j,4] for j in 1:size(solt1,1)]
v4plot=[solt1[j,5] for j in 1:size(solt1,1)]
v5plot=[solt1[j,6] for j in 1:size(solt1,1)]
v6plot=[solt1[j,7] for j in 1:size(solt1,1)]
Plots.plot(v1plot)
Plots.plot(v2plot)
Plots.plot(v3plot)
Plots.plot(v4plot)
Plots.plot(v5plot)
Plots.plot(v6plot)

function avarr(A)
    dim1=size(A,1)
    dim2=size(A,2)
    sum=0
    for i=1:dim1
            for j=1:dim2
                sum=sum+norm(A[i,j])
            end
    end
    return(sum/(dim1*dim2))
end

avarr([1 2; 3 4])
(1/2.5)*[1 2; 3 4]


norm(-1)
avarr([1 2+4*im; 3 4])
