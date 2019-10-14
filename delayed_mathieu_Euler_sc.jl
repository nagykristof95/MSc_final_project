
using Printf
using Statistics
using LinearAlgebra
using MDBM
using PyPlot
pygui(true);
using Interpolations
using Plots

function evplot(a,k)
U=a[:,k]
eigRe1=[real(U[n,1]) for n in 1:size(U,1)]
eigIm1=[imag(U[n,1]) for n in 1:size(U,1)]
plot1=Plots.plot(eigRe1,eigIm1,seriestype=:scatter)
Re_circle=[cos(n) for n in 0:0.05:2*pi]
Im_circle=[sin(n) for n in 0:0.05:2*pi]
return(Plots.plot!(plot1,Re_circle,Im_circle))
end


function check(a)
    v=false
    s=maximum(size(a))
    dist=zeros(s)
    aRe=real(a)
    aIm=broadcast(abs,imag(a))
    for j=1:s
        if abs(aIm[j])<mar
        else
            for k=1:s
                if k==j
                dist[k]=1
                else
                dist[k]=norm((aRe[j]+im*aIm[j])-(aRe[k]+im*aIm[k]))
                end
            end
            if minimum(dist)<mar
            else v=true
            end
        end
    end
    return(v)
end

function it(A)
    matrdim=size(A,2)-1
    scaleitp=real(A[:,1])
    knots=(scaleitp,)
    ARe=real(A)
    AIm=imag(A)

    imRe=[interpolate(knots, ARe[:,2],Gridded(Linear()))]
    for j = 3:matrdim+1
        imRetemp=[interpolate(knots, ARe[:,j], Gridded(Linear()))]
        imRe=vcat(imRe,imRetemp)
    end
    imIm=[interpolate(knots, AIm[:,2],Gridded(Linear()))]
    for j = 3:matrdim+1
        imImtemp=[interpolate(knots, AIm[:,j], Gridded(Linear()))]
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

function findnormmax(a)
    s=maximum(size(a))
    norma=zeros(s)
    for j = 1:s
        norma[j]=norm(a[j])
    end
    return(a[findmax(norma)[2]])
end

function control(a,r)
    s=maximum(size(a))
    norma=zeros(s)
    for j = 1:s
        norma[j]=norm(a[j])
    end
        max1=findmax(norma)[2]
        norma[max1]=0
        max2=findmax(norma)[2]
        norma[max2]=0
    var=true
    for j=1:s
        if norma[j]<r
        else var=false
        end
    end
    if (abs(imag(a[max1]))<mar && abs(imag(a[max2]))<mar) || abs((real(a[max1])-real(a[max2])))<mar
    else var=false
    end
    return(var)
end



dim=2 #DoF of the system (after Cauchy transcription)
mult=16 #Ns=dim*mult
n=40 #timestep number
mar=0.0001
#choosing numerical simulation method
#RK4: 4th order Runge-Kutta
#EI: Implicit Euler
method="RK4"
i=50 #number of iteration
r=n-1

#parameters of the delayed Mathieu equation
tau=2*pi
omega=1
kappa=0.1
delta=0.5
epsilon=1
b=1.2
T=2*pi/omega

#timestep
dt=tau/(n-1)

#matrix coefficients
function A(t,vv)
    Amatr=zeros(dim,dim)
    #vv[1]=kappa
    #vv[2]=delta
    #vv[3]=epsilon
    #vv[4]=omega
    Amatr=[-vv[1] -(vv[2]+vv[3]*cos(vv[4]*t)); 1 0]
    return(Amatr)
end

function B(bv)
    Bmatr=zeros(dim*dim)
    Bmatr=[0 bv; 0 0]
    return(Bmatr)
end

function simulation(kappa1,delta1,epsilon1,omega1,b1)
        v=[kappa1, delta1, epsilon1, omega1]
        solreturn=zeros(ComplexF64,1,1+mult*dim) #matrix for time integration solution (test)
        eval=zeros(ComplexF64,mult,1)
        S1=zeros(ComplexF64,n,mult*dim)
        V1=zeros(ComplexF64,n,mult*dim)
        S=zeros(ComplexF64,n*dim,mult)
        V=zeros(ComplexF64,n*dim,mult)
        Vj=zeros(ComplexF64,n*dim,mult)
        Vjnorm=zeros(ComplexF64,n*dim,mult)
        H=zeros(ComplexF64,mult,mult)
        tvect0=collect(-tau:dt:tau+dt/2) #discretized time vector
        sol=zeros(ComplexF64,(n-1),mult*dim) #empty solution vector
        sol00=rand(n,mult*dim)*(-2)+2*ones(n,mult*dim) #S0 inital random array
        g=1
        #iteration loop
        for g = 1:i
        sol0=hcat(tvect0,vcat(sol00,sol)) #constructing solution array
            #time integration
            for m = 2:dim:(dim*mult-(dim-2)) #choosing "one physical system"
            sol0m=hcat(sol0[:,1],sol0[:,m:m+dim-1])
                #time integration steps
            if method == "EE"
                for j = 0:n-2

                    Y=sol0m[n+j,2:2+(dim-1)]+dt*(A(sol0m[n+j,1],v)*sol0m[n+j,2:2+(dim-1)]+B(b1)*sol0m[1+j,2:2+(dim-1)])
                    sol0m[n+j+1,2:2+(dim-1)]=transpose(Y) #filling solution array of the physical system

                end

            elseif method == "RK4_int"
                for j = 0:n-2
                    sol0mit=it(sol0m)

                    Y1=sub(sol0mit,j*dt)
                    Y2=Y1+(dt/2)*(A(j*dt,v)*Y1+B(b1)*sub(sol0mit,-tau+j*dt))
                    Y3=Y1+(dt/2)*(A((j+0.5)*dt,v)*Y2+B(b1)*sub(sol0mit,-tau+(j+0.5)*dt))
                    Y4=Y1+dt*(A((j+0.5)*dt,v)*Y3+B(b1)*sub(sol0mit,-tau+(j+0.5)*dt))

                    Y=Y1+(dt/6)*((A(j*dt,v)*Y1+B(b1)*sub(sol0mit,-tau+j*dt))+2*(A((j+0.5)*dt,v)*Y2+B(b1)*sub(sol0mit,-tau+(j+0.5)*dt))+2*(A((j+0.5)*dt,v)*Y3+B(b1)*sub(sol0mit,-tau+(j+0.5)*dt))+(A((j+1)*dt,v)*Y4+B(b1)*sub(sol0mit,-tau+(j+1)*dt)))

                    sol0m[n+j+1,2:2+(dim-1)]=transpose(Y) #filling solution array of the physical system

                end

            elseif method == "RK4"
                for j = 0:n-2

                    Ydtau=0.5*(sol0m[j+1,2:2+(dim-1)]+sol0m[j+2,2:2+(dim-1)])

                    Y1=sol0m[n+j,2:2+(dim-1)]
                    Y2=Y1+(dt/2)*(A(j*dt,v)*Y1+B(b1)*sol0m[j+1,2:2+(dim-1)])
                    Y3=Y1+(dt/2)*(A((j+0.5)*dt,v)*Y2+B(b1)*Ydtau)
                    Y4=Y1+dt*(A((j+0.5)*dt,v)*Y3+B(b1)*Ydtau)

                    Y=Y1+(dt/6)*((A(j*dt,v)*Y1+B(b1)*sol0m[j+1,2:2+(dim-1)])+2*(A((j+0.5)*dt,v)*Y2+B(b1)*Ydtau)+2*(A((j+0.5)*dt,v)*Y3+B(b1)*sol0m[j+1,2:2+(dim-1)])+(A((j+1)*dt,v)*Y4+B(b1)*sol0m[j+2,2:2+(dim-1)]))

                    sol0m[n+j+1,2:2+(dim-1)]=transpose(Y) #filling solution array of the physical system

                end
            end
                sol0[:,m:m+dim-1]=sol0m[:,2:2+(dim-1)] #filling iterational solution array
            end
            #solreturn=vcat(solreturn,sol0) #filling global solution array (test)
            S1=sol0[1:n,2:mult*dim+1]
            V1=sol0[n:2*n-1,2:mult*dim+1]
                for p=1:n
                    for s=1:mult
                        for q=1:dim
                            S[1+(p-1)*dim+(q-1),s]=S1[p,1+(s-1)*dim+(q-1)]
                            V[1+(p-1)*dim+(q-1),s]=V1[p,1+(s-1)*dim+(q-1)]
                        end
                    end
                end
            H=pinv(S)*V
            Am=eigvecs(H) #calculation eigenvectors of Hj, creating matrix Am (Ns x Ns)
            for o = 1:n*dim
                Vj[o,:]=Am*V[o,:]
            end
            for h = 1:mult
                Vjnorm[:,h]=normalize(Vj[:,h])
            end
            sol00=zeros(ComplexF64,n,mult*dim)
            for p=1:n
                for s=1:mult
                    for q=1:dim
                             sol00[p,1+(s-1)*dim+(q-1)]=Vjnorm[1+dim*(p-1)+(q-1),s]
                    end
                end
            end
            #sol00=sol1
            evm=eigvals(H) #obtaining dominant eigenvalues in the iteration step
            eval=hcat(eval,evm)      #filling up the eigenvalue array
            #if (control(evm,mar) && g>3) break  end
            if check(evm) print(g-1); break end
        end
        #evm=eigvals(H)
         return(norm(findnormmax(eval[:,end-1]))-1)
        #return(norm(findnormmax(evm))-1)
        #return(eval[:,2:end]) #returning the global eigenvalue array
        #return(solreturn)
end

Sim=simulation(kappa,delta,epsilon,omega,b)

evplot(Sim,1)
evplot(Sim,2)
evplot(Sim,3)
evplot(Sim,4)
evplot(Sim,5)
evplot(Sim,6)
evplot(Sim,7)
evplot(Sim,8)
evplot(Sim,10)
evplot(Sim,20)
evplot(Sim,15)
evplot(Sim,20)
evplot(Sim,30)
evplot(Sim,100)
evplot(Sim,200)
evplot(Sim,400)

check(Sim[:,4])


function foo(x,y)
    v=[kappa, x, epsilon, omega]
    return(simulation(v[1],v[2],v[3],v[4],y))
end

ax1=Axis([-1.0,0,1,2,3,4,5],"delta") # initial grid in x direction
ax2=Axis(-1.5:0.5:1.5,"b") # initial grid in y direction

mymdbm=MDBM_Problem(foo,[ax1,ax2])
iteration=4 #number of refinements (resolution doubling)
@time solve!(mymdbm,iteration)

x_eval,y_eval=getevaluatedpoints(mymdbm)
x_sol,y_sol=getinterpolatedsolution(mymdbm)
fig = figure(1);clf()
PyPlot.scatter(x_eval,y_eval,s=5)
PyPlot.scatter(x_sol,y_sol,s=5)

[1,2,5,-4,5][end-1]
