using Printf
using Statistics
using LinearAlgebra
using MDBM
using PyPlot
pygui(true);
using Interpolations
using Plots

#creating iteration array
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

#function for finding the greatest eigenvalue
function normmax(a)
    s=maximum(size(a))
    norma=zeros(s)
    for j = 1:s
        norma[j]=norm(a[j])
    end
    return(a[findmax(norma)[2]])
end

#function for finding the smallest eigenvalue
function normmin(a)
    s=maximum(size(a))
    norma=zeros(s)
    for j = 1:s
        norma[j]=norm(a[j])
    end
    return(a[findmin(norma)[2]])
end

#checking the difference between the greatest and the smallest eigenvalue
function checkconv(a,b,err)
    var=false
    amax=norm(normmax(a))
    amin=norm(normmin(a))
    bmax=norm(normmax(b))
    bmin=norm(normmin(b))
    if abs((amax-amin)-(bmax-bmin))/(amax-amin)<err
        var=true
    end
    return(var)
end

#matrix exponential
function expM(M)
    p=8
    n=size(M)[1]
    Mval=zeros(Float64,n,n)
    for j = 0:p
        Mval=Mval+(M^j)/(factorial(j))
    end
    return(Mval)
end

###################### ISIM #######################
mult=20 #Ns=dim*mult
n=30 #timestep number
#EE: Explicit Euler
#RK4: 4th order Runge-Kutta
#RK4_int: 4th order Runge-Kutta with interpolation function
#IE: Implicit Euler
#RK2: 2nd order Runge-Kutta
method="RK4" #choosing numerical simulation method
i=10 #number of iteration

###################### SDM ########################
m=40 #discretization number
#dtSDM=tau/m


#parameters of the delayed Mathieu equation
dim=2 #DoF of the system (after Cauchy transcription)
const tau=2*pi
const kappa=0.1
const delta=0.5
const epsilon=1.0
const b=0.2
const k=1
#const T=k*tau
#const omega=2*pi/T
#timestep
#dt=tau/(n-1)

#matrix coefficients
function A(t,vv)
    Amatr=zeros(Float64,dim,dim)
    #vv[1]=kappa #vv[2]=delta #vv[3]=epsilon #vv[4]=omega
    Amatr=[-vv[1] -(vv[2]+vv[3]*cos(vv[4]*t)); 1 0]
    return(Amatr)
end

function B(t,bv)
    Bmatr=zeros(dim*dim)
    Bmatr=[0 bv; 0 0]
    return(Bmatr)
end

###################### ISIM #######################
function ISIM(kappa1,delta1,epsilon1,b1,tau1,k1)
        dt=tau1/(n-1)
        v=[kappa1, delta1, epsilon1, 2*pi/(k1*tau1)]
        eval=zeros(ComplexF64,mult,1)
        S=zeros(ComplexF64,n*dim,mult)
        V=zeros(ComplexF64,n*dim,mult)
        Vjnorm=zeros(ComplexF64,n*dim,mult)
        tvect0=collect(-tau1:dt:tau1*k1+dt/2) #discretized time vector
        sol=zeros(ComplexF64,k1*(n-1),mult*dim) #empty solution vector
        sol00=rand(n,mult*dim)*(-2)+2*ones(n,mult*dim) #S0 inital random array
        for g = 1:i #iteration
        sol0=hcat(tvect0,vcat(sol00,sol)) #constructing solution array
            #time integration
            for m = 2:dim:(dim*mult-(dim-2)) #choosing "one physical system"
            sol0m=hcat(sol0[:,1],sol0[:,m:m+dim-1])
            #choosing numerical method
            if method == "EE"
                for j = 0:k1*(n-2)
                    Y=sol0m[n+j,2:2+(dim-1)]+dt*(A(sol0m[n+j,1],v)*sol0m[n+j,2:2+(dim-1)]+B(sol0m[n+j,1],b1)*sol0m[1+j,2:2+(dim-1)])
                    sol0m[n+j+1,2:2+(dim-1)]=transpose(Y)
                end
            elseif method == "RK4_int"
                for j = 0:k1*(n-2)
                    sol0mittau=it(sol0m[j+1:j+3,:])

                    Ytau0=sub(sol0mittau,-tau1+j*dt)
                    Ytau05=sub(sol0mittau,-tau1+(j+0.5)*dt)
                    Ytau1=sub(sol0mittau,-tau1+(j+1)*dt)

                    Y1=sol0m[j+n,2:dim+1]
                    Y2=Y1+(dt/2)*(A(j*dt,v)*Y1+B(j*dt,b1)*Ytau0)
                    Y3=Y1+(dt/2)*(A((j+0.5)*dt,v)*Y2+B((j+0.5)*dt,b1)*Ytau05)
                    Y4=Y1+dt*(A((j+0.5)*dt,v)*Y3+B((j+0.5)*dt,b1)*Ytau05)

                    Y=Y1+(dt/6)*((A(j*dt,v)*Y1+B(j*dt,b1)*Ytau0)+2*(A((j+0.5)*dt,v)*Y2+B((j+0.5)*dt,b1)*Ytau05)+2*(A((j+0.5)*dt,v)*Y3+B((j+0.5)*dt,b1)*Ytau05)+(A((j+1)*dt,v)*Y4+B((j+1)*dt,b1)*Ytau1))

                    sol0m[n+j+1,2:2+(dim-1)]=transpose(Y)
                end
            elseif method == "RK4"
                for j = 0:k1*(n-2)
                    Ydtau=0.5*(sol0m[j+1,2:2+(dim-1)]+sol0m[j+2,2:2+(dim-1)])

                    Y1=sol0m[n+j,2:2+(dim-1)]
                    Y2=Y1+(dt/2)*(A(j*dt,v)*Y1+B(j*dt,b1)*sol0m[j+1,2:2+(dim-1)])
                    Y3=Y1+(dt/2)*(A((j+0.5)*dt,v)*Y2+B((j+0.5)*dt,b1)*Ydtau)
                    Y4=Y1+dt*(A((j+0.5)*dt,v)*Y3+B((j+0.5)*dt,b1)*Ydtau)

                    Y=Y1+(dt/6)*((A(j*dt,v)*Y1+B(j*dt,b1)*sol0m[j+1,2:2+(dim-1)])+2*(A((j+0.5)*dt,v)*Y2+B((j+0.5)*dt,b1)*Ydtau)+2*(A((j+0.5)*dt,v)*Y3+B((j+0.5)*dt,b1)*sol0m[j+1,2:2+(dim-1)])+(A((j+1)*dt,v)*Y4+B(j*dt,b1)*sol0m[j+2,2:2+(dim-1)]))

                    sol0m[n+j+1,2:2+(dim-1)]=transpose(Y)
                end
            elseif method == "IE"
                for j = 0:k1*(n-2)
                    X=sol0m[j+n,2:2+(dim-1)]
                    Yt=sol0m[j+n,2:2+(dim-1)]
                    Ytau=sol0m[j+1,2:2+(dim-1)]
                    for l=1:3
                        X=X-inv(dt*A(j*dt,v)-Matrix{Float64}(I,dim,dim))*(Yt+dt*(A((j+1+n)*dt,v)*X+B((j+1+n)*dt,b1)*Ytau)-X)
                    end
                    Y=X
                   sol0m[n+j+1,2:2+(dim-1)]=transpose(Y) #filling solution array of the physical system
                end
            elseif method == "RK2"
                for j = 0:k1*(n-2)

                    Ydtau=0.5*(sol0m[j+1,2:2+(dim-1)]+sol0m[j+2,2:2+(dim-1)])

                    Y1=sol0m[n+j,2:2+(dim-1)]
                    Y2=Y1+(dt/2)*(A((j+0.5)*dt,v)*Y1+B((j+0.5)*dt,b1)*Ydtau)

                    Y=Y1+dt*(A((j+0.5)*dt,v)*Y2+B((j+0.5)*dt,b1)*Ydtau)

                    sol0m[n+j+1,2:2+(dim-1)]=transpose(Y)
                end
            end

                sol0[:,m:m+dim-1]=sol0m[:,2:2+(dim-1)] #filling solution array
            end
            #retrieving the solution of last and second-to-last periods
            S1=sol0[1:n,2:mult*dim+1]
            V1=sol0[end-(n-1):end,2:mult*dim+1]
            #restructuring the solution
                for p=1:n
                    for s=1:mult
                        for q=1:dim
                            S[1+(p-1)*dim+(q-1),s]=S1[p,1+(s-1)*dim+(q-1)]
                            V[1+(p-1)*dim+(q-1),s]=V1[p,1+(s-1)*dim+(q-1)]
                        end
                    end
                end
            H=pinv(S)*V #pseudo-inverse calculation
            evm=eigen(H).values #eigenvalue calculation
            Vj=V*eigen(H).vectors #calculating of new set of eigenvectors
            #normalizing the results
            for h = 1:mult
                Vjnorm[:,h]=normalize(Vj[:,h])
            end
            sol00=zeros(ComplexF64,n,mult*dim) #creating new initial solution array
            for p=1:n
                for s=1:mult
                    for q=1:dim
                             sol00[p,1+(s-1)*dim+(q-1)]=Vjnorm[1+dim*(p-1)+(q-1),s]
                    end
                end
            end
            eval=hcat(eval,evm)     #filling up the eigenvalue array
            #checking convergence
            if g>2 && checkconv(eval[:,end],eval[:,end-1],0.05)
                 break
            end
            #if check(evm) print(g-1); break end
        end
        #print(size(eval)[2]-1)
        return(norm(normmax(eval[:,end]))-1)
        #return(eval[:,2:end]) #returning the global eigenvalue array
end

###################### SDM ########################
function SDM(kappa1, delta1, epsilon1, b1, tau1, k1)
        dt=tau1/m #definition of time-step
        m2=trunc(Int,m/2)
        dim2=trunc(Int,dim/2)
        dimg=dim*(m2+1)
        v=[kappa1, delta1, epsilon1, 2*pi/(k1*tau1)]
        P=zeros(Float64,dim,dim*m) #construction of Pi matrices
        for j=1:m
            P[:,1+dim*(j-1):dim*j]=expM(A((j-1)*dt,v)*dt)
        end
        R=zeros(Float64,dim,dim*m) #construction of Ri matrices
        for j=1:m
            R[:,1+dim*(j-1):dim*j]=0.5*((expM(A((j-1)*dt,v)*dt)-Matrix{Float64}(I, dim, dim))*inv(A((j-1)*dt,v))*B((j-1)*dt,b1))
        end
        Zw=zeros(Float64,dimg,dimg*m) #construction of Zi matrices
        for j=1:m
            P1=P[:,1+(j-1)*dim:j*dim]
            R1=R[:,1+(j-1)*dim:j*dim]
            Zw[:,1+(j-1)*dim*(m2+1):j*dim*(m2+1)]=vcat(hcat(P1,zeros(dim,dim*(m2-1)),R1[:,dim2+1:dim],R1[:,dim2+1:dim]),hcat(zeros(m2*dim,dim2),Matrix{Float64}(I, m2*dim,m2*dim),zeros(m2*dim,dim2)))
        end
        #construction of final Z matrix
        Z=Zw[:,1+dimg*(m-2):dimg*(m-1)]*Zw[:,1+dimg*(m-1):dimg*m]
        for j=1:(m-2)
            Z=Z*Zw[:,1+dimg*(m-2-j):dimg*(m-1-j)]
        end
        #evaluating stability based on largest eigenvalue
        return(norm(normmax(eigvals(Z)))-1)
    end

####### Multi-Dimensional Bisection Method #########
function foo(x,y)
    return(ISIM(kappa, x, epsilon, y , tau, k))
end

ax1=Axis(-0.95:0.1:0.95,"delta") # initial grid in x direction
ax2=Axis(-1.5:0.5:1.5,"b") # initial grid in y direction

ax1=Axis(0.5:0.25:2.0,"tau") # initial grid in x direction
ax2=Axis(-1.5:0.5:1.5,"b") # initial grid in y direction

mymdbm=MDBM_Problem(foo,[ax1,ax2])
iteration=4 #number of refinements (resolution doubling)
@time solve!(mymdbm,iteration)

x_eval,y_eval=getevaluatedpoints(mymdbm)
x_sol,y_sol=getinterpolatedsolution(mymdbm)

fig = figure(1);clf()
PyPlot.scatter(x_eval,y_eval,s=5)
PyPlot.scatter(x_sol,y_sol,s=5)
