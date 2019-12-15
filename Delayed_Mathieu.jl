using Printf
using Statistics
using LinearAlgebra
using MDBM
using Interpolations
using Plots
pyplot()
using PyPlot
pygui(true);
using DelimitedFiles
using Random
using DifferentialEquations

rng = MersenneTwister(1234)

#creating iteration array
function it(A)
    inttyp=BSpline(Quadratic(Line(OnGrid())))
    #inttyp=BSpline(Constant())
    matrdim=size(A,2)-1
    step=abs(real(A[end,1]-A[1,1]))/(size(A,1)-1)
    scaleitp=real(A[1,1]):step:real(A[end,1])
    ARe=real(A)
    AIm=imag(A)
    imRe=[scale(interpolate(ARe[:,2],inttyp),scaleitp)]
    for j = 3:matrdim+1
        imRetemp=[scale(interpolate(ARe[:,j],inttyp),scaleitp)]
        imRe=vcat(imRe,imRetemp)
    end
    imIm=[scale(interpolate(AIm[:,2],inttyp),scaleitp)]
    for j = 3:matrdim+1
        imImtemp=[scale(interpolate(AIm[:,j],inttyp),scaleitp)]
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
    #if abs((amax-amin)-(bmax-bmin))/(amax-amin)<err
    #   var=true
    #end
    if abs(amax-bmax)/amax<err
       var=true
    end
    return(var)
end

function err(a,b)
    amax=norm(normmax(a))
    bmax=norm(normmax(b))

    err=abs(amax-bmax)/amax
    return(err)
end

#matrix exponential
function expM(M)
    p=12
    n=size(M)[1]
    Mval=zeros(Float64,n,n)
    for j = 0:p
        Mval=Mval+(M^j)/(factorial(j))
    end
    return(Mval)
end

function infnancheck(a)
    var=false
    if any(isnan,a)==true || any(isinf,a)==true
    var=true
    end
return(var)
end

function evplot(a,k)
    U=a[:,k]
    eigRe1=[real(U[n,1]) for n in 1:size(U,1)]
    eigIm1=[imag(U[n,1]) for n in 1:size(U,1)]
    plot1=Plots.plot(eigRe1,eigIm1,seriestype=:scatter)
    Re_circle=[cos(n) for n in 0:0.05:2*pi]
    Im_circle=[sin(n) for n in 0:0.05:2*pi]
    return(Plots.plot!(plot1,Re_circle,Im_circle))
end


function convtime(nmin,nmax,ntest,rep,method1,i1,i2,mult1)
    global n=1500; global i=10; global method="RK4"; global mult=16;
    valref=normmax(ISIM(kappa,delta,epsilon,b,tau,k)[:,end])
    valref=real(valref)+im*abs(imag(valref))
    global method=method1; global i=i1; global mult=mult1;
    imatrix=collect(i1:i2)
    ni=size(imatrix,1)
    nstep=trunc(Int,(nmax-nmin)/ntest)
    nmatrix=collect(nmin:nstep:nmax)
    nst=size(nmatrix,1)
    errortemp=Array{String}(undef,(nst*ni,2+2*rep))
    for in=1:nst
        global n=nmatrix[in]
        for ii=1:ni
            for j=1:2:2*rep-1
                global i=imatrix[ii]
                errortemp[1+(in-1)*ni+(ii-1),1]=string(nmatrix[in])
                errortemp[1+(in-1)*ni+(ii-1),2]=string(imatrix[ii])
                valtemp=@timed ISIM(kappa,delta,epsilon,b,tau,k)
                errortemp[1+(in-1)*ni+(ii-1),j+2]=string(normmax((valtemp[1])[:,end]))
                errortemp[1+(in-1)*ni+(ii-1),j+3]=string(valtemp[2])
            end
        end
    end
    #println(errortemp)
    errormatrix=zeros(Float64,nst*ni,4)
    for in=1:nst
        for ii=1:ni
            errormatrix[1+(in-1)*ni+(ii-1),1]=parse.(Float64,errortemp[1+(in-1)*ni+(ii-1),1])
            errormatrix[1+(in-1)*ni+(ii-1),2]=parse.(Float64,errortemp[1+(in-1)*ni+(ii-1),2])
                if any(x->x=="?", errortemp[ii,:])
                errormatrix[1+(in-1)*ni+(ii-1),2]=0
                else
                vecttemp=parse.(ComplexF64,errortemp[1+(in-1)*ni+(ii-1),3:2:end])
                vecttempabs=real.(vecttemp)+im*abs.(imag.(vecttemp))

                errormatrix[1+(in-1)*ni+(ii-1),3]=mean(abs.(vecttempabs-valref*ones(ComplexF64,rep,1)))
                errormatrix[1+(in-1)*ni+(ii-1),4]=mean(parse.(Float64,errortemp[1+(in-1)*ni+(ii-1),4:2:end]))
                end
        end
    end
    return(errormatrix)
end

A_Jul_3_4=convtime(10,1000,30,4,"Julia",3,3,4)
A_Jul_3_8=convtime(10,1000,30,20,"Julia",3,3,8)
A_Jul_3_16=convtime(10,1000,30,4,"Julia",3,3,16)

A_RK4q_3_4=convtime(10,1000,30,20,"RK4_int",3,3,4)
A_RK4q_3_8=convtime(10,1000,30,20,"RK4_int",3,3,8)
A_RK4q_3_16=convtime(10,1000,30,20,"RK4_int",3,3,16)

A_RK4_3_4=convtime(10,1000,30,20,"RK4",3,3,4)
A_RK4_3_8=convtime(10,1000,30,20,"RK4",3,3,8)
A_RK4_3_16=convtime(10,1000,30,20,"RK4",3,3,16)

A_RK3_3_4=convtime(10,1000,30,20,"RK3",3,3,4)
A_RK3_3_8=convtime(10,1000,30,20,"RK3",3,3,8)
A_RK3_3_16=convtime(10,1000,30,20,"RK3",3,3,16)

A_RK2_3_4=convtime(10,1000,30,20,"RK2",3,3,4)
A_RK2_3_8=convtime(10,1000,30,20,"RK2",3,3,8)
A_RK2_3_16=convtime(10,1000,30,20,"RK2",3,3,16)

A_EE_3_4=convtime(10,1000,30,20,"EE",3,3,4)
A_EE_3_8=convtime(10,1000,30,20,"EE",3,3,8)
A_EE_3_16=convtime(10,1000,30,20,"EE",3,3,16)

Plots.plot(A_Jul_3_4[:,1],hcat(A_Jul_3_4[:,3],A_Jul_3_8[:,3],A_Jul_3_16[:,3]),xscale=:log10, yscale=:log10,title="Eigval error",label=["mult=4" "mult=8" "mult=16"])
Plots.plot(A_Jul_3_4[:,1],hcat(A_Jul_3_4[:,4],A_Jul_3_8[:,4],A_Jul_3_16[:,4]),xscale=:log10, yscale=:log10,title="Time",label=["mult=4" "mult=8" "mult=16"])
Plots.plot(hcat(A_Jul_3_4[:,4],A_Jul_3_8[:,4],A_Jul_3_16[:,4]),hcat(A_Jul_3_4[:,3],A_Jul_3_8[:,3],A_Jul_3_16[:,3]),xscale=:log10, yscale=:log10,title="Time vs precision",label=["mult=4" "mult=8" "mult=16"])

Plots.plot(A_RK4q_3_4[:,1],hcat(A_RK4q_3_4[:,3],A_RK4q_3_8[:,3],A_RK4q_3_16[:,3]),xscale=:log10, yscale=:log10,title="Eigval error",label=["mult=4" "mult=8" "mult=16"])
Plots.plot(A_RK4q_3_4[:,1],hcat(A_RK4q_3_4[:,4],A_RK4q_3_8[:,4],A_RK4q_3_16[:,4]),xscale=:log10, yscale=:log10,title="Time",label=["mult=4" "mult=8" "mult=16"])
Plots.plot(hcat(A_RK4q_3_4[:,4],A_RK4q_3_8[:,4],A_RK4q_3_16[:,4]),hcat(A_RK4q_3_4[:,3],A_RK4q_3_8[:,3],A_RK4q_3_16[:,3]),xscale=:log10, yscale=:log10,title="Time vs precision",label=["mult=4" "mult=8" "mult=16"])


Plots.plot(A_RK4_3_4[:,1],hcat(A_RK4_3_4[:,3],A_RK3_3_4[:,3],A_RK2_3_4[:,3],A_EE_3_4[:,3]),xscale=:log10, yscale=:log10,title="Eigval error mult=4 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])

Plots.plot(A_RK4_3_8[:,1],hcat(A_RK4_3_8[:,3],A_RK3_3_8[:,3],A_RK2_3_8[:,3],A_EE_3_8[:,3]),xscale=:log10, yscale=:log10,title="Eigval error mult=8 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])

Plots.plot(A_RK4_3_16[:,1],hcat(A_RK4_3_16[:,3],A_RK3_3_16[:,3],A_RK2_3_16[:,3],A_EE_3_16[:,3]),xscale=:log10, yscale=:log10,title="Eigval error mult=16 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])


Plots.plot(A_RK4_3_4[:,1],hcat(A_RK4_3_4[:,4],A_RK3_3_4[:,4],A_RK2_3_4[:,4],A_EE_3_4[:,4]),xscale=:log10, yscale=:log10,title="Time mult=4 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])

Plots.plot(A_RK4_3_8[:,1],hcat(A_RK4_3_8[:,4],A_RK3_3_8[:,4],A_RK2_3_8[:,4],A_EE_3_8[:,4]),xscale=:log10, yscale=:log10,title="Time mult=8 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])

Plots.plot(A_RK4_3_16[:,1],hcat(A_RK4_3_16[:,4],A_RK3_3_16[:,4],A_RK2_3_16[:,4],A_EE_3_16[:,4]),xscale=:log10, yscale=:log10,title="Time mult=16 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])



Plots.plot(hcat(A_RK4_3_4[:,4],A_RK3_3_4[:,4],A_RK2_3_4[:,4],A_EE_3_4[:,4]),hcat(A_RK4_3_4[:,3],A_RK3_3_4[:,3],A_RK2_3_4[:,3],A_EE_3_4[:,3]),xscale=:log10, yscale=:log10,title="Time vs precision mult=4 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])

Plots.plot(hcat(A_RK4_3_8[:,4],A_RK3_3_8[:,4],A_RK2_3_8[:,4],A_EE_3_8[:,4]),hcat(A_RK4_3_8[:,3],A_RK3_3_8[:,3],A_RK2_3_8[:,3],A_EE_3_8[:,3]),xscale=:log10, yscale=:log10,title="Time vs precision mult=8 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])

Plots.plot(hcat(A_RK4_3_16[:,4],A_RK3_3_16[:,4],A_RK2_3_16[:,4],A_EE_3_16[:,4]),hcat(A_RK4_3_16[:,3],A_RK3_3_16[:,3],A_RK2_3_16[:,3],A_EE_3_16[:,3]),xscale=:log10, yscale=:log10,title="Time vs precision mult=16 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])



###################### ISIM #######################
mult=10 #Ns=dim*mult
n=40 #timestep number
#EE: Explicit Euler
#RK4: 4th order Runge-Kutta
#RK4_int: 4th order Runge-Kutta with interpolation function
#IE: Implicit Euler
#RK2: 2nd order Runge-Kutta
method="Julia" #choosing numerical simulation method
i=4 #number of iteration

###################### SDM ########################
m=40 #discretization number


#parameters of the delayed Mathieu equation
dim=2 #DoF of the system (after Cauchy transcription)
const tau=2*pi
const kappa=0.1
const delta=0.5
const epsilon=1.0
const b=0.2
const k=2
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

function B(t,v)
    Bmatr=zeros(dim*dim)
    Bmatr=[0 v[5]; 0 0]
    return(Bmatr)
end

function Ai(v)
    return(t->A(t,v))
end

function Bi(v)
    return(t->B(t,v))
end

function bc_model(du,u,h,p,t)
    kappa,delta,epsilon,omega,b,tau,AA,BB = p
    hist = h(p, t-tau)
    dutemp=zeros(ComplexF64,dim)
    for j=1:dim
        for jj=1:dim
        dutemp[j] = dutemp[j]+AA(t)[j,jj]*u[jj]+BB(t)[j,jj]*hist[jj]
        end
    end
    for j=1:dim
    du[j] = dutemp[j]
    end
end

function fsol(v,tau,k,IF,dt1)
    alg = MethodOfSteps(Tsit5())
    lags = [tau]
    p = (v[1],v[2],v[3],v[4],v[5],tau,Ai(v),Bi(v))
    interp=it(IF)
    h(p,t)=sub(interp,t)
    tspan = (0,k*tau+0.5*dt1)
    u0 = IF[n,2:2+(dim-1)]
    prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)
    return(solve(prob,alg,abstol=1e-13,reltol=1e-13,maxiters=Int(1e9)))
end

method="Julia"
aaa=ISIM(kappa,delta,epsilon,b,tau,k)
method="RK4"
bbb=ISIM(kappa,delta,epsilon,b,tau,k)
evplot(aaa,6)
evplot(bbb,6)

###################### ISIM #######################
function ISIM(kappa1,delta1,epsilon1,b1,tau1,k1)

        dt=tau1/(n-1)
        v=[kappa1, delta1, epsilon1, 2*pi/(k1*tau1),b1]

        eval=zeros(ComplexF64,mult,1)
        S=zeros(ComplexF64,n*dim,mult)
        V=zeros(ComplexF64,n*dim,mult)
        Vjnorm=zeros(ComplexF64,n*dim,mult)
        #tvect0=collect(-tau1:dt:tau1*k1+0.0000001*dt)
        tvect0=collect(-tau1:dt:tau1*k1) #discretized time vector
        sol=zeros(ComplexF64,k1*(n-1),mult*dim) #empty solution vector
        sol00=randn!(rng, zeros(n,mult*dim))#S0 inital random array
           #solreturn=zeros(ComplexF64,1,1+mult*dim)
        for g = 1:i #iteration
        sol0=hcat(tvect0,vcat(sol00,sol)) #constructing solution array
            #time integration
            for m = 2:dim:(dim*mult-(dim-2)) #choosing "one physical system"
            sol0m=hcat(sol0[:,1],sol0[:,m:m+dim-1])
            #choosing numerical method

            if method == "IE"
                for j = 0:k1*(n-1)-1
                    X=sol0m[j+n,2:2+(dim-1)]
                    Yt=sol0m[j+n,2:2+(dim-1)]
                    Ytau=sol0m[j+1,2:2+(dim-1)]
                    for l=1:3
                        X=X-inv(dt*A(j*dt,v)-Matrix{Float64}(I,dim,dim))*(Yt+dt*(A((j+1+n)*dt,v)*X+B((j+1+n)*dt,v)*Ytau)-X)
                    end
                    Y=X
                   sol0m[n+j+1,2:2+(dim-1)]=transpose(Y) #filling solution array of the physical system
                end
            elseif method == "EE"
                    for j = 0:k1*(n-1)-1
                        Y=sol0m[n+j,2:2+(dim-1)]+dt*(A(j*dt,v)*sol0m[n+j,2:2+(dim-1)]+B(j*dt,v)*sol0m[1+j,2:2+(dim-1)])
                        sol0m[n+j+1,2:2+(dim-1)]=transpose(Y)
                    end
            elseif method == "RK2"
                for j = 0:k1*(n-1)-1

                    Ydtau=0.5*(sol0m[j+1,2:2+(dim-1)]+sol0m[j+2,2:2+(dim-1)])

                    Y1=sol0m[n+j,2:2+(dim-1)]
                    Y2=Y1+(dt/2)*(A((j+0.5)*dt,v)*Y1+B((j+0.5)*dt,v)*Ydtau)

                    Y=Y1+dt*(A((j+0.5)*dt,v)*Y2+B((j+0.5)*dt,v)*Ydtau)

                    sol0m[n+j+1,2:2+(dim-1)]=transpose(Y)
                end
            elseif method == "RK3"
                for j = 0:k1*(n-1)-1

                    Ydtau1=sol0m[j+1,2:2+(dim-1)]
                    Ydtau12=0.5*(sol0m[j+1,2:2+(dim-1)]+sol0m[j+2,2:2+(dim-1)])
                    Ydtau2=sol0m[j+2,2:2+(dim-1)]

                    Y1=sol0m[n+j,2:2+(dim-1)]

                    F1=A(j*dt,v)*Y1+B(j*dt,v)*Ydtau1
                    F2=A((j+0.5)*dt,v)*(Y1+dt*F1/2)+B((j+0.5)*dt,v)*Ydtau12
                    F3=A((j+1)*dt,v)*(Y1+dt*(-F1+2*F2))+B((j+0.5)*dt,v)*Ydtau2

                    Y=Y1+dt/6*(F1+4*F2+F3)

                    sol0m[n+j+1,2:2+(dim-1)]=transpose(Y)
                end
            elseif method == "RK4"
                for j = 0:k1*(n-1)-1
                    Ydtau=0.5*(sol0m[j+1,2:2+(dim-1)]+sol0m[j+2,2:2+(dim-1)])

                    Y1=sol0m[n+j,2:2+(dim-1)]
                    Y2=Y1+(dt/2)*(A(j*dt,v)*Y1+B(j*dt,v)*sol0m[j+1,2:2+(dim-1)])
                    Y3=Y1+(dt/2)*(A((j+0.5)*dt,v)*Y2+B((j+0.5)*dt,v)*Ydtau)
                    Y4=Y1+dt*(A((j+0.5)*dt,v)*Y3+B((j+0.5)*dt,v)*Ydtau)

                    Y=Y1+(dt/6)*((A(j*dt,v)*Y1+B(j*dt,v)*sol0m[j+1,2:2+(dim-1)])+2*(A((j+0.5)*dt,v)*Y2+B((j+0.5)*dt,v)*Ydtau)+2*(A((j+0.5)*dt,v)*Y3+B((j+0.5)*dt,v)*sol0m[j+1,2:2+(dim-1)])+(A((j+1)*dt,v)*Y4+B(j*dt,v)*sol0m[j+2,2:2+(dim-1)]))

                    sol0m[n+j+1,2:2+(dim-1)]=transpose(Y)
                end
            elseif method == "RK4_int"
                for jk=1:k1
                    sol0mittau=it(sol0m[1+(jk-1)*(n-1):n+(jk-1)*(n-1),:])
                    for j = 0:(n-1)-1
                        jj=j+(jk-1)*(n-1)

                        Ytau0=sub(sol0mittau,-tau1+jj*dt)
                        Ytau05=sub(sol0mittau,-tau1+(jj+0.5)*dt)
                        Ytau1=sub(sol0mittau,-tau1+(jj+1)*dt)

                        Y1=sol0m[jj+n,2:dim+1]
                        Y2=Y1+(dt/2)*(A(jj*dt,v)*Y1+B(jj*dt,v)*Ytau0)
                        Y3=Y1+(dt/2)*(A((jj+0.5)*dt,v)*Y2+B((jj+0.5)*dt,v)*Ytau05)
                        Y4=Y1+dt*(A((jj+0.5)*dt,v)*Y3+B((jj+0.5)*dt,v)*Ytau05)

                        Y=Y1+(dt/6)*((A(jj*dt,v)*Y1+B(jj*dt,v)*Ytau0)+2*(A((jj+0.5)*dt,v)*Y2+B((jj+0.5)*dt,v)*Ytau05)+2*(A((jj+0.5)*dt,v)*Y3+B((jj+0.5)*dt,v)*Ytau05)+(A((jj+1)*dt,v)*Y4+B((jj+1)*dt,v)*Ytau1))

                        sol0m[n+jj+1,2:2+(dim-1)]=transpose(Y)
                    end
                end
            elseif method == "Julia"
                    solJulia=fsol(v,tau1,k1,sol0m[1:n,:],dt)
                    for j = 1:k1*(n-1)
                        sol0m[n+j,2:2+(dim-1)]=transpose(solJulia(dt*j)[1:dim,1:end])
                    end
                end

                sol0[:,m:m+dim-1]=sol0m[:,2:2+(dim-1)] #filling solution array
            end
            #solreturn=vcat(solreturn,sol0)
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

            #H=pinv(S,rtol = sqrt(eps(real(float(one(eltype(S)))))))*V

                    if infnancheck(H)==true
                        return("?")
                    end
                    eigH=eigen(H)
                    evm=eigH.values #eigenvalue calculation
                    Vjnorm=V*eigH.vectors #calculating of new set of eigenvectors
                    #normalizing the results
                    #for h = 1:mult
                    #    Vjnorm[:,h]=normalize(Vj[:,h])
                    #end
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
                    #if g>2 && checkconv(eval[:,end],eval[:,end-1],0.01)
                    #     break
                    #end
                    #if check(evm) print(g-1); break end
                end
                print(size(eval)[2]-1)
                return(eval)
                #return(eval[:,2:end]) #returning the global eigenvalue array
                #return(solreturn)

                #return(err(eval[:,end-1],eval[:,end]))
end



###################### SDM ########################
function SDM(kappa1, delta1, epsilon1, b1, tau1, k1)
        dt=tau1/m #definition of time-step
        m2=trunc(Int,m/2)
        dim2=trunc(Int,dim/2)
        dimg=dim*(m2+1)
        v=[kappa1, delta1, epsilon1, 2*pi/(k1*tau1)]
        P=zeros(Float64,dim,dim*m*k1) #construction of Pi matrices
        for j=1:m*k1
            P[:,1+dim*(j-1):dim*j]=expM(A((j-1)*dt,v)*dt)
        end
        R=zeros(Float64,dim,dim*m*k1) #construction of Ri matrices
        for j=1:m*k1
            R[:,1+dim*(j-1):dim*j]=0.5*((expM(A((j-1)*dt,v)*dt)-Matrix{Float64}(I, dim, dim))*inv(A((j-1)*dt,v))*B((j-1)*dt,b1))
        end
        Zw=zeros(Float64,dimg,dimg*m*k1) #construction of Zi matrices
        for j=1:m*k1
            P1=P[:,1+(j-1)*dim:j*dim]
            R1=R[:,1+(j-1)*dim:j*dim]
            Zw[:,1+(j-1)*dim*(m2+1):j*dim*(m2+1)]=vcat(hcat(P1,zeros(dim,dim*(m2-1)),R1[:,dim2+1:dim],R1[:,dim2+1:dim]),hcat(zeros(m2*dim,dim2),Matrix{Float64}(I, m2*dim,m2*dim),zeros(m2*dim,dim2)))
        end
        #construction of final Z matrix
        Z=Zw[:,1+dimg*(m-2)+m*dimg*(k1-1):dimg*(m-1)+m*dimg*(k1-1)]*Zw[:,1+dimg*(m-1)+m*dimg*(k1-1):dimg*m+m*dimg*(k1-1)]
        for j=1:(m-2)+(k1-1)*m
            Z=Z*Zw[:,1+dimg*(m-2-j)+m*dimg*(k1-1):dimg*(m-1-j)+m*dimg*(k1-1)]
        end
        #evaluating stability based on largest eigenvalue
        return(eigvals(Z))
        #return(P)
    end

####### Multi-Dimensional Bisection Method #########
function foo(x,y)
    return(norm(normmax(ISIM(kappa,x,epsilon,y,tau,k)[:,end]))-1)
end

ax1=Axis(-1.0:0.5:5,"delta") # initial grid in x direction
ax2=Axis(-1.5:0.5:1.5,"b") # initial grid in y direction

ax1=Axis(0.1:0.25:5.0,"tau") # initial grid in x direction
ax2=Axis(0.05:0.2:1.5,"b") # initial grid in y direction

mymdbm=MDBM_Problem(foo,[ax1,ax2])
iteration=3 #number of refinements (resolution doubling)
@time MDBM.solve!(mymdbm,iteration)

x_eval,y_eval=getevaluatedpoints(mymdbm)
x_sol,y_sol=getinterpolatedsolution(mymdbm)

fig = figure(1);clf()
PyPlot.scatter(x_eval,y_eval,s=5)
PyPlot.scatter(x_sol,y_sol,s=5)
