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

rng = MersenneTwister(1234)

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


function convtime(nmin,nmax,ntest,rep,method1,i1,i2,mult1)
    global n=1500; global i=10; global method="RK4"; global mult=16;
    valref=normmax(ISIM(omega,w,z,ae,d,strat)[:,end])
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
                valtemp=@timed ISIM(omega,w,z,ae,d,strat)
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

Plots.plot(A_RK4_3_4[:,1],hcat(A_RK4_3_4[:,3],A_RK3_3_4[:,3],A_RK2_3_4[:,3],A_EE_3_4[:,3]),xscale=:log10, yscale=:log10,title="Eigval error mult=4 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])

Plots.plot(A_RK4_3_8[:,1],hcat(A_RK4_3_8[:,3],A_RK3_3_8[:,3],A_RK2_3_8[:,3],A_EE_3_8[:,3]),xscale=:log10, yscale=:log10,title="Eigval error mult=8 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])

Plots.plot(A_RK4_3_16[:,1],hcat(A_RK4_3_16[:,3],A_RK3_3_16[:,3],A_RK2_3_16[:,3],A_EE_3_16[:,3]),xscale=:log10, yscale=:log10,title="Eigval error mult=16 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])


Plots.plot(A_RK4_3_4[:,1],hcat(A_RK4_3_4[:,4],A_RK3_3_4[:,4],A_RK2_3_4[:,4],A_EE_3_4[:,4]),xscale=:log10, yscale=:log10,title="Time mult=4 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])

Plots.plot(A_RK4_3_8[:,1],hcat(A_RK4_3_8[:,4],A_RK3_3_8[:,4],A_RK2_3_8[:,4],A_EE_3_8[:,4]),xscale=:log10, yscale=:log10,title="Time mult=8 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])

Plots.plot(A_RK4_3_16[:,1],hcat(A_RK4_3_16[:,4],A_RK3_3_16[:,4],A_RK2_3_16[:,4],A_EE_3_16[:,4]),xscale=:log10, yscale=:log10,title="Time mult=16 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])



Plots.plot(hcat(A_RK4_3_4[:,4],A_RK3_3_4[:,4],A_RK2_3_4[:,4],A_EE_3_4[:,4]),hcat(A_RK4_3_4[:,3],A_RK3_3_4[:,3],A_RK2_3_4[:,3],A_EE_3_4[:,3]),xscale=:log10, yscale=:log10,title="Time vs precision mult=4 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])

Plots.plot(hcat(A_RK4_3_8[:,4],A_RK3_3_8[:,4],A_RK2_3_8[:,4],A_EE_3_8[:,4]),hcat(A_RK4_3_8[:,3],A_RK3_3_8[:,3],A_RK2_3_8[:,3],A_EE_3_8[:,3]),xscale=:log10, yscale=:log10,title="Time vs precision mult=8 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])

Plots.plot(hcat(A_RK4_3_16[:,4],A_RK3_3_16[:,4],A_RK2_3_16[:,4],A_EE_3_16[:,4]),hcat(A_RK4_3_16[:,3],A_RK3_3_16[:,3],A_RK2_3_16[:,3],A_EE_3_16[:,3]),xscale=:log10, yscale=:log10,title="Time vs precision mult=16 10000/0.02",label=["RK4" "RK3" "RK2" "EE"])

function evplot(a,k)
    U=a[:,k]
    eigRe1=[real(U[n,1]) for n in 1:size(U,1)]
    eigIm1=[imag(U[n,1]) for n in 1:size(U,1)]
    plot1=Plots.plot(eigRe1,eigIm1,seriestype=:scatter)
    Re_circle=[cos(n) for n in 0:0.05:2*pi]
    Im_circle=[sin(n) for n in 0:0.05:2*pi]
    return(Plots.plot!(plot1,Re_circle,Im_circle))
end

# A_RK3_3_8=convtime(10,300,20,5,"RK3",3,8)
# A_RK3_3_16=convtime(10,300,20,5,"RK3",3,16)
# A_RK3_4_8=convtime(10,300,20,5,"RK3",4,8)
# A_RK3_4_16=convtime(10,300,20,5,"RK3",4,16)
# A_RK3_5_8=convtime(10,300,20,5,"RK3",5,8)
# A_RK3_5_16=convtime(10,300,20,5,"RK3",5,16)
#
# A_RK2_3_8=convtime(10,300,20,5,"RK2",3,8)
# A_RK2_3_16=convtime(10,300,20,5,"RK2",3,16)
# A_RK2_4_8=convtime(10,300,20,5,"RK2",4,8)
# A_RK2_4_16=convtime(10,300,20,5,"RK2",4,16)
# A_RK2_5_8=convtime(10,300,20,5,"RK2",5,8)
# A_RK2_5_16=convtime(10,300,20,5,"RK2",5,16)

open("A_RK2_5_16.txt", "w") do io
            writedlm(io, A_RK2_5_16)
       end


###################### ISIM #######################
mult=12 #Ns=dim*mult
n=40 #timestep number
#EE: Explicit Euler
#RK4: 4th order Runge-Kutta
#RK4_int: 4th order Runge-Kutta with interpolation function
#IE: Implicit Euler
#RK2: 2nd order Runge-Kutta
method="RK4" #choosing numerical simulation method
i=4 #number of iteration

###################### SDM ########################
m=50 #discretization number

#parameters of the milling model
dim=4 #DoF of the system (after Cauchy transcription)

#machine and tool fixation parameters
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

omega=10000
w=0.002
#tool coefficient matrices
Mm=[mx 0; 0 my];
Km=[kx 0; 0 ky];
Sm=[sx 0; 0 sy];
strat="Down"


function g(phi,phi_in1,phi_out1)
    g_val=0
    if phi < phi_in1
        g_val=0
    elseif (phi >= phi_in1) && (phi <= phi_out1)
        g_val=1
    #elseif phi > phi_out1
    else
        g_val=0
    end
    return(g_val)
end



function W(t,v1)
    omega1=v1[1]
    w1=v1[2]
    z1=Int(v1[3])
    phi_in1=v1[4]
    phi_out1=v1[5]

    phi=zeros(z1)
    for j=0:(z1-1)
    phi[j+1]=mod(t*omega1+j*2*pi/z1,2*pi)
    end
    W_f=zeros(Float64,2,2)
    for j = 0:(z1-1)
        W_f=W_f+g(phi[j+1],phi_in1,phi_out1)*[-sin(2*phi[j+1])-kr+kr*cos(2*phi[j+1]) -1-cos(2*phi[j+1])-kr*sin(2*phi[j+1]); 1-cos(2*phi[j+1])-kr*sin(2*phi[j+1]) sin(2*phi[j+1])-kr-kr*cos(2*phi[j+1])]
    end
    return ((w1*kt/2)*W_f)
end

function A(t,v1)
    vcat(hcat(-inv(Mm)*Km,-inv(Mm)*(Sm-W(t,v1))),hcat(Matrix{Float64}(I, 2, 2),zeros(2,2)))
end

function B(t,v1)
    vcat(hcat(zeros(2,2),-inv(Mm)*W(t,v1)),hcat(zeros(2,2),zeros(2,2)))
end

################### W(t) matric plot ###############
# dim0=2
# if strat=="Up"
#     aepu=ae/d; phi_in=0; phi_out=acos(1-2*aepu);
# elseif strat=="Down"
#     aepd=2-ae/d; phi_in=acos(1-2(aepd-1)); phi_out=pi;
# end
# tau=2*pi/(omega*z)
#
# v=[omega,w,z,phi_in,phi_out]
#
# t0=0
# t1=2*tau
# dtp=(t1-t0)/200
# tvec=collect(t0:dtp:t1)
# stvec=size(tvec,1)
#
# plotarr=zeros(Float64,dim0*dim0,stvec)
# for ind1=1:dim0
#     for ind2=1:dim0
#         for ind3=1:stvec
#             plotarr[1+(ind1-1)*dim0+(ind2-1),ind3]=W(tvec[ind3],v)[ind1,ind2]
#         end
#     end
# end
# fig = figure(1);clf()
# PyPlot.plot(tvec,transpose(plotarr))
########################################################

###################### ISIM #######################
function ISIM(omega1,w1,z1,ae1,d1,strat1)
        k1=1
        if strat1=="Up"
            aepu=ae1/d1; phi_in=0; phi_out=acos(1-2*aepu);
        elseif strat=="Down"
            aepd=2-ae1/d1; phi_in=acos(1-2(aepd-1)); phi_out=pi;
        end
        tau1=2*pi/(omega1*z1)
        dt=tau1/(n-1)

        v=[omega1,w1,z1,phi_in,phi_out]

        eval=zeros(ComplexF64,mult,1)
        S=zeros(ComplexF64,n*dim,mult)
        V=zeros(ComplexF64,n*dim,mult)
        Vjnorm=zeros(ComplexF64,n*dim,mult)
        tvect0=collect(-tau1:dt:tau1*k1+0.00001*dt) #discretized time vector
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
                for j = 0:k1*(n-1)-1
                    sol0mittau=it(sol0m[j+1:j+3,:])

                    Ytau0=sub(sol0mittau,-tau1+j*dt)
                    Ytau05=sub(sol0mittau,-tau1+(j+0.5)*dt)
                    Ytau1=sub(sol0mittau,-tau1+(j+1)*dt)

                    Y1=sol0m[j+n,2:dim+1]
                    Y2=Y1+(dt/2)*(A(j*dt,v)*Y1+B(j*dt,v)*Ytau0)
                    Y3=Y1+(dt/2)*(A((j+0.5)*dt,v)*Y2+B((j+0.5)*dt,v)*Ytau05)
                    Y4=Y1+dt*(A((j+0.5)*dt,v)*Y3+B((j+0.5)*dt,v)*Ytau05)

                    Y=Y1+(dt/6)*((A(j*dt,v)*Y1+B(j*dt,v)*Ytau0)+2*(A((j+0.5)*dt,v)*Y2+B((j+0.5)*dt,v)*Ytau05)+2*(A((j+0.5)*dt,v)*Y3+B((j+0.5)*dt,v)*Ytau05)+(A((j+1)*dt,v)*Y4+B((j+1)*dt,v)*Ytau1))

                    sol0m[n+j+1,2:2+(dim-1)]=transpose(Y)
                end
            end

                sol0[:,m:m+dim-1]=sol0m[:,2:2+(dim-1)] #filling solution array
            end
            #solreturn=vcat(solreturn,sol0)
            #return(solreturn)
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
function SDM(omega1,w1,z1,ae1,d1,strat1)
        k1=1
        if strat1=="Up"
            aepu=ae1/d1; phi_in=0; phi_out=acos(1-2*aepu);
        elseif strat1=="Down"
            aepd=2-ae1/d1; phi_in=acos(1-2(aepd-1)); phi_out=pi;
        end
        tau1=2*pi/(omega1*z1)
        dt=tau1/m

        v=[omega1,w1,z1,phi_in,phi_out]

        m2=trunc(Int,m/2)
        dim2=trunc(Int,dim/2)
        dimg=dim*(m2+1)

        P=zeros(Float64,dim,dim*m*k1) #construction of Pi matrices
        for j=1:m*k1
            P[:,1+dim*(j-1):dim*j]=expM(A((j-1)*dt,v)*dt)
        end
        R=zeros(Float64,dim,dim*m*k1) #construction of Ri matrices
        for j=1:m*k1
            R[:,1+dim*(j-1):dim*j]=0.5*((expM(A((j-1)*dt,v)*dt)-Matrix{Float64}(I, dim, dim))*inv(A((j-1)*dt,v))*B((j-1)*dt,v))
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
    return(norm(normmax(ISIM(x,y,z,ae,d,strat)[:,end]))-1)
end

ax1=Axis(2000.0:50.0:3000.0,"omega") # initial grid in x direction
ax2=Axis(0.000:0.0005:0.005,"w") # initial grid in y direction

mymdbm=MDBM_Problem(foo,[ax1,ax2])
iteration=3 #number of refinements (resolution doubling)
@time solve!(mymdbm,iteration)

x_eval,y_eval=getevaluatedpoints(mymdbm)
x_sol,y_sol=getinterpolatedsolution(mymdbm)

fig = figure(1);clf()
PyPlot.scatter(x_eval,y_eval,s=5)
PyPlot.scatter(x_sol,y_sol,s=5)
