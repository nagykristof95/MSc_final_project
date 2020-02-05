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

#ISIM multi-numerikus szimuláció
#konvergencia figyelés és finomítás "valós időben"
#lehetőség beépített numerikus szimulációra
#tekintet nélkül tau és T viszonyára (?)


################### SYSTEM DEFINITION ###################
#function definiton : f(t)=A(t)y(t)+B(t)y(t-tau)
dim=4 #SYSTEM DIMENSION (Post-Cauchy)

#SYSTEM PARAMETERS
mx=0.01986;my=0.02008;kx=1.60312;ky=1.155697;sx=408866;sy=413445; #structural parameters
kt=644*10^6;kr=0.368; #cutting parameters
d=0.008;ae=0.0004;z=1;strat=1.0; #tool/strategy parameters


omega=5000 #spindle speed
w=0.002 #axial immersion
#v=[1 mx  2 my  3 kx  4 ky  5 sx  6 sy  7 kt  8 kr  9 d  10 ae  11 z  12 strat  13 omega  14 w ]
v=[mx my kx ky sx sy kt kr d ae z strat omega w]

#time delay
function tau(v1)
    (2*pi)/(v1[13]*v1[11])
end

#period
function T(v1)
    (2*pi)/v1[13]
end

function Mm(t,v1) #structural matrix Mm
    return([v1[1] 0; 0 v1[2]])
end

function Km(t,v1) #structural matrix Km
    return([v1[3] 0; 0 v1[4]])
end

function Sm(t,v1) #structural matrix Sm
    return([v1[5] 0; 0 v1[6]])
end

function g(phi,phi_in1,phi_out1) #sign function (working teeth)
    g_val=0
    if phi < phi_in1
        g_val=0
    elseif (phi >= phi_in1) && (phi <= phi_out1)
        g_val=1
    else
        g_val=0
    end
    return(g_val)
end

function W(t,v1) #periodic matrix
    if v1[12]==1.0
        aepu=v1[10]/v1[9]; phi_in=0; phi_out=acos(1-2*aepu);
    else
        aepd=2-v1[10]/v1[9]; phi_in=acos(1-2(aepd-1)); phi_out=pi;
    end
    zint=trunc(Int,v1[11])
    phi=zeros(zint)
    for j=0:(zint-1)
    phi[j+1]=mod(t*v1[13]+j*2*pi/v1[11],2*pi)
    end
    W_f=zeros(Float64,2,2)
    for j = 0:(zint-1)
        W_f=W_f+g(phi[j+1],phi_in,phi_out)*[-sin(2*phi[j+1])-v1[8]+v1[8]*cos(2*phi[j+1]) -1-cos(2*phi[j+1])-v1[8]*sin(2*phi[j+1]); 1-cos(2*phi[j+1])-v1[8]*sin(2*phi[j+1]) sin(2*phi[j+1])-v1[8]-v1[8]*cos(2*phi[j+1])]
    end
    return ((v1[14]*v1[7]/2)*W_f)
end

function A(t,v1) #coefficient matrix A
    vcat(hcat(-inv(Mm(t,v1))*Km(t,v1),-inv(Mm(t,v1))*(Sm(t,v1)-W(t,v1))),hcat(Matrix{Float64}(I, 2, 2),zeros(2,2)))
end

function B(t,v1) #coefficient matrix B (delay)
    vcat(hcat(zeros(2,2),-inv(Mm(t,v1))*W(t,v1)),hcat(zeros(2,2),zeros(2,2)))
end

###################################################################################

################### SYSTEM DEFINITION ###################

dim=2 #DoF of the system (after Cauchy transcription)

omega=1; kappa=0.2; delta=0.5; epsilon=1.0; b=0.2; tau0=2*pi;

#v=[1 omega  2 kappa  3 delta  4 epsilon  5 b  6 tau0]
v=[omega kappa delta epsilon b tau0]

#time delay
function tau(v1)
    v1[6]
end

#period
function T(v1)
    (2*pi)/v1[1]
end

#matrix coefficients
function A(t,v1)
    [-v1[2] -(v1[3]+v1[4]*cos(v1[1]*t)); 1 0]
end

function B(t,v1)
    [0 v1[5]; 0 0]
end
#########################################################

################### COMPUTATION PARAMETERS DEFINITION ###################

###################################################################################

function f(t,y,ytau,v1) #right-hand side for TNS
    A(t,v1)*y+B(t,v1)*ytau
end

function Amult(t,v1,mult1)
    Atemp=zeros(ComplexF64,mult1*dim,mult1*dim)
    for i1=0:dim:(mult1-1)*dim
        for j=1:dim
            for jj=1:dim
            Atemp[i1+j,i1+jj] =  A(t,v1)[j,jj]
            end
        end
    end
    return(Atemp)
end

function Bmult(t,v1,mult1)
    Btemp=zeros(ComplexF64,mult1*dim,mult1*dim)
    for i1=0:dim:(mult1-1)*dim
        for j=1:dim
            for jj=1:dim
            Btemp[i1+j,i1+jj] =  B(t,v1)[j,jj]
            end
        end
    end
    return(Btemp)
end

function fmult(t,y,ytau,v1,mult1) #right-hand side for TNS
    Amult(t,v1,mult1)*y+Bmult(t,v1,mult1)*ytau
end



#auxiliar functions definitions for DE.jl
function Ai(v1) #coefficient matrix A DE.jl
    return(t->A(t,v1))
end

function Bi(v1) #coefficient matrix B DE.jl
    return(t->B(t,v1))
end

function bc_model(du,u,h,p,t) #DE.jl problem definiton
    tau,AA,BB,mult1 = p
    hist = h( out1,p, t-tau)
    dutemp=zeros(ComplexF64,mult*dim)
    for i1=0:dim:(mult-1)*dim
        for j=1:dim
            for jj=1:dim
            dutemp[j+i1] = dutemp[j+i1]+AA(t)[j,jj]*u[jj+i1]+BB(t)[j,jj]*out1[jj+i1]
            end
        end
    end
    for j=1:mult*dim
    du[j] = dutemp[j]
    end
end



function fsol(tau,IF,AA,BB,tend,mult1,dt1) #DE.jl solution defintion from 0 to tend
    alg =  MethodOfSteps(BS3())
    #alg =  MethodOfSteps(AutoTsit5(Rosenbrock23(autodiff=false)))
    #alg=MethodOfSteps(RK4())
    lags = [tau]
    p = (tau,AA,BB,mult1)
    interp=it(IF)
    global out1=zeros(ComplexF64,dim*mult1)  #define a cache variable
    h(out1,p,t)=(out1.=sub(interp,t))
    tspan = (0.0,tend)
    u0 = sub(interp,0.0)
    prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)
    return(solve(prob,alg,dtmax=dt1,dtmin=dt1,force_dtmin=true,saveat=collect(0.0:dt1:tend)))
    #return(solve(prob,alg,dt=dt1,dtmax=dt1,dtmin=dt1,force_dtmin=true,saveat=collect(0.0:dt1:tend)))
end


function butcher(t,inttau,y,dt1,(Ba1,Bb1,Bc1),v1,tau1,mult1) #one step by Butcher table (explicit only!)
    s=size(Bb1)[1]
    kvec=zeros(ComplexF64,dim*mult1,s)
    Svec=zeros(ComplexF64,dim*mult1,s)
    for j=1:s
        for jj=1:s
            Svec[:,j]=Svec[:,j]+Ba1[j,jj]*kvec[:,jj]
        end
        kvec[:,j]=dt1*fmult(t+Bc1[j]*dt1,y+Svec[:,j],sub(inttau,t+Bc1[j]*dt1-tau1),v1,mult1)
    end
    yn=y
    for j=1:s
        yn=yn+Bb1[j]*kvec[:,j]
    end
    return(yn)
end


BaRK4=[0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0]
BbRK4=[1/6, 1/3, 1/3, 1/6]
BcRK4=[0, 0.5, 0.5, 1]
BRK4=(BaRK4,BbRK4,BcRK4)

BaRK3=[0 0 0; 0.5 0 0;-1.0 2.0 0;]
BbRK3=[1/6, 2/3, 1/6]
BcRK3=[0, 0.5, 1.0]
BRK3=(BaRK3,BbRK3,BcRK3)

BaRK2=[0 0; 0.5 0]
BbRK2=[0, 1.0]
BcRK2=[0, 0.5]
BRK2=(BaRK2,BbRK2,BcRK2)

BaEE=[0.0]
BbEE=[1.0]
BcEE=[0.0]
BEE=(BaEE,BbEE,BcEE)

################## General functions ###############
function it(A) #creating complex iteration array
    #inttyp=BSpline(Linear(Line(OnGrid())))
    inttyp=BSpline(Quadratic(Reflect(OnCell())))
    #inttyp=BSpline(Linear())
    matrdim=size(A,2)-1
    step=abs(real(A[end,1]-A[1,1]))/(size(A,1)-1)
    scaleitp=real(A[1,1]):step:real(A[end,1])
    ARe=real(A); AIm=imag(A);
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


function sub(it,t) #substitution in interpolating function
    subdim=size(it,1)
    out=zeros(ComplexF64,subdim)
    for j = 1:subdim
        out[j]=it[j,1](t)+it[j,2](t)*im
    end
    return(out)
end

function normmax(a) #finding the greatest eigenvalue
    s=maximum(size(a))
    norma=zeros(s)
    for j = 1:s
        norma[j]=norm(a[j])
    end
    return(a[findmax(norma)[2]])
end

function normmin(a) #finding the smallest eigenvalue
    s=maximum(size(a))
    norma=zeros(s)
    for j = 1:s
        norma[j]=norm(a[j])
    end
    return(a[findmin(norma)[2]])
end

function infnancheck(a) #checking inversion feasability
    var=false
    if any(isnan,a)==true || any(isinf,a)==true
    var=true
    end
return(var)
end

function evplot(a,k) #eigenvalue plot with circle
    U=a[:,k]
    eigRe1=[real(U[n,1]) for n in 1:size(U,1)]
    eigIm1=[imag(U[n,1]) for n in 1:size(U,1)]
    plot1=Plots.plot(eigRe1,eigIm1,seriestype=:scatter)
    Re_circle=[cos(n) for n in 0:0.05:2*pi]
    Im_circle=[sin(n) for n in 0:0.05:2*pi]
    return(Plots.plot!(plot1,Re_circle,Im_circle))
end

function test(gmaxtest1,ntest1,multtest1,rep)
    valref=valrefA
    gmaxn=size(gmaxtest1)[1]
    ntestn=size(ntest1)[1]
    multtestn=size(multtest1)[1]
    table=zeros(Float64,multtestn*ntestn*gmaxn,5)

    tabtempcomp=zeros(Float64,rep)
    tabtemperr=zeros(ComplexF64,rep)

    for j1=1:multtestn
        global mult=multtest1[j1]
        for j2=1:gmaxn
            global gmax=gmaxtest1[j2]
            for j3=1:ntestn
                global n=ntest1[j3]
                for j4=1:rep
                    valtemp=@timed ISIM(v)
                    tabtemperr[j4]=normmax((valtemp[1])[:,end])
                    tabtempcomp[j4]=valtemp[2]
                end
                vecttempabs=real.(tabtemperr)+im*abs.(imag.(tabtemperr))
                tableerr=real(mean(abs.(vecttempabs-valref*ones(ComplexF64,rep,1))))
                tablecomp=mean(tabtempcomp)
                table[(j1-1)*ntestn*gmaxn+(j2-1)*ntestn+j3,:]=[convert(Float64,multtest1[j1]) convert(Float64,gmaxtest1[j2]) convert(Float64,ntest1[j3]) tablecomp tableerr]
            end
        end
    end
    return(table)
end

function iter(S1,V1)
    n1=size(S1)[1]
    mult1=trunc(Int,size(S1)[2]/dim)
    S=zeros(ComplexF64,n1*dim,mult1)
    V=zeros(ComplexF64,n1*dim,mult1)
    for p=1:n1
           for s=1:mult1
               for q=1:dim
                   S[1+(p-1)*dim+(q-1),s]=S1[p,1+(s-1)*dim+(q-1)]
                   V[1+(p-1)*dim+(q-1),s]=V1[p,1+(s-1)*dim+(q-1)]
               end
           end
       end
       H=pinv(S)*V #pseudo-inverse calculation
       eigH=eigen(H)
       Hval=eigH.values #eigenvalue calculation
       Vj=V*eigH.vectors #calculating of new set of eigenvectors
       Sj=zeros(ComplexF64,n1,mult1*dim) #creating new initial solution array
       for p=1:n1
           for s=1:mult1
               for q=1:dim
                        Sj[p,1+(s-1)*dim+(q-1)]=Vj[1+dim*(p-1)+(q-1),s]
               end
           end
       end
       return((Sj,Hval))
end

###################### ISIM #######################

function ISIM(v1)
        dt=tau(v1)/(n-1) #timestep
        nmax=floor(Int,round((T(v1)/dt)))
        kint=floor(Int,(T(v1)/tau(v1)))
        nrest=nmax-kint*(n-1)
        tvec=collect(-tau(v1):dt:(nmax*dt))
        sol00=randn!(rng, zeros(ComplexF64,n,mult*dim))
        sol=zeros(ComplexF64,nmax,mult*dim)  #empty solution matrix
        sol0m=zeros(ComplexF64,n,mult*dim)
        sol0=zeros(ComplexF64,size(tvec)[1],mult*dim+1)
        Hval0=zeros(ComplexF64,mult)

        for g=1:gmax
            sol0=hcat(tvec,vcat(sol00,sol))

            if method == "Julia"
                solarr=fsol(tau(v1),sol0[1:n+1,:],Ai(v1),Bi(v1),nmax*dt,mult,dt)

                for tv=0:n-1
                    for j=1:mult*dim
                        sol0m[tv+1,j]=solarr((nmax*dt-tau(v1))+tv*dt)[j]
                    end
                end

            elseif method == "RK4"
                for k=1:kint
                interp=it(sol0[1+(k-1)*(n-1):n+(k-1)*(n-1)+1,:])
                    for j=1:(n-1)
                        sol0[n+j+(k-1)*(n-1),2:end]=transpose(butcher(real(sol0[n+(j-1)+(k-1)*(n-1),1]),interp,sol0[n+(j-1)+(k-1)*(n-1),2:end],dt,BRK4,v1,tau(v1),mult))
                    end
                end
                if nrest>0
                    interp=it(sol0[1+(kint-1)*(n-1):n+(kint-1)*(n-1),:])
                    for j=1:nrest
                        sol0[n+j+(kint-1)*(n-1),2:end]=transpose(butcher(real(sol0[n+j+(kint-1)*(n-1),1]),interp,sol0[n+j-1+(kint-1)*(n-1),2:end],dt,BRK4,v1,tau(v1),mult))
                    end
                end
                sol0m=sol0[end-(n-1):end,2:end]

            end
            resit=iter(sol0[1:n,2:mult*dim+1],sol0m)
            sol00=resit[1]
            Hval0=hcat(Hval0,resit[2])
        end
        print(gmax)
        return(Hval0[:,2:end])
end

gmax=16
mult=32
n=2002
method="RK4"
valrefAn=normmax((ISIM(v))[:,end])
valrefA=real(valrefAn)+im*abs(imag.(valrefAn))

gmaxtest=collect(2:1:8)
ntest=collect(100:100:1500)
multtest=collect(4:1:12)

text=test(gmaxtest,ntest,multtest,2)

open("temp.txt", "w") do io
            writedlm(io,text)
       end

n=50
mult=6
gmax=10
####### Multi-Dimensional Bisection Method #########
ISIM(v)[end,:]

fitc=ISIM(v)
fitc[2,2]

function foo(x,y)
    return(norm(normmax(ISIM([omega kappa x epsilon y tau0])[:,end]))-1)
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
