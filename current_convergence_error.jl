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
using LsqFit
@. expmodel(x,p)=p[1]+p[2]*exp(p[3]*x)

xvec=vec([1.0 2.0 3.0 4.0])
yvec=vec([3.0 3.1 3.15 3.16])

@time curve_fit(expmodel,xvec,yvec,[0,0.0,0.0]).param


rng = MersenneTwister(1234)

function sortev(Hval)
    s=size(Hval)[1]
    normHval=norm.(Hval)
    sortHval=hcat(Hval,normHval)
    sort!(sortHval, dims=2)
        return(sortHval)
end


function fitev(eigval0,mmax)
    eigval1=sort(norm.(eigval0),dims=1)
    (m1,g1)=size(eigval1)
    fitarray=Array{ExpFit{Float64}}(undef,(mmax,2))
    g1half=convert(Int,floor(g1/2))
    xvec=collect(1.0:g1)
    for m=1:mmax
        eigvsort=sortslices(transpose(hcat(xvec,eigval1[end-m+1,:])),dims=1,lt=(x,y)->isless(x[2],y[2]))
        eigvsortmin=eigvsort[2,1:g1half]
        eigvsortming=eigvsort[1,1:g1half]
        eigvsortmax=eigvsort[2,1+g1half:g1]
        eigvsortmaxg=eigvsort[1,1+g1half:g1]

        fitarray[m,1]=curve_fit(expmodel,eigvsortming,eigvsortmin,[0,0.0,0.0])
        fitarray[m,2]=curve_fit(expmodel,eigvsortmaxg,eigvsortmax,[0,0.0,0.0])
    end
    return(fitarray)
end

rndtest=randn!(rng, zeros(ComplexF64,4,6))

fitev(rndtest,4)
