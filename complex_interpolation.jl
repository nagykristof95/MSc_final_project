using Printf
using Plots
using Statistics
using LinearAlgebra
using MDBM
using PyPlot
pygui(true);
using Interpolations

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

B=transpose([0 0.5 1 1.5 2; 3*im 2+4*im -5+1*im 6*im 5*im-4; 7*im -8+3*im -5+2*im 7*im 3*im-2])
real(B)
typeof(real(B))

iterB=it(B)

interpolate((real(B[:,1]),), real(B[:,2]),Gridded(Linear()))
#substitution in interpolating function
function sub(it,t)
    subdim=size(it,1)
    out=zeros(ComplexF64,subdim)
    for j = 1:subdim
        out[j]=it[j,1](t)+it[j,2](t)*im
    end
    return(out)
end

sub(iterB,1.25)
