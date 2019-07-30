"""
    loc_max(arr::Array{Float64,1})

Shorcut to Julia's findmax.
"""
function loc_max(arr::Array{Float64,1})
    return findmax(arr)[2]
end


"""
    Complex_NFT(x_arr::Array{Float64,1},y_arr::Array{Float64,1},t_::Float64)

Fastest implementation of type III NUDFT that I could think of right now.

Calculates the square of the magnitude of the fourier transform of data secuence (x_arr,y_arr) at PERIOD t_

I have disabled all possible security checks or assumptions on arrays, namely:

- function is inlined
- loops can be performed in arbitrary order, so every cicle is independent of the others.
- never waits on a previous iteration to make forward progress
- Does not check for memory dependences on iterations
- only unit increments on loop index
- non nested loops

And have the following numeric assumptions
- no nan values
- no inf values
- no signed zeros
- no division, using reciprocal multiplication  (a/b -> a*1/b)
- fused operations allowed
- common math functions are approximated when convenient (sin, exp, log, etc)
- re-association of associative operations is allowed (!)

Benchmark results with same data (300 random numbers each) on 1000 random periods, 100 evals/period:

    memory estimate:  0 bytes
    allocs estimate:  0
    --------------
    minimum time:     9.149  μs (0.00% GC)
    median time:      9.389  μs (0.00% GC)
    mean time:        10.179 μs (0.00% GC)
    maximum time:     36.571 μs (0.00% GC)

and mapping over an array of 1000 random periods, 5 runs, 100 evals/run:

    memory estimate:  7.95 KiB
    allocs estimate:  2
    --------------
    minimum time:     10.768 ms (0.00% GC)
    median time:      11.237 ms (0.00% GC)
    mean time:        11.173 ms (0.00% GC)
    maximum time:     11.529 ms (0.00% GC)

The benchmraks does count random number generation time.

"""
@inline @fastmath function Complex_NFT(x_arr::Array{Float64,1},y_arr::Array{Float64,1},t_::Float64)
    s::ComplexF64=0.0+0.0im
    @simd ivdep for i=1:length(x_arr)
        @inbounds s+=y_arr[i]*exp(-6.283185307179586im*x_arr[i]/t_)
    end
    return abs2(s)
end


"""
    find_period_ephemeris(X_::Array{Float64,1},Y_::Array{Float64,1},T_::{Float64:Float64:Float64})

Attempt to find a period and ephemeris of the signal data (X_,Y_) on the Julia range T_ via maximum of fourier transform modulus.

The ephemeris is found with the formula: period*(arg(F)/2pi+0.25), where the 0.25 counts for displaying the phase as a sin wave.

Returns a tuple (period,ephemeris)
"""
@inline function find_period_ephemeris(X_::Array{Float64,1},Y_::Array{Float64,1},T_::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}})
    F::Array{ComplexF64}=map(t->Complex_NFT(X_,Y_,t),T_)
    i=loc_max(abs2.(F))
    return T_[i], T_[i]*(angle(F[i])/2pi+0.25)
end


"""
    normalize101!(arr_::Array{Float64})

Normalize the array in-place to mean=0 and std=1

Is used to feed find_period_ephemeris function; the normalized magnitude eliminates period aliasing due to cuasi-constant sample on the time series due to phase cancellation on the fourier transform sum.
"""
function normalize101!(arr_::Array{Float64})
    arr_.= (arr_.-mean(arr_))./std(arr_)
end






