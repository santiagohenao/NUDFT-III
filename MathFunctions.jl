"""
Shorcut to Julia's findmax.
"""
function loc_max(arr::Array{Float64,1})
    return findmax(arr)[2]
end


"""
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

Benchmark results with same data on 1000 random periods, 100 evals/period:

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
    s=0.0+0.0im
    @simd ivdep for i=1:length(x_arr)
        @inbounds s+=y_arr[i]*exp(-6.283185307179586im*x_arr[i]/t_)
    end
    return abs2(s)
end

