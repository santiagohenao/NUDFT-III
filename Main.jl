using Statistics
using Dates

include("FileOperations.jl")
include("MathFunctions.jl")

function run(star_id::Int64,dt::Float64=10e-4,refine::Float64=10e-3)
    local start::Float64=time()
    local JD::Array{Float64,1}
    local M::Array{Float64,1}
    local err::Array{Float64,1}
    local per::Float64
    local eph::Float64
    local Q::Float64 # period quality factor
    
    # import data
    JD,M,err=star_id |> OGLE_file_string |> read_float_table |> permutedims |> x->[x[:,i] for i in 1:3]
    
    # chech for gibberish
    if JD==[0,0,0]
        return (star_id,0,0)
    end
    
    # normalize magnitude
    normalize101!(M)
    
    # rough search
    per,eph,Q=find_period_ephemeris(JD,M,0.1:dt:20.0)
    
    # refined search
    per,eph,Q=find_period_ephemeris(JD,M,per-2dt:dt*refine:per+2dt)
    
    return (star_id,r5(per),r5(eph),r5(Q),r5(time()-start))
end

F_numbers=parse.(Int64,(x->x[14:17]).(readdir("./data")))



for i=F_numbers[1:100]
    println(run(i))
end

