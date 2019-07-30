using Statistics
using Dates

include("FileOperations.jl")
include("MathFunctions.jl")



JD,M,err=5 |> OGLE_file_string |> read_float_table |> permutedims |> x->[x[:,i] for i in 1:3]
normalize101!(M)
dt=10e-4
per0,eph0=find_period_ephemeris(JD,M,0.1:dt:10.0)
per,eph=find_period_ephemeris(JD,M,per0-2dt:dt/10^3:per0+2dt)

F_numbers=parse.(Int64,(x->x[14:17]).(readdir("./data")))



for i=F_numbers[1:100]
    start=time()
    JD,M,err=i |> OGLE_file_string |> read_float_table |> permutedims |> x->[x[:,i] for i in 1:3]
    if JD==[0,0,0]
        continue
    end
    normalize101!(M)
    dt=10e-4
    per0,eph0=find_period_ephemeris(JD,M,0.1:dt:10.0)
    per,eph=find_period_ephemeris(JD,M,per0-2dt:dt/10^3:per0+2dt)
    println("$i\t$(round(per;digits=5))\t\t$(round(time()-start;digits=5))")
end

