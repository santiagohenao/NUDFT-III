# to run in parallel julia has to start with "julia -p auto" and the "JULIA_NUM_THREADS" env variable on the terminal has to be setted properly to the desired number of threads to run this program
# on powershell $env:JULIA_NUM_THREADS=4 is enough

@everywhere using Statistics
@everywhere using Dates

@everywhere include("FileOperations.jl")
@everywhere include("MathFunctions.jl")

@everywhere function run(star_id::Int64,dt::Float64=10e-4,refine::Float64=10e-3)
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
    if JD==[0.0,0.0,0.0]
        return (0.0,0.0,0.0)
    end
    
    # normalize magnitude
    normalize101!(M)
    
    # rough search
    per,eph,Q=find_period_ephemeris(JD,M,0.1:dt:20.0)
    
    # refined search
    per,eph,Q=find_period_ephemeris(JD,M,per-2dt:dt*refine:per+2dt)
    
    #println(star_id,"\tdone in\t",r5(time()-start))
    return (r5(per),r5(eph),r5(Q))
end

# read stars id from data filenames
const F_numbers=parse.(Int64,(x->x[14:17]).(readdir("./data")))
####### Compiling run function #########
#run(2)                                ##
##### does not take this seriously #####
global_start=time()
# declare an (empty) Channel for threads to look which star work on
const stars = RemoteChannel(()->Channel{Int}(32));
# declare an (empty) Channel of the type of run(star_id) to put results on
const results = RemoteChannel(()->Channel{Tuple}(32));
# function to find a period for all the stars
@everywhere function calculate_stars(stars,results)
    while true
        star_id = take!(stars)
        io = open("$(myid()).dat", append=true);
        s1=time()
        period,ephemeris,quality=run(star_id)
        write(io,"$(myid())\t$star_id\t$(round(time()-s1; digits=6))\t$period\t$ephemeris\t$quality\n")
        close(io)
        if period==0.0
            put!(results,("$(myid()), $star_id no data",1))
        else
            put!(results, ("$(myid()), $star_id done",0))
        end
    end
end

# fill (define) the stars Channel
function make_stars(n::Int64)
    for i in F_numbers[1:n]
        put!(stars,i)
    end
end

@everywhere n=400

# creates a task on the general queue
@async make_stars(n)


for p in workers() # start tasks on the workers to process requests in parallel
    remote_do(calculate_stars, p, stars, results)
end

start=time()

while n > 0
    s_=take!(results)
    println(s_[1])
    global n = n - 1
end

println("elapsed: ",time()-start)
#println("CPU: ",sum(times))
























