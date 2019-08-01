function OGLE_file_string(index_::Int64)::String
    return "./data/OGLE-LMC-CEP-"*lpad(index_,4,"0")*".dat"
end

"""
Read a table of floats from file. If dimensions mismatch, returns zeros(3,3). Can deal with "nan" and "inf" on the float table.
"""
function read_float_table(fname_::String)::Array{Float64,2}
    try 
        return parse.(Float64,(hcat(split.(readlines(fname_))...)))
    catch err
        if isa(err,DimensionMismatch)
            #println("Read error at file $fname_: not a float table. Returning `zeros(3,3)`")
            return zeros(3,3)
        elseif isa(err,ArgumentError)
            println("ArgumentError: "*err.msg*". Returning `zeros(3,3)`")
            return zeros(3,3)
        else
            println("Unhandled error on read_float_table(\"$fname_\"):")
            println(err)
        end
    end
end

















