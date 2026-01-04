"""
    Par(n::Integer) -> Vector{Vector{Int}}
    Par(n::Integer, col::Integer) -> Vector{Vector{Int}}
    Par(n::Integer, row::Integer, col::Integer) -> Vector{Vector{Int}}

Generate all partitions of an integer `n`.

A partition of `n` is a non-increasing sequence of positive integers that sum to `n`.

# Arguments
- `n`: The integer to partition
- `row`: (optional) Maximum number of rows (default: `n`)
- `col`: (optional) Maximum number of columns (default: `n`)

# Examples
```jldoctest
julia> MOPS.Par(4)
5-element Vector{Vector{Int64}}:
 [4]
 [3, 1]
 [2, 2]
 [2, 1, 1]
 [1, 1, 1, 1]

julia> MOPS.Par(3)
3-element Vector{Vector{Int64}}:
 [3]
 [2, 1]
 [1, 1, 1]
```
"""
function Par(n::Integer)
    return Par_sub(n, n, n)
end

function Par(n::Integer, col::Integer)
    return Par_sub(n, n, col)
end

function Par(n::Integer, row::Integer, col::Integer)
    return Par_sub(n, row, col)
end

"""
    Par_sub(n::Integer, row::Integer, col::Integer) -> Vector{Vector{Int}}

Internal recursive function to generate partitions.

# Arguments
- `n`: The integer to partition
- `row`: Maximum number of rows
- `col`: Maximum number of columns
"""
function Par_sub(n::Integer, row::Integer, col::Integer)
    if n == 0
        return [Int[]]
    end
    
    if col == 0
        return Vector{Vector{Int}}()
    end
    
    result = Vector{Vector{Int}}()
    
    # Maple: i = -min(row,n) .. -iquo(col+n-1,col)
    # iquo(col+n-1, col) is integer division: (col+n-1) ÷ col
    # So i goes from -min(row,n) to -iquo(col+n-1,col)
    # Then we use -i as the actual part size
    
    i_min = -min(row, n)
    i_max = -div(col + n - 1, col)
    
    for i in i_min:i_max
        part_size = -i  # The actual part size (since i is negative)
        rest_n = n + i  # n - part_size, but i is negative so n + i = n - part_size
        
        # Recursively get partitions of the remainder
        rest_partitions = Par_sub(rest_n, -i, col - 1)
        
        # Prepend part_size to each partition
        for rest_part in rest_partitions
            new_partition = vcat([part_size], rest_part)
            push!(result, new_partition)
        end
    end
    
    return result
end

