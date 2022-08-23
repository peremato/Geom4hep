using Virtual, BenchmarkTools
abstract type Animal end

struct Dog{Race} <: Animal end
struct Tiger <: Animal end
struct Duck <: Animal end
struct Dog2 <: Animal end
#struct Tiger2 <: Animal end
#struct Duck2 <: Animal end
#struct Dog3 <: Animal end
#struct Tiger3 <: Animal end
#struct Duck3 <: Animal end

const UAnimal = Union{Dog, Tiger, Duck, Dog2}

@virtual fast_func(x::Animal, y::Int) = error("No default method for score!")
@override fast_func(x::Dog{R}, y::Int) where R  = 2 + y
@override fast_func(x::Tiger, y::Int) = 3 + y
@override fast_func(x::Duck, y::Int) = 4 + y
fast_func(x::Dog2, y::Int) = 2 + y
#@override fast_func(x::Tiger2, y::Int) = 3 + y
#@override fast_func(x::Duck2, y::Int) = 4 + y
#@override fast_func(x::Dog3, y::Int) = 2 + y
#@override fast_func(x::Tiger3, y::Int) = 3 + y
#@override fast_func(x::Duck3, y::Int) = 4 + y

dyn_func(x::Animal, y::Int) = error("No default method for score!")
dyn_func(x::Dog, y::Int) = 2 + y
dyn_func(x::Tiger, y::Int) = 3 + y
dyn_func(x::Duck, y::Int) = 4 + y
dyn_func(x::Dog2, y::Int) = 2 + y
#dyn_func(x::Tiger2, y::Int) = 3 + y
#dyn_func(x::Duck2, y::Int) = 4 + y
#dyn_func(x::Dog3, y::Int) = 2 + y
#dyn_func(x::Tiger3, y::Int) = 3 + y
#dyn_func(x::Duck3, y::Int) = 4 + y

manual_func(x::Animal, y::Int) =
    if x isa Dog
        2 + y
    elseif x isa Tiger
        3 + y
    elseif x isa Duck
        4 + y
    elseif x isa Dog2
        2 + y
 #   elseif x isa Tiger2
 #       3 + y
 #   elseif x isa Duck2
 #       4 + y 
 #   elseif x isa Dog3
 #       2 + y
 #   elseif x isa Tiger3
 #       3 + y
 #   elseif x isa Duck3
 #       4 + y 
    else
        error("No default method for score!")
    end

const samples = (Dog(), Duck(), Tiger())
                # Dog2(), Duck2(), Tiger2(),
                # Dog3(), Duck3(), Tiger3())
animals = Animal[samples[rand(1:length(samples))] for i = 1:100]
uanimals = UAnimal[samples[rand(1:length(samples))] for i = 1:100]

function sum_score(score_func, xs::Vector{S}) where S
    s = 0
    for x in xs
        s += score_func(x, 3)
    end
    return s
end

@info "fast_func via Virtual.jl"
display(@benchmark(sum_score(fast_func, animals)))
@info "manual split by hand"
display(@benchmark(sum_score(manual_func, animals)))
@info "dyn_func by dynamic multiple dispatch"
display(@benchmark(sum_score(dyn_func, animals)))
@info "dyn_func by dynamic multiple dispatch using a Union"
display(@benchmark(sum_score(dyn_func, uanimals)))
