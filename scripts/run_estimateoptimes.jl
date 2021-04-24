#!/usr/bin/env julia


using Random: rand
using SpecialFunctions: besselj0, besselj1


function plus(a::Float64, b::Float64)::Float64
    return a + b
end

function mult(a::Float64, b::Float64)::Float64
    return a * b
end

function divi(a::Float64, b::Float64)::Float64
    return a / b
end

function root(a::Float64)::Float64
    return sqrt(a)
end


function timing(N::Int, print::Bool=true)::Nothing

    data = rand(Float64, N) * 1e4

    plustime = @elapsed plus.(data, Float64(pi))
    multtime = @elapsed mult.(data, Float64(pi))
    divitime = @elapsed divi.(data, Float64(pi))
    expntime = @elapsed exp.(data)
    roottime = @elapsed root.(data)
    sinetime = @elapsed sin.(data)
    bes0time = @elapsed besselj0.(data)
    bes1time = @elapsed besselj1.(data)

    if print
        println("="^40)
        println("Addition time: $(plustime / N)")
        println("Multiplication time: $(multtime / N)")
        println("Division time: $(divitime / N)")
        println("Exponentiation time: $(expntime / N)")
        println("Square root time: $(roottime / N)")
        println("Sine time: $(sinetime / N)")
        println("Bessel 0 time: $(bes0time / N)")
        println("Bessel 1 time: $(bes1time / N)")
        println("-"^20)
        println("Mult cost: $(multtime / multtime)")
        println("Divi cost: $(divitime / multtime)")
        println("Expn cost: $(expntime / multtime)")
        println("Root cost: $(roottime / multtime)")
        println("Sine cost: $(sinetime / multtime)")
        println("Bes0 cost: $(bes0time / multtime)")
        println("Bes1 cost: $(bes1time / multtime)")
        println("="^40)
    end

end


timing(4,false)
for i in 1:3
    timing(10000)
end


