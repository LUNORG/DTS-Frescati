using Dates
using LibPQ
using IniFile
using DataFrames
using Plots
using Dierckx
using LinearAlgebra
using JSON
using CSV
using Statistics

x = [10, 20, 30, 40, 50, 60];
s = 10;
e = 60;

function mirror_bh(x, s, e)
    y_vec = []
    for i in x
        if i-s<(e-s)/2
            y=s-i
        else
            y=i-e
        end
        push!(y_vec, y)
    end
return y_vec
end

testy = mirror_bh(x, s, e)

