#!/bin/bash

julia -t 24 chen_LE_detailed.jl
julia -t 24 rand_LE.jl
julia -t 24 rand_LE_detailed.jl
