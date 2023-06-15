module Rotations
export rotate2d

function rotate2d(v::Vector{Float64}, rad::Float64)
    if length(v) != 2
        error("This function is only intended for rotations of 2d vectors.")
    end
    rot = [cos(rad) -sin(rad); sin(rad) cos(rad)]
    rotv = rot * v
    return rotv
end

end