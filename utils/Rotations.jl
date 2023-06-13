module Rotations
using Rotations

function rotate2d(v::Vector, rad::Float64)
    if length(v2d) != 2
        error("This function is only intended for rotations of 2d vectors.")
    end
    rot = [cos(rad) -sin(rad); sin(rad) cos(rad)]
    rotv = rot * v
    return rotv
end

end