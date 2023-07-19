module Reciprocal
using FromFile
using LinearAlgebra: norm
using Plots
@from "Rotations.jl" using Rotations

export reciprocalgraphene, 
       discretehexagon,
       getgraphenepath

function reciprocalgraphene()
    # Assuming lattice constant of 1
    k1 = zeros(2)
    k1[1] = 1
    k1[2] = sqrt(3)
    k1 = 2 * pi / sqrt(3) * k1
    k2 = zeros(2)
    k2[1] = 1
    k2[2] = -sqrt(3)
    k2 = 2 * pi / sqrt(3) * k2

    return (k1, k2)
end

# CONSTRUCTING HIGH SYMMETRY PATH
function getgraphenepath()
    k1, k2 = reciprocalgraphene()
  
   # from bottom face to origin
    start = k2 / 2
    stop = zeros(2)
    steps = 300
    delta = (stop - start) / steps
    path1 = [start + i*delta for i in 1:steps]
  
    # from origin to corner
    start = stop
    rot90 = [0 -1; 1 0]
    stop = 1 / sqrt(3) * rot90 * k2
    steps = 300
    delta = (stop - start) / steps
    path2 = [start + i * delta for i in 1:steps]
  
    # from corner to corner
    start = stop
    stop = start + 1 / sqrt(3) * norm(k2) * [0; -1]
    steps = 300
    delta = (stop - start) / steps
    path3 = [start + i * delta for i in 1:steps]
  
    # from corner back to origin
    start = stop
    stop = zeros(2)
    steps = 300
    delta = (stop - start) / steps
    path4 = [start + i * delta for i in 1:steps]
  
    # from origin to top face
    start = stop
    stop = k1 / 2
    steps = 300
    delta = (stop - start) / steps
    path5 = [start + i * delta for i in 1:steps]
  
    points = cat(path1, path2, path3, path4, path5, dims = 1)
    return points
  end

function discretehexagon(eps::Float64, method::String)

    _, k2 = reciprocalgraphene()
    points = []
    tpoints = []

    # first calculate the bottom most corner points of the rightmost equlitaral triangle
    # from the reciprocal lattice vectors.
    corner = rotate2d(k2 / 2, pi / 6)

    #=
    Now we identify discrete and more or less evenly spaced points within this region.
    it is important here to make sure that the points are more or less evenly spaced
    so that we do not oversample any one region in the brillouin zone, this would harm
    the convergence of riemann sums in this region to their integral values.
    =#

    if lowercase(method) == "bisect"
        lowerbound = recursivebisect([0.0; 0.0], corner, eps)
    elseif lowercase(method) == "constantstep"
        lowerbound = constantstep([0.0; 0.0], corner, eps)
    end
    upperbound = collect(Iterators.map(point -> rotate2d(point, pi / 3), lowerbound))

    push!(tpoints, lowerbound[2:end]...)
    for (ub, lb) in zip(upperbound[2:end], lowerbound[2:end])
        if lowercase(method) == "bisect"
            ps = recursivebisect(ub, lb, eps)
        elseif lowercase(method) == "constantstep"
            ps = constantstep(ub, lb, eps)
        end
        push!(tpoints, ps[2:end-1]...) # avoiding adding lowerbound and upperbound twice
    end

    # Rotating every point in this triangle by pi / 3 radians gives the northeast triangle
    # repeating this four more times will give us the discretized points within the regions
    # encompassed by the NW, W, SW, SE triangles.
    for n in 0:5
        for point in tpoints
            push!(points, rotate2d(point, n * pi / 3))
        end
    end

    push!(points, [0.0; 0.0])
    return points
end

function constantstep(startpoint, endpoint, stepsize)
    if norm(endpoint - startpoint) < stepsize
        return [startpoint, endpoint]
    end

    points = [startpoint]
    connector = (endpoint - startpoint) / norm(endpoint - startpoint)

    i = 1
    while true
        newpoint = points[i] + stepsize * connector
        norm(newpoint - endpoint) > stepsize || break
        push!(points, newpoint)
        i+=1
    end
    push!(points, endpoint)
    return points
end

function recursivebisect(startpoint, endpoint, epsilon)
    delta = norm(startpoint - endpoint)
    pointsin = [startpoint, endpoint]
    pointsout = []
    while delta > epsilon
        pointleft = [0; 0]
        pointright = pop!(pointsin)
        while length(pointsin) != 0
            pointleft = pop!(pointsin)
            pointmid =  pointleft + (pointright - pointleft) / 2
            pushfirst!(pointsout, pointright)
            pushfirst!(pointsout, pointmid)
            pointright = pointleft
        end
        pushfirst!(pointsout, pointleft)

        pointsin, pointsout = pointsout, pointsin
        delta = norm(pointsin[1] - pointsin[2])
    end
    return pointsin
end

end