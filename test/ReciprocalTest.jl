module ReciprocalTest
using FromFile
using Plots
@from "../utils/Reciprocal.jl" using Reciprocal

# THIS IS NOT A TRUE TEST MODULE AND AS SUCH DOES NOT APPEAR IN runtests.jl
function discretehexagontest(eps::Float64, method::String)
    points = discretehexagon(eps, method)
    x = []
    y = []
    for point in points
      push!(x, point[1])
      push!(y, point[2])
    end
    plot = scatter(x, y)
    savefig(plot, "test.png")
end

discretehexagontest(0.05, "constantstep")
end