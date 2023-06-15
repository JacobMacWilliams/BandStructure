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
    savepath = joinpath("test", "plots", "discretehexagon.png")
    savefig(plot, savepath)
end

function plotpath()
  points = getgraphenepath()
  x = []
  y = []
  for p in points
      push!(x, p[1])
      push!(y, p[2])
  end
  scatter(x, y)
  savepath = joinpath("test", "plots", "graphenepath.png")
  savefig(savepath)
end

discretehexagontest(0.05, "constantstep")
end