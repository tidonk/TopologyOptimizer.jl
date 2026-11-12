# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                           #
#  Copyright 2026 Operations Research, RWTH Aachen University               #
#                                                                           #
#  Licensed under the Apache License, Version 2.0 (the "License");          #
#  you may not use this file except in compliance with the License.         #
#  You may obtain a copy of the License at                                  #
#                                                                           #
#      http://www.apache.org/licenses/LICENSE-2.0                           #
#                                                                           #
#  Unless required by applicable law or agreed to in writing, software      #
#  distributed under the License is distributed on an "AS IS" BASIS,        #
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. #
#  See the License for the specific language governing permissions and      #
#  limitations under the License.                                           #
#                                                                           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

using LinearAlgebra: eigen, Symmetric, I

"classical MDS: project a distance matrix to 2D coordinates."
function classical_mds_2d(dist_matrix::Matrix{Float64})::Matrix{Float64}
    n = size(dist_matrix, 1)
    D2 = dist_matrix .^ 2
    H = Matrix{Float64}(I, n, n) .- (1.0/n)
    B = -0.5 .* H * D2 * H
    B_sym = Symmetric((B .+ B') ./ 2)  # ensure symmetry
    eig = eigen(B_sym)
    # take two largest eigenvalues (at end since eigen sorts ascending)
    idx = sortperm(eig.values, rev=true)
    coords = zeros(n, 2)
    for d in 1:2
        λ = max(eig.values[idx[d]], 0.0)  # clamp negative eigenvalues
        coords[:, d] = eig.vectors[:, idx[d]] .* sqrt(λ)
    end
    return coords
end

"write an interactive HTML scatter plot of the cut clustering to `path`."
function write_cut_cluster_html(path::String,
                                dist_matrix::Matrix{Float64},
                                clusters::Vector{Int},
                                violations::Vector{Float64},
                                coeff_matrix::Matrix{Float64},
                                chosen::BitVector;
                                casename::String=_CASENAME)
    coords = classical_mds_2d(dist_matrix)
    n = size(coeff_matrix, 1)
    n_vars = size(coeff_matrix, 2)

    # build per-point data
    xs = coords[:, 1]
    ys = coords[:, 2]
    texts = String[]
    for i in 1:n
        nz = count(abs.(coeff_matrix[i, :]) .> SOLVER_PRECISION)
        pct = round(100.0 * nz / n_vars, digits=1)
        push!(texts, "Cut $i<br>Nonzero: $nz / $n_vars ($pct%)<br>Violation: $(round(violations[i], digits=4))<br>Cluster: $(clusters[i])<br>Chosen: $(chosen[i])")
    end

    # group by (cluster, chosen) for separate Plotly traces
    unique_clusters = sort(unique(clusters))
    traces_json = String[]
    palette = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"]
    for cl in unique_clusters
        for ch in [true, false]
            mask = (clusters .== cl) .& (chosen .== ch)
            if !any(mask)
                continue
            end
            idxs = findall(mask)
            color = palette[mod1(cl, length(palette))]
            marker_sym = ch ? "circle" : "x"
            marker_size = ch ? 10 : 7
            opacity = ch ? 1.0 : 0.5
            name = "Cluster $cl" * (ch ? " (chosen)" : "")
            x_arr = join(["$(xs[i])" for i in idxs], ",")
            y_arr = join(["$(ys[i])" for i in idxs], ",")
            t_arr = join(["\"$(texts[i])\"" for i in idxs], ",")
            push!(traces_json, """{
                x: [$x_arr], y: [$y_arr],
                text: [$t_arr], mode: 'markers',
                marker: {color: '$color', symbol: '$marker_sym', size: $marker_size, opacity: $opacity},
                name: '$name', hoverinfo: 'text', type: 'scatter'
            }""")
        end
    end

    title_text = isempty(casename) ? "Cut Clustering (MDS projection)" : "Cut Clustering — $casename (MDS projection)"

    html = """<!DOCTYPE html><html><head><meta charset="utf-8">
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script></head>
    <body><div id="plot" style="width:100%;height:95vh;"></div><script>
    Plotly.newPlot('plot', [$(join(traces_json, ",\n"))], {
        title: '$title_text',
        xaxis: {title: 'MDS dim 1'}, yaxis: {title: 'MDS dim 2'},
        hovermode: 'closest'
    });</script></body></html>"""

    open(path, "w") do f
        write(f, html)
    end
    @info("Cut cluster visualization written to $path")
end
