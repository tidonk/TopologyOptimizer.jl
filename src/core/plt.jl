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

using Colors, ColorSchemes, Plots, PowerPlots, StatsPlots, VegaLite, PlotlyJS

function get_load_percentage(branch)
  if get(branch, "br_status_new", branch["status"]) > 0
    if ! haskey(branch, "pf_new")
      @warn("Branch $(branch["name"]) has no pf_new: \n $branch")
      return 0.0
    end
    p = abs(branch["pf_new"])

    p_max = branch["rate_a"]
    # due to slacks on rate_a, this can be larger than 1
    return p/p_max
  else
    return 0.0
  end
end

function make_flow_plot(pm_data)
  scheme = palette(ColorSchemes.:Reds_9, convert(Int64,round(length(pm_data["branch"])*2)))

  branch_contingencies = [branch_indices for branch_indices in [get(contingency, :branch, []) for contingency ∈ pm_data["contingencies"]]]
  #branches_in_contingencies = reduce(append!, [branch_indices for branch_indices in [get(contingency, :branch, []) for contingency ∈ pm_data["contingencies"]]])

  # define branches to observe
  b_ids = []
  # fbus and tbus of these branches
  b_buses = union([pm_data["branch"][string(b_id)]["f_bus"] for b_id in b_ids], [pm_data["branch"][string(b_id)]["t_bus"] for b_id in b_ids])

  if ! isempty(b_ids)
    # print all branches 
    for (_, branch_data) ∈ pm_data["branch"]
      if branch_data["index"] ∈ b_ids
        @show(branch_data)
      end
    end
    
    # print all gens that are connected to these branches
    for (_, gen_data) ∈ pm_data["gen"]
      if gen_data["gen_bus"] ∈ b_buses
        @show(gen_data)
      end
    end
  
    # print all loads that are connected to these branches
    for (_, load_data) ∈ pm_data["load"]
      if load_data["load_bus"] ∈ b_buses
        @show(load_data)
      end
    end
  end

  # calculate additional branch data
  for (br_id, branch) ∈ pm_data["branch"]
    if get(branch, "br_status_new", 0) > 0 || pm_data["BUSBARS"] 
      if haskey(branch, "pf_new")
        branch["load_percentage"] = get_load_percentage(branch)
        branch["overloaded"] = branch["load_percentage"] == 1 ? 1 : 0
        if haskey(branch, "pf")
          branch["orig_sol_diff"] = abs(abs(branch["pf_new"]) - abs(branch["pf"]))
          branch["orig_sol_diff_binary"] = branch["orig_sol_diff"] > 0.1 ? 1 : 0
        else
          branch["orig_sol_diff"] = 0.0
          branch["orig_sol_diff_binary"] = 0
        end
      end
    end
    #branch["in_contingencies"] = branch["index"] ∈ branches_in_contingencies
    branch["contingency"] = [contingency for contingency in pm_data["contingencies"] if branch["index"] ∈ get(contingency, :branch, [])]
    branch["in_contingencies"] = length(branch["contingency"]) ≠ 0 ? 1 : 0
    branch["monitored"] = Dict(:branch => [branch["index"]]) ∈ get(pm_data, "monitoring", [])
    if branch["index"] ∈ b_ids
      branch["marking"] = 1
    else
      branch["marking"] = 0
    end

    if :ratea in pm_data["slack_types"]
      if haskey(branch, "raf_s_new") && haskey(branch, "rat_s_new")
        branch["slack_used"] = branch["raf_s_new"] < SOLVER_PRECISION && branch["rat_s_new"] < SOLVER_PRECISION ? 0 : 1
      end
    end
  end
  
  for (l_id, load) ∈ pm_data["load"]
    if :pd in pm_data["slack_types"]
      load["slack_used"] = load["pd_slack_lower"] < SOLVER_PRECISION && load["pd_slack_upper"] < SOLVER_PRECISION ? 0 : 1
    end
    if load["index"] ∈ []
      load["marking"] = 1
    else
      load["marking"] = 0
    end
  end

  for (g_id, gen) ∈ pm_data["gen"]
    if :pg in pm_data["slack_types"]
      gen["slack_used"] = gen["pmin_slack"] < SOLVER_PRECISION && gen["pmax_slack"] < SOLVER_PRECISION ? 0 : 1
    end
    if gen["index"] ∈ []
      gen["marking"] = 1
    else
      gen["marking"] = 0
    end
  end

  for (b_id, bus) ∈ pm_data["bus"]
    bus["connected_lines"] = Dict(branch_data["index"] => round(get(branch_data, "status", -1), digits=2) for branch_data in [get_branch_by_index(pm_data, branch) for branch in get_connected_branches(pm_data, bus["index"])])
    if bus["index"] ∈ b_buses
      bus["marking"] = 1
    else
      bus["marking"] = 0
    end
  end


  plt = powerplot(pm_data, width=3840, height=2160,
              bus_size=5, gen_size=0, load_size=0, branch_size=1, connector_size=0,
              bus_data="marking",
              bus_data_type="nominal",
              bus_color = ["purple", "green"],
              branch_data="orig_sol_diff_binary", # "marking", 
              branch_data_type="quantitative",
              #branch_color= scheme,
              branch_color=colorscheme2array(ColorSchemes.:Reds_9)[2:end],
              gen_data="pg",
              gen_data_type="quantitative",
              gen_color=colorscheme2array(ColorSchemes.:Purples)[4:end],
              load_data="pd",
              load_data_type="quantitative",
              load_color=colorscheme2array(ColorSchemes.:Reds)[1:end-3],
              show_flow=false,
              fixed=true)
  return plt
end

function make_switching_plot(pm_data)
  plt = powerplot(pm_data, width=2480/2, height=3508/2.5,
              bus_size=50, gen_size=40, load_size=40, branch_size=1.7, connector_size=2,
              branch_data="br_status_new", 
              branch_data_type="nominal",
              branch_color= ["yellow", "grey", "orange"],
              gen_data="pg_new",
              gen_data_type="quantitative",
              gen_color=colorscheme2array(ColorSchemes.:Purples)[4:end],
              load_data="pd",
              load_data_type="quantitative",
              load_color=colorscheme2array(ColorSchemes.:Reds)[4:end],
              show_flow=false,
              fixed=true
  )
  return plt
end

function make_layout_plot(pm_data)
  for (b_id, bus) ∈ pm_data["bus"]
    bus["connected_lines"] = Dict(branch_data["index"] => round(get(branch_data, "rate_a", -1), digits=2) for branch_data in [get_branch_by_index(pm_data, branch) for branch in get_connected_branches(pm_data, bus["index"])])
    if bus["index"] ∈ []
      bus["marking"] = 1
    else
      bus["marking"] = 0
    end
  end


  plt = powerplot(pm_data, width=2480, height=3408,
              bus_size=10, gen_size=0, load_size=0, branch_size=5.0, connector_size=0,
              bus_data="marking",
              bus_data_type="quantitative",
              bus_color = ["green", "red"],
              branch_data="type", 
              branch_data_type="nominal",
              branch_color= ["blue", "orange", "purple", "orange", "purple"],
              gen_data="pg",
              gen_data_type="quantitative",
              gen_color=colorscheme2array(ColorSchemes.:Purples)[4:end],
              load_data="pd",
              load_data_type="quantitative",
              load_color=colorscheme2array(ColorSchemes.:Reds)[4:end],
              show_flow=false,
              fixed=true
  )
  return plt
end


function make_plots(data, result, case; dir=CONSTS[:USER_ENVIRONMENT][:DIRECTORIES][:PLOT_DIR])
# make plots
  if haskey(data, "nw")
    data = data["nw"]["1"]
    #result["solution"] = result["solution"]["nw"]["1"] 
  end
  # create dir to save plots to
  mkpath(dir)

  # fix layout to coordinates

  if any([haskey(bus, "Lat") && haskey(bus, "Lon") for (b_id, bus) ∈ data["bus"]])
    for (b_id, bus) ∈ data["bus"]
      if bus["Lat"] == 0.0 && bus["Lon"] == 0.0
        #@error("Bus $(bus["index"]) has coordinates (0,0)")
        lat = 0
        lon = 0
        n_adjacent_buses = 0

        # match to closest bus by interpolated coordinates of all adjacent buses
        for (br_id, branch) ∈ data["branch"]
          if branch["f_bus"] == bus["index"]
            if data["bus"][string(branch["t_bus"])]["Lat"] == 0.0 && data["bus"][string(branch["t_bus"])]["Lon"] == 0.0
              continue
            end
            n_adjacent_buses += 1
            lat += data["bus"][string(branch["t_bus"])]["Lat"]
            lon += data["bus"][string(branch["t_bus"])]["Lon"]
          elseif branch["t_bus"] == bus["index"]
            if data["bus"][string(branch["f_bus"])]["Lat"] == 0.0 && data["bus"][string(branch["f_bus"])]["Lon"] == 0.0
              continue
            end
            n_adjacent_buses += 1
            lat += data["bus"][string(branch["f_bus"])]["Lat"]
            lon += data["bus"][string(branch["f_bus"])]["Lon"]
          end
        end
        
        bus["xcoord_1"] = lat / n_adjacent_buses
        bus["ycoord_1"] = lon / n_adjacent_buses
      else
        bus["xcoord_1"] = bus["Lat"]
        bus["ycoord_1"] = bus["Lon"]
      end
    end
  else
    data = layout_network(data; iterations=200)
  end
  # update to result data (flows, etc.)
  #update_data_safe!(data, result)
  
  # show angle deviations
  #println("Dict{String, Float64}(")
  #for (bus_id, bus) in data["bus"]
  #  println("$bus_id => $(bus["va_new"] - bus["va_start"]),")
  #  #@show(bus_id, bus["va_new"], bus["va_start"],  bus["va_new"] - bus["va_start"])
  #end
  #println(")")

  #plt = make_switching_plot(data)
  # make plot with coloring according to switching status determined in optimal solution
  if "br_status_new" in keys(data["branch"][first(keys(data["branch"]))])
    #plt = make_switching_plot(data)
    #VegaLite.save(dir * "$(case)_switching.pdf", plt)
  end
  # make plot with coloring according to flow
  #plt = make_flow_plot(data)
  #VegaLite.save(dir * "$(case)_flow.pdf", plt)

  #PowerPlots.Experimental.add_zoom!(plot)


  # make plot illustrating the layout
  plt = make_flow_plot(data)
  #plt = make_layout_plot(data)
  VegaLite.save(dir * "$(case)_layout.html", plt)

  return plt
end

function runtime_comparison_plot(runtimes::Dict{String, Vector{Float64}})
  nam = repeat(collect(keys(runtimes)), outer = length(runtimes[collect(keys(runtimes))[begin]]))
  ctg = repeat(["SCOTS"], inner = length(collect(keys(runtimes))))
  runtimes = reduce(vcat,transpose.([values(runtimes)...]))
  groupedbar(nam, runtimes, group = ctg, xlabel = "Case", ylabel = "Runtime [s]",
          title = "Runtimes of DC Formulations", bar_width = 0.67,
          lw = 0, framestyle = :box, yaxis=:log)
end