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

# Regression test: partial decomposition with contingency PFs only (no base PF) in master.
# Mode = :ETA, N_MAX_CONTINGENCIES = 1, INCLUDE_BASE_PF = false. Warmstart enabled so
# η-scores are available for the selection step.

if ! isdefined(Main, :TopologyOptimizer) || ! isdefined(Main, :Test)
    include("../../common.jl")
end

function test_dc_scots_slack_benders_partial_eta_no_base(path_prefix)
    @testset verbose = true "LEQ Subproblem (PARTIAL_DECOMP: ETA, N=1, no base PF)" begin
        additional_params = Dict{Symbol, Any}(
            :WARMSTARTING => Dict{Symbol, Any}(:ENABLED => true,
                                               :WARMSTART_METHOD => :scopf),
            :BENDERS      => Dict{Symbol, Any}(
                :SUBPROBLEM_SOLVER_TYPE => TopologyOptimizer.LinEqSysSolver,
                :PARTIAL_DECOMP => Dict{Symbol, Any}(
                    :ENABLED              => true,
                    :INCLUDE_BASE_PF      => false,
                    :MODE                 => :OVERLOAD,
                    :N_MAX_CONTINGENCIES  => 1,
                    :N_MAX_CONTINGENCIES_PCT => 0.0,
                    :MIN_VIOLATION        => 0,
                    :WEIGHTS => Dict{Symbol, Any}(:RND => 0.0, :FULL => 0.0, :SENS => 0.0,
                                                  :CONGPROX => 0.0, :ETA => 1.0)
                )
            )
        )

        run_multi_test(CONSTS[:TEST_SETS][:PIPELINE_TESTSET][1:1],
                       path_prefix       = path_prefix,
                       application       = "DC_SCOTS_TIME_SLACK_BENDERS",
                       additional_params = additional_params,
                       test_default      = true,
                       test_busbars      = false,
                       test_reconf       = false,
                       test_gen_sw       = false,
                       test_name         = "DC_SCOTS_BENDERS_PARTIAL_ETA_NO_BASE",
                       validate          = false)
    end
end

test_dc_scots_slack_benders_partial_eta_no_base(joinpath(@__DIR__, "..", "..", ".."))
