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

## Default TopologyOptimizer Package Test ##
# Execution script to read largest possible contingency sets and then run SCOTS on it
# Execute from within folder `test` using `julia runtests.jl`

include("common.jl")

@testset verbose=true "TopologyOptimizer: Literature Data Test" begin
    @testset verbose=true "AC" begin # ACR formulation
        include("ac/ac_scopf_slack.jl")
        include("ac/ac_opf_slack.jl")
        #include("ac/ac_scots_slack.jl") # IMPLEMENT ME
    end
    @testset verbose=true "DC" begin
        @testset verbose=true "Compact" begin
            include("dc/dc_scots_slack.jl")
            include("dc/dc_scots_time_slack.jl") 
        end

        @testset verbose=true "Benders" begin
            @testset verbose=true "Normal" begin
                @testset verbose=true "BPF in Master" begin
                    include("dc/benders/dc_scots_benders.jl")
                    #include("dc/benders/dc_scots_benders_lp.jl") # FIXME
                    include("dc/benders/dc_scots_benders_leq.jl")
                    include("dc/benders/dc_scots_benders_leq_rf.jl")
                    include("dc/benders/dc_scots_benders_leq_rf_tal2.jl")
                    include("dc/benders/dc_scots_benders_leq_with_warmstart.jl")
                end

                @testset verbose=true "BPF in Subproblem" begin
                    include("dc/benders/dc_scots_benders_bfsp.jl")
                    #include("dc/benders/dc_scots_benders_bfsp_lp.jl") # FIXME
                    include("dc/benders/dc_scots_benders_bfsp_leq.jl")
                    include("dc/benders/dc_scots_benders_bfsp_leq_rf.jl")
                    include("dc/benders/dc_scots_benders_bfsp_leq_with_warmstart.jl")
                end
            end
            
            @testset verbose=true "Nested" begin
                @testset verbose=true "BPF in Intermediate" begin
                    include("dc/benders/dc_scots_benders_nested.jl")
                    #include("dc/benders/dc_scots_benders_nested_lp.jl") # FIXME
                    include("dc/benders/dc_scots_benders_nested_bfip_leq.jl")
                    include("dc/benders/dc_scots_benders_nested_bfip_leq_gmp.jl")
                    include("dc/benders/dc_scots_benders_nested_bfmp_leq_gmp.jl")
                    #include("dc/benders/dc_scots_benders_leq_with_warmstart.jl") # IMPLEMENT ME
                end

                @testset verbose=true "BPF in Subproblem" begin
                    include("dc/benders/dc_scots_benders_nested_bfsp.jl")
                    #include("dc/benders/dc_scots_benders_nested_lp.jl") # FIXME
                    include("dc/benders/dc_scots_benders_nested_bfsp_leq.jl")
                    include("dc/benders/dc_scots_benders_nested_bfsp_leq_gmp.jl")
                    include("dc/benders/dc_scots_benders_nested_bfsp_leq_gmp_csp.jl")
                    #include("dc/benders/dc_scots_benders_leq_with_warmstart.jl") # IMPLEMENT ME
                end
            end
        end
    end
end