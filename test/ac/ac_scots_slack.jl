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

if ! isdefined(Main, :TopologyOptimizer) || ! isdefined(Main, :Test)
    include("../common.jl")
end

function test_ac_scots_slack_benders(path_prefix)
    @testset verbose = true "AC-SCOTS (ACRModel)" begin
        run_multi_test(CONSTS[:TEST_SETS][:PIPELINE_TESTSET_SMALL], path_prefix = path_prefix, application = "AC_SCOTS_TIME_SLACK", test_busbars = false, test_reconf = false, test_name = "AC_SCOTS")
    end
end
                    
test_ac_scots_slack_benders(joinpath(@__DIR__, "..", ".."))