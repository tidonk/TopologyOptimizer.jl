if ! isdefined(Main, :TopologyOptimizer) || ! isdefined(Main, :Test)
    include("../common.jl")
end

function test_ac_opf_slack_benders(path_prefix)
    @testset verbose = true "AC-OPF (ACRModel)" begin
        run_multi_test(CONSTS[:TEST_SETS][:PIPELINE_TESTSET_SMALL], path_prefix = path_prefix, application = "AC_OPF_TIME_SLACK", test_busbars = false, test_reconf = false, test_name = "AC_OPF")
    end
end
                    
test_ac_opf_slack_benders(joinpath(@__DIR__, "..", ".."))