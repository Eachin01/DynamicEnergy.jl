using Test
include("simpleHydraulic_dassl.jl")
using ModelingToolkit

@testset "DASSL hydraulic time integration" begin
    mktempdir() do output_dir
        run = HydraulicDASSLDemo.run_hydraulic_dassl(; output_dir, display_plot = false)
        @test ModelingToolkit.SciMLBase.successful_retcode(run.sol)
        @test first(run.sol.t) == 0.0
        @test last(run.sol.t) == 600.0
        @test all(diff(run.sol.t) .> 0)
        @test maximum(diff(run.sol.t)) <= run.settings.dtmax + 1e-10
        @test length(equations(run.sys)) == 2
        @test length(unknowns(run.sys)) == 2
        @test !isempty(observed(run.sys))
        @test count(run.prob.differential_vars) == 1
        @test maximum(abs, run.initial_residual) < 1e-8
        @test run.sol[run.Tank.T][1] ≈ 293.15
        @test run.sol[run.Tank.T][end] ≈ 313.10130954 atol = 1e-4
        @test run.max_temperature_error < 1e-4
        @test all(diff(run.sol[run.Tank.T]) .>= -1e-8)
        @test all(isapprox.(run.sol[run.Pipe.port1.mFlow], 10.0; atol = 1e-7))
        @test all(isapprox.(run.sol[run.Tank.port_2.mFlow], -10.0; atol = 1e-7))
        @test all(isapprox.(run.sol[run.Pump.port2.p], 201325.0; atol = 1e-6))
        @test all(isapprox.(run.sol[run.Tank.port_1.p], 101325.0; atol = 1e-6))
        @test isfile(run.report_path)
        @test occursin("observed(sys)", read(run.report_path, String))
        @test occursin("equations(sys)", read(run.report_path, String))
        @test read(run.plot_path)[1:8] == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]
    end
end
