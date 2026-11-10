using Test

# Keep this example isolated from definitions in simpleHydraulic.jl.
module FixedHydraulicExample
    include("simpleHydraulic_fixed.jl")
end
using .FixedHydraulicExample: WaterPort, Boundary_pT, Pump_dp, PressureDrop, MixingVolume
using ModelingToolkit
using ModelingToolkit: t_nounits as t

# Reuse MTK's installed solver dependency in these checks; the example itself only imports MTK.
const hydraulic_test_algorithm = ModelingToolkit.SimpleNonlinearSolve.SimpleNewtonRaphson()

function check_network(dp)
    @named source = Boundary_pT(T = 313.15)
    @named sink = Boundary_pT(T = 293.15)
    @named pump = Pump_dp(; dp)
    @named pipe = PressureDrop()
    @named network = System([
        connect(source.port, pump.port1),
        connect(pump.port2, pipe.port1),
        connect(pipe.port2, sink.port),
    ], t, [], []; systems = [source, sink, pump, pipe])
    network = mtkcompile(network)
    expected_flow = sign(dp) * sqrt(abs(dp) / 1000.0)
    guess = iszero(dp) ? 0.0 : sign(dp)
    prob = NonlinearProblem(network, [pipe.port1.mFlow => guess])
    sol = solve(prob, hydraulic_test_algorithm)
    @test ModelingToolkit.SciMLBase.successful_retcode(sol)
    @test maximum(abs, prob.f(sol.u, prob.p)) < 1e-6
    @test sol[pipe.port1.mFlow] ≈ expected_flow atol = 1e-8
    @test sol[pump.port2.p] ≈ 101325.0 + dp
    @test sol[source.port.mFlow] ≈ -expected_flow atol = 1e-8
    @test sol[sink.port.mFlow] ≈ expected_flow atol = 1e-8
    if !iszero(dp)
        incoming_h = 4186.0 * (dp > 0 ? 313.15 : 293.15)
        outgoing_h = dp > 0 ? sol[pipe.port2.h_outflow] : sol[pump.port1.h_outflow]
        @test outgoing_h ≈ incoming_h
    end
end

function check_volume(dp)
    @named source = Boundary_pT(T = 313.15)
    @named sink = Boundary_pT(T = 303.15)
    @named pump = Pump_dp(; dp)
    @named pipe = PressureDrop()
    @named tank = MixingVolume(T_init = 293.15)
    @named network = System([
        connect(source.port, pump.port1),
        connect(pump.port2, pipe.port1),
        connect(pipe.port2, tank.port_1),
        connect(tank.port_2, sink.port),
    ], t, [], []; systems = [source, sink, pump, pipe, tank])
    network = mtkcompile(network)
    # Algebraic flow gets a guess, not an extra initial-condition equation.
    flow = sign(dp) * sqrt(abs(dp) / 1000.0)
    prob = ODEProblem(network, [], (0.0, 1.0);
        guesses = [pipe.port1.mFlow => flow], fully_determined = true)
    initial = ModelingToolkit.SciMLBase.build_solution(prob, nothing, [0.0], [prob.u0])
    @test initial[tank.T][1] ≈ 293.15
    @test initial[pipe.port1.mFlow][1] ≈ flow atol = 1e-8
    @test initial[tank.port_1.mFlow][1] + initial[tank.port_2.mFlow][1] ≈ 0 atol = 1e-8
    @test initial[tank.port_1.h_outflow][1] ≈ 4186.0 * 293.15
    @test initial[tank.port_2.h_outflow][1] ≈ 4186.0 * 293.15
    rhs = similar(prob.u0)
    prob.f(rhs, prob.u0, prob.p, 0.0)
    temperature_index = findfirst(isequal(tank.T), unknowns(network))
    incoming_T = dp > 0 ? 313.15 : 303.15
    @test rhs[temperature_index] ≈ abs(flow) / 997.0 * (incoming_T - 293.15) atol = 1e-10
    algebraic_indices = setdiff(eachindex(rhs), [temperature_index])
    @test all(abs.(rhs[algebraic_indices]) .< 1e-6)
end

@testset "Minimal hydraulic correction" begin
    @test length(equations(FixedHydraulicExample.model)) == 1
    @testset "Steady network dp=$dp" for dp in (100000.0, -100000.0, 0.0)
        check_network(dp)
    end
    @testset "Mixing volume dp=$dp" for dp in (100000.0, -100000.0, 0.0)
        check_volume(dp)
    end
    @test_throws AssertionError MixingVolume(name = :invalid_volume, volume = 0.0)
    @test_throws AssertionError MixingVolume(name = :invalid_ports, nPorts = 0)
    @test_throws AssertionError PressureDrop(name = :invalid_flow, nominalFlowRate = 0.0)
end
