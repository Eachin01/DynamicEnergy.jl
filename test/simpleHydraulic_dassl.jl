module HydraulicDASSLDemo

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DASSL: DAEProblem, dassl
using Plots

module Components
    include("simpleHydraulic_fixed.jl")
end
using .Components: Boundary_pT, MixingVolume, Pump_dp, PressureDrop

"""
    run_hydraulic_dassl(; output_dir = @__DIR__, display_plot = true)

Integrate Source -> Pump -> Pipe -> Tank -> Sink with DASSL for 600 seconds.
Print/save MTK equations and observed equations, and save a three-panel PNG.
Run from the repository root:
    julia --project=. test\\simpleHydraulic_dassl.jl

`observed(sys)` contains eliminated variables reconstructed from the solution;
`parameters(sys)` contains fixed symbolic coefficients, not time-series states.
"""
function run_hydraulic_dassl(; output_dir = @__DIR__, display_plot = true)
    settings = (
        inlet_T = 313.15, initial_T = 293.15, volume = 1.0,
        density = 997.0, specific_heat = 4186.0,
        boundary_p = 101325.0, pump_dp = 100000.0,
        nominal_flow = 1.0, nominal_dp = 1000.0,
        tspan = (0.0, 600.0), reltol = 1e-8, abstol = 1e-8, dtmax = 2.0,
    )
    @named Source = Boundary_pT(p = settings.boundary_p, T = settings.inlet_T)
    @named Sink = Boundary_pT(p = settings.boundary_p, T = settings.initial_T)
    @named Pump = Pump_dp(dp = settings.pump_dp)
    @named Pipe = PressureDrop(nominalFlowRate = settings.nominal_flow,
        nominalPressureDrop = settings.nominal_dp)
    @named Tank = MixingVolume(volume = settings.volume, T_init = settings.initial_T)
    connections = [
        connect(Source.port, Pump.port1),
        connect(Pump.port2, Pipe.port1),
        connect(Pipe.port2, Tank.port_1),
        connect(Tank.port_2, Sink.port),
    ]
    @named raw_model = System(connections, t, [], [];
        systems = [Source, Sink, Pump, Pipe, Tank])
    sys = mtkcompile(raw_model)

    expected_flow = settings.nominal_flow * sqrt(settings.pump_dp / settings.nominal_dp)
    tau = settings.volume * settings.density / expected_flow
    analytical_T(time) = settings.inlet_T +
        (settings.initial_T - settings.inlet_T) * exp(-time / tau)

    # Only T is differential. The steady algebraic flow has derivative zero;
    # MTK solves its value and D(T) using the initialization equations.
    prob = DAEProblem(sys, [D(Pipe.port1.mFlow) => 0.0], settings.tspan;
        guesses = [Pipe.port1.mFlow => expected_flow], fully_determined = true)
    # MTK has already computed consistent u0 and du0. CheckInit verifies F(du0,u0)=0
    # without rerunning DASSL's default nonlinear initialization solver selection.
    sol = solve(prob, dassl(); initializealg = ModelingToolkit.SciMLBase.CheckInit(),
        reltol = settings.reltol, abstol = settings.abstol, dtmax = settings.dtmax)
    @assert ModelingToolkit.SciMLBase.successful_retcode(sol) "DASSL integration failed"
    @assert sol.t[end] == settings.tspan[end] "DASSL stopped before the requested end time"
    max_temperature_error = maximum(abs.(sol[Tank.T] .- analytical_T.(sol.t)))

    # DASSL returns adaptive internal steps; use sol.t instead of assuming saveat support.
    temperature_plot = plot(sol.t, sol[Tank.T] .- 273.15;
        label = "Tank (DASSL)", ylabel = "Temperature (degC)", linewidth = 2)
    plot!(temperature_plot, sol.t, analytical_T.(sol.t) .- 273.15;
        label = "Tank (analytical)", linestyle = :dash, linewidth = 2)
    plot!(temperature_plot, sol.t, fill(settings.inlet_T - 273.15, length(sol.t));
        label = "Source", linestyle = :dot)

    flow_plot = plot(sol.t, sol[Pipe.port1.mFlow];
        label = "Pipe inlet", ylabel = "Mass flow (kg/s)", linewidth = 2)
    plot!(flow_plot, sol.t, sol[Tank.port_2.mFlow];
        label = "Tank outlet (negative = leaving)", linestyle = :dash)
    ylims!(flow_plot, (-1.2 * expected_flow, 1.2 * expected_flow))

    pressure_plot = plot(sol.t, sol[Pump.port1.p] ./ 1000;
        label = "Pump inlet", ylabel = "Pressure (kPa)", linewidth = 2)
    plot!(pressure_plot, sol.t, sol[Pump.port2.p] ./ 1000;
        label = "Pump outlet", linewidth = 2)
    plot!(pressure_plot, sol.t, sol[Tank.port_1.p] ./ 1000;
        label = "Tank", linestyle = :dash)
    figure = plot(temperature_plot, flow_plot, pressure_plot;
        layout = (3, 1), size = (1000, 900), xlabel = "Time (s)",
        plot_title = "Hydraulic mixing volume - DASSL", legend = :right)

    initial_residual = similar(prob.u0)
    prob.f(initial_residual, prob.du0, prob.u0, prob.p, first(settings.tspan))
    report = sprint() do io
        println(io, "Julia: ", VERSION, "; ModelingToolkit: ", pkgversion(ModelingToolkit),
            "; DASSL: ", pkgversion(parentmodule(dassl)))
        println(io, "Topology: Source -> Pump -> Pipe -> Tank -> Sink")
        println(io, "Settings: ", settings)
        for (title, expressions) in (
            ("unknowns(sys)", unknowns(sys)),
            ("equations(sys)", equations(sys)),
            ("observed(sys) [observations]", observed(sys)),
            ("parameters(sys)", parameters(sys)),
        )
            println(io, "\n", title, " (", length(expressions), "):")
            for (i, expression) in enumerate(expressions)
                println(io, "  [", i, "] ", expression)
            end
        end
        println(io, "\nu0 (unknown order): ", prob.u0)
        println(io, "du0 (unknown order): ", prob.du0)
        println(io, "differential_vars: ", prob.differential_vars)
        println(io, "Initial residual infinity norm: ", maximum(abs, initial_residual))
        println(io, "Solver: dassl(); retcode: ", sol.retcode)
        println(io, "Time interval: ", (first(sol.t), last(sol.t)),
            "; adaptive output points: ", length(sol.t))
        println(io, "Thermal time constant (s): ", tau)
        println(io, "Final tank temperature (degC): ", sol[Tank.T][end] - 273.15)
        println(io, "Maximum analytical temperature error (K): ", max_temperature_error)
        println(io, "Plot series: Tank.T, analytical T, Source T; Pipe.port1.mFlow,",
            " Tank.port_2.mFlow; Pump.port1.p, Pump.port2.p, Tank.port_1.p")
        println(io, "Temperature plotted in degC, pressure in kPa; equations use SI units.")
    end
    print(report)
    mkpath(output_dir)
    report_path = joinpath(output_dir, "simpleHydraulic_dassl_equations.txt")
    plot_path = joinpath(output_dir, "simpleHydraulic_dassl.png")
    write(report_path, report)
    savefig(figure, plot_path)
    display_plot && display(figure)
    println("Report: ", report_path, "\nPlot: ", plot_path)
    return (; sys, prob, sol, figure, settings, tau, max_temperature_error,
        initial_residual, Source, Sink, Pump, Pipe, Tank, report_path, plot_path)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    hydraulic_run = HydraulicDASSLDemo.run_hydraulic_dassl()
end
