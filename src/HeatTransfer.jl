"""
Library of heat transfer
"""

using ModelingToolkit: @parameters, @variables, @component, @connector, @named
using ModelingToolkit: System, Equation, connect, Flow
using ModelingToolkit: t_nounits as t, D_nounits as D

@connector function HeatPin(;
    name
)
    vars = @variables begin
        T(t), [guess = 273.15 + 20.0]   # Temperature[K]
        q_flow(t), [guess = 0.0, connect = Flow]    # Energy flow rate [W]
    end
    pars = @parameters begin
    end
    System(Equation[], t, vars, pars; name)
end

@component function Resistor(;
    name,
    G_value:: Float64)  # G = UA [W/k]
    vars = @variables begin
    end
    pars = @parameters begin
        G_value = G_value
    end
    @named pin_1 = HeatPin()
    @named pin_2 = HeatPin()
    deltaT = pin_1.T - pin_2.T
    eqs = [
        pin_1.q_flow = G_value*deltaT,
        pin_1.q_flow + pin_2.q_flow ~ 0,
    ]
    System(eqs, t, vars, pars; systems = [pin_1, pin_2], name)
end

@component function Capacitor(;
    name,
    C:: Float64)    # C = mc []
    vars = @variables begin
    end
    pars = @parameters begin
        C = C
    end
    @named pin_1 = HeatPin()
    @named pin_2 = HeatPin()
    deltaT = pin_1.T - pin_2.T
    eqs = [
        D(deltaT) ~ pin_1.q_flow,
        pin_1.q_flow + pin_2.q_flow ~ 0,
    ]
    System(eqs, t, vars, pars; systems = [pin_1, pin_2], name)
end