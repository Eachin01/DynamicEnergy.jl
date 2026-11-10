using ModelingToolkit, Plots, OrdinaryDiffEq, Setfield
using ModelingToolkit: t_nounits as t, D_nounits as D

@connector function WaterPort(;
    name,
    DENSITY:Float64 = 997.0,
    SPECIFIC_HEAT:Float64 = 4186.0,
)
"""
    Connect for water flow, including:
        1. p: [Pa] pressure - Potential
        2. mFlow: [kg/s] mass flow - Flow
        3. h_outflow: [J/kg] specific enthalpy - Stream
"""
    pars = @parameters begin
        DENSITY = DENSITY
        SPECIFIC_HEAT = SPECIFIC_HEAT
    end
    vars = @variables begin
        p(t), [guess = 0]
        mFlow(t), [guess = 0, connect = Flow]
        h_outflow(t), [guess = 0, connect = Stream]
    end

    System(Equation[], t, vars, pars; name)
end

@component function Boundary_pT(;
    name,
    p::Float64 = 101325,
    T::Float64 = 273.15 + 20,
)
    @named port = WaterPort()
    eqs = [
        port.p ~ p
        port.h_outflow ~ port.SPECIFIC_HEAT * T
    ]
    System(eqs, t, [], []; systems = [port], name)
end

@component function MixingVolume(;
    name,
    volume::Float64 = 1.0,
    nPorts::Int = 2,
    T_init::Float64 = 273.15 + 20
)
    @named ports = [WaterPort(name = "port_$i") for i in 1:nPorts]
    vars = @variables begin
        T(t), [guess = T_init]
    end
    h = cp * (T - 273.15)
    eqs = [
        volume * 997.0 * 4186.0 * D(T) ~ sum(
            port.mFlow > 0 ?
            port.mFlow * instream(port.h_outflow) :
            port.mFlow < 0 ?
            port.mFlow * h : 0
            for port in ports
        ),
        sum(port.mFlow for port in ports) ~ 0,
        ports[1].p ~ ports[n].p for n in 2:nPorts
    ]
    System(eqs, t, vars, []; systems = [ports...], name)
end

@component function PressureDrop(;
    name,
    nominalFlowRate::Float64 = 1.0,
    nominalPressureDrop::Float64 = 1000.0,
)
"""
    Pressure drop across a component, with a first order lag on the pressure drop.
    Dp = k * mFlow^2
"""
    @named port1 = WaterPort()
    @named port2 = WaterPort()
    k = nominalPressureDrop / nominalFlowRate^2

    eqs = [
        port1.mFlow + port2.mFlow ~ 0,

        port1.p - port2.p ~ k * port1.mFlow^2,

        port1.h_outflow ~ instream(port2.h_outflow),
        port2.h_outflow ~ instream(port1.h_outflow),
    ]
    System(eqs, t, [], []; systems = [port1, port2], name)
end

@component function Pump_dp(;
    name,
    dp::Float64 = 100000.0,
)
    @named port1 = WaterPort()
    @named port2 = WaterPort()
    eqs = [
        port1.mFlow + port2.mFlow ~ 0,
        port1.p - port2.p ~ dp,
        port1.h_outflow ~ instream(port2.h_outflow),
        port2.h_outflow ~ instream(port1.h_outflow),
    ]
    System(eqs, t, [], []; systems = [port1, port2], name)
end

@named WaterSource = Boundary_pT()
@named WaterSink = Boundary_pT()
@named Pump = Pump_dp(dp = 100000.0)
@named Pipe = PressureDrop(nominalFlowRate = 1.0, nominalPressureDrop = 1000.0)

hydr_eqs = [
    connect(WaterSource.port, Pump.port1),
    connect(Pump.port2, Pipe.port1),
    connect(Pipe.port2, WaterSink.port),
]

@named model = System(hydr_eqs, t, [], []; systems = [WaterSource, WaterSink, Pump, Pipe])
model = mtkcompile(model)