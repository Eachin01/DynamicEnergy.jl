using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

# Minimal incompressible, constant-property model (SI units).
# mFlow > 0 means flow into a component; T is in kelvin and h = cp*T everywhere.
@connector function WaterPort(;
    name,
    DENSITY::Float64 = 997.0,
    SPECIFIC_HEAT::Float64 = 4186.0,
)
    pars = @parameters begin
        DENSITY = DENSITY
        SPECIFIC_HEAT = SPECIFIC_HEAT
    end
    vars = @variables begin
        p(t), [guess = 101325.0]
        mFlow(t), [guess = 1.0, connect = Flow]
        h_outflow(t), [guess = 4186.0 * 293.15, connect = Stream]
    end
    System(Equation[], t, vars, pars; name)
end

@component function Boundary_pT(;
    name,
    p::Float64 = 101325.0,
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
    T_init::Float64 = 273.15 + 20,
)
    @assert volume > 0 "volume must be positive"
    @assert nPorts >= 1 "nPorts must be at least one"
    ports = [WaterPort(; name = Symbol("port_$i")) for i in 1:nPorts]
    vars = @variables begin
        T(t) = T_init
    end
    rho = ports[1].DENSITY
    cp = ports[1].SPECIFIC_HEAT
    h = cp * T
    eqs = Equation[
        volume * rho * cp * D(T) ~ sum(
            port.mFlow * ifelse(port.mFlow > 0, instream(port.h_outflow), h)
            for port in ports
        ),
        sum(port.mFlow for port in ports) ~ 0,
    ]
    append!(eqs, [ports[1].p ~ ports[i].p for i in 2:nPorts])
    append!(eqs, [port.h_outflow ~ h for port in ports])
    System(eqs, t, vars, []; systems = ports, name)
end

# Memoryless resistance: dp = k*mFlow*abs(mFlow), valid in either flow direction.
@component function PressureDrop(;
    name,
    nominalFlowRate::Float64 = 1.0,
    nominalPressureDrop::Float64 = 1000.0,
)
    @assert nominalFlowRate > 0 "nominalFlowRate must be positive"
    @assert nominalPressureDrop > 0 "nominalPressureDrop must be positive"
    @named port1 = WaterPort()
    @named port2 = WaterPort()
    k = nominalPressureDrop / nominalFlowRate^2
    eqs = [
        port1.mFlow + port2.mFlow ~ 0,
        port1.p - port2.p ~ k * port1.mFlow * abs(port1.mFlow),
        port1.h_outflow ~ instream(port2.h_outflow),
        port2.h_outflow ~ instream(port1.h_outflow),
    ]
    System(eqs, t, [], []; systems = [port1, port2], name)
end

# Positive dp raises pressure from port1 to port2.
# Ideal hydraulic driver only: pump work/heating is not included in the enthalpy balance.
@component function Pump_dp(;
    name,
    dp::Float64 = 100000.0,
)
    @named port1 = WaterPort()
    @named port2 = WaterPort()
    eqs = [
        port1.mFlow + port2.mFlow ~ 0,
        port2.p - port1.p ~ dp,
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

# This topology has no storage/differential state: it is a steady algebraic network.
# Expected: Pipe.port1.mFlow = 10 kg/s, Pump.port2.p = 201325 Pa.
# MixingVolume is fixed above, but deliberately not added to the original topology.
# See simpleHydraulic_dassl.jl for a tank-connected DAE, DASSL integration and plots.
