"""
Library of chilled water models.
"""

using ModelingToolkit: @parameters, @variables, @component, @connector, @named
using ModelingToolkit: System, Equation, connect, instream, Flow, Stream
using ModelingToolkit: t_nounits as t, D_nounits as D

DENSITY::Float64 = 997.0
SPECIFIC_HEAT::Float64 = 4186.0
ENTHALPY_REF_TEMPERATURE::Float64 = 273.15       # Reference temperature for enthalpy calculations (in Kelvin) 

function enthalpy2TemperatureK(h)
    return h / SPECIFIC_HEAT + ENTHALPY_REF_TEMPERATURE
end

function temperature2Enthalpy(T)
    return (T - ENTHALPY_REF_TEMPERATURE) * SPECIFIC_HEAT
end

# Minimal incompressible, constant-property model (SI units).
# mFlow > 0 means flow into a component; T is in kelvin and h = cp*(T - ENTHALPY_REF_TEMPERATURE) everywhere.
@connector function WaterPort(;
    name
)
    vars = @variables begin
        p(t), [guess = 101325.0]
        mFlow(t), [guess = 1.0, connect = Flow]
        h_outflow(t), [guess = 4186.0 * (273.15 + 20 - ENTHALPY_REF_TEMPERATURE), connect = Stream]
    end
    pars = @parameters begin
    end
    System(Equation[], t, vars, pars; 
    assertions = [(p >= 0) => "Pressure must be non-negative, currently $(p) Pa",
                    (h_outflow >= 0) => "Temperature must be non-negative, currently $(enthalpy2TemperatureK(h_outflow) - 273.15) °C"], name)
end

@component function Boundary_pT(;
    name,
    p::Float64 = 101325.0,
    T::Float64 = 273.15 + 20,
)
    vars = @variables begin
    end
    pars = @parameters begin
    end
    cp = SPECIFIC_HEAT
    @named port = WaterPort()
    eqs = [
        port.p ~ p
        port.h_outflow ~ temperature2Enthalpy(T)
    ]
    System(eqs, t, vars, pars; systems = [port], name)
end

@component function MixingVolume(;
    name,
    volume::Float64 = 1.0,
    nPorts::Int = 2,
    T_init::Float64 = 273.15 + 20,
)
    @assert volume > 0 "volume must be positive"
    @assert nPorts >= 1 "nPorts must be at least one"
    vars = @variables begin
        T(t) = T_init
    end
    pars = @parameters begin
    end
    ports = [WaterPort(; name = Symbol("port_$i")) for i in 1:nPorts]
    rho = DENSITY
    cp = SPECIFIC_HEAT
    h = temperature2Enthalpy(T)
    eqs = Equation[
        volume * rho * cp * D(T) ~ sum(
            port.mFlow * ifelse(port.mFlow > 0, instream(port.h_outflow), h)
            for port in ports
        ),
        sum(port.mFlow for port in ports) ~ 0,
    ]
    append!(eqs, [ports[1].p ~ ports[i].p for i in 2:nPorts])
    append!(eqs, [port.h_outflow ~ h for port in ports])
    System(eqs, t, vars, pars; systems = ports, name)
end

# Memoryless resistance: dp = k*mFlow*abs(mFlow), valid in either flow direction.
@component function PressureDrop(;
    name,
    nominalFlowRate::Float64 = 1.0,
    nominalPressureDrop::Float64 = 1000.0,
)
    @assert nominalFlowRate > 0 "nominalFlowRate must be positive"
    @assert nominalPressureDrop > 0 "nominalPressureDrop must be positive"
    vars = @variables begin
    end
    pars = @parameters begin
    end
    @named port1 = WaterPort()
    @named port2 = WaterPort()
    k = nominalPressureDrop / nominalFlowRate^2
    eqs = [
        port1.mFlow + port2.mFlow ~ 0,
        port1.p - port2.p ~ k * port1.mFlow * abs(port1.mFlow),
        port1.h_outflow ~ instream(port2.h_outflow),
        port2.h_outflow ~ instream(port1.h_outflow),
    ]
    System(eqs, t, vars, pars; systems = [port1, port2], name)
end

@component function Pump_dp(;
    name,
    dp::Float64 = 100000.0,
)
    vars = @variables begin
    end
    pars = @parameters begin
    end
    @named port1 = WaterPort()
    @named port2 = WaterPort()
    eqs = [
        port1.mFlow + port2.mFlow ~ 0,
        port2.p - port1.p ~ dp,
        port1.h_outflow ~ instream(port2.h_outflow),
        port2.h_outflow ~ instream(port1.h_outflow),
    ]
    System(eqs, t, vars, pars; systems = [port1, port2], name)
end