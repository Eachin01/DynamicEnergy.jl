"""
Library of chilled water models.
"""

using ModelingToolkit: @parameters, @variables, @component, @connector, @named
using ModelingToolkit: System, Equation, connect, instream, Flow, Stream
using ModelingToolkit: t_nounits as t, D_nounits as D
using DataInterpolations: PCHIPInterpolation, ExtrapolationType

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
        mFlow(t), [guess = 0.0, connect = Flow]
        h_outflow(t), [guess = 4186.0 * (273.15 + 20 - ENTHALPY_REF_TEMPERATURE), connect = Stream] # 20C water enthalpy with reference temperature 0C
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
    use_heat_flux::Bool = false
)
    @assert volume > 0 "volume must be positive"
    @assert nPorts >= 1 "nPorts must be at least one"
    vars = @variables begin
        T(t) = T_init
    end
    if use_heat_flux
        @variables flux_input(t)
        append!(vars, flux_input)
        heat_flux = flux_input
    else
        heat_flux = 0
    end
    pars = @parameters begin
    end
    # port = [WaterPort(; name = Symbol("port_$i")) for i in 1:nPorts]
    @named port[1:nPorts] = WaterPort()
    rho = DENSITY
    cp = SPECIFIC_HEAT
    h = temperature2Enthalpy(T)
    eqs = Equation[
        volume * rho * cp * D(T) ~ sum(
            pt.mFlow * ifelse(pt.mFlow > 0, instream(pt.h_outflow), h)
            for pt in port
        ) + ifelse(use_heat_flux, heat_flux, 0),
        sum(pt.mFlow for pt in port) ~ 0,
        [port[1].p ~ port[i].p for i in 2:nPorts]...,
        [port.h_outflow ~ h for port in port]...

    ]
    System(eqs, t, vars, pars; systems = use_heat_flux ? [port...] : port, name)
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
    @named port_1 = WaterPort()
    @named port_2 = WaterPort()
    k = nominalPressureDrop / nominalFlowRate^2
    eqs = [
        port_1.mFlow + port_2.mFlow ~ 0,
        port_1.p - port_2.p ~ k * port_1.mFlow * abs(port_1.mFlow),
        port_1.h_outflow ~ instream(port_2.h_outflow),
        port_2.h_outflow ~ instream(port_1.h_outflow),
    ]
    System(eqs, t, vars, pars; systems = [port_1, port_2], name)
end

@component function Pump_dp(;
    name,
    dp::Float64 = 100000.0,
)
    vars = @variables begin
    end
    pars = @parameters begin
    end
    @named port_1 = WaterPort()
    @named port_2 = WaterPort()
    eqs = [
        port_1.mFlow + port_2.mFlow ~ 0,
        port_1.p - port_2.p ~ - dp,
        port_1.h_outflow ~ instream(port_2.h_outflow),
        port_2.h_outflow ~ instream(port_1.h_outflow),
    ]
    System(eqs, t, vars, pars; systems = [port_1, port_2], name)
end

@component function Pump_curve(;
    name,
    m_flows_data::Vector{Float64},
    dPs_data::Vector{Float64},
    power_data::Vector{Float64})
    @assert length(m_flows_data) == length(dPs_data) == length(power_data) "Three vectors must be with the same length"
    dPcurve = PCHIPInterpolation(dPs_data, m_flows_data, extrapolation = ExtrapolationType.Linear)
    PowerCurve = PCHIPInterpolation(power_data, m_flows_data, extrapolation = ExtrapolationType.Linear)

    vars = @variables begin
        u(t), [guess = 1.0, input = true]
        dP(t), [guess = 0.0, output = true]
        Power(t), [guess = 0.0, output = true]
    end
    pars = @parameters begin
    end
    @named port_1 = WaterPort()
    @named port_2 = WaterPort()
    u_bounded = max(min(u, 1.0), 1e-6)
    m_flow_norm = port_1.mFlow / u_bounded
    eqs = [
        - dP ~ port_1.p - port_2.p,
        dP ~ u_bounded^2 * dPcurve(m_flow_norm),
        Power ~ u_bounded^3 * PowerCurve(m_flow_norm),
        port_1.mFlow + port_2.mFlow ~ 0,
        port_1.h_outflow ~ instream(port_2.h_outflow),
        port_2.h_outflow ~ instream(port_1.h_outflow),
    ]
    System(eqs, t, vars, pars; systems = [port_1, port_2], name)
end

@component function CheckValve_equalPercentage(;
    name,
    nominalFlowRate::Float64 = 1.0,
    nominalPressureDrop::Float64 = 1000.0,
    R::Float64 = 50.0,
    eps_dP::Float64 = 1e-3)
    @assert nominalFlowRate > 0 "nominalFlowRate must be positive"
    @assert nominalPressureDrop > 0 "nominalPressureDrop must be positive"
    vars = @variables begin
        u(t), [guess = 1, input = true]
        dP(t), [guess = 0, output = true]
    end
    pars = @parameters begin
    end
    @named port_1 = WaterPort()
    @named port_2 = WaterPort()
    u_bound = max(min(u, 1.0), 1e-6)
    K_max = nominalPressureDrop / nominalFlowRate^2
    K = K_max * R^(2*(1 - u_bound))
    eqs = [
        dP ~ port_1.p - port_2.p,
        max(dP, 0) ~ K * port_1.mFlow * sqrt(port_1.mFlow^2 + eps_dP^2),
        port_1.mFlow + port_2.mFlow ~ 0,
        port_1.h_outflow ~ instream(port_2.h_outflow),
        port_2.h_outflow ~ instream(port_1.h_outflow),
    ]
    System(eqs, t, vars, pars; systems = [port_1, port_2], name)
end

@component function PolyCOP(;
    name,
    COPatFull::Float64, 
    # COP at full load and nominal evaporator fluid outlet temperature 
    # and condensor fluid inlet temperature
    COPatPeak:: Float64,
    # COP at peak Load ratio
    LRpeakCOP:: Float64,
    # The load ratio where the peak COP occurs
    TevapCOPratio:: Float64,
    # The effect of evaporator's temperature on COP, a number larger than 1
    TcondCOPratio:: Float64,
    # The effect of condensor's temperature on COP, a number larger than 1
    TevapOutNml:: Float64,
    # Nominal evaporator fluid outlet temperature, in K
    TcondInNml:: Float64,
    # Nominal condensor fluid inlet temperature, in K
    COPmin:: Float64,
    # COP min
    COPmax:: Float64,
    # COP max
    )
    @assert 0 < LRpeakCOP < 1 "LRpeakCOP must lie in (0, 1)"
    @assert COPatPeak > COPatFull "COPatPeak must exceed COPatFull, otherwise the COP-LR curve opens upward"
    vars = @variables begin
        LR(t), [input = true]
        TcondIn(t), [input = true]
        TevapOut(t), [input = true]
        COP(t), [output = true]
    end
    pars = @parameters begin
        COPatFull = COPatFull
        COPatPeak = COPatPeak
        LRpeakCOP = LRpeakCOP
        TevapCOPratio = TevapCOPratio
        TcondCOPratio = TcondCOPratio
        TevapOutNml = TevapOutNml
        TcondInNml = TcondInNml
        COPmin = COPmin
        COPmax = COPmax
    end
    a = (COPatPeak - COPatFull)/((1 - LRpeakCOP)^2)
    COP_nml = COPatPeak - a*((LR - LRpeakCOP)^2)
    COP_nml_fixed = COP_nml*
                    ((TevapCOPratio)^(TevapOut - TevapOutNml))*
                    (TcondCOPratio^(TcondInNml - TcondIn))
    eqs = [
        COP ~ max(COPmin, min(COPmax, COP_nml_fixed))
    ]
    System(eqs, t, vars, pars; name)
end

"""
This is chiller cooling output model. 
Equation:
    x = TevapOut - TevapOutNml      
    y = TcondIn - TcondInNml     
    CAPFT = 1 + c1·x + c2·x² + c3·y + c4·y² + c5·x·y
    Q_cap = QClgNominal · max(CAPFT, CAP_min)
Constraint:
    dCAPFT/dx|bound TevapOut = fx ≥ 0   
    dCAPFT/dy|bound TcondIn = fy ≤ 0 
"""
@component function ChillerQcapPoly(;
    name,
    TevapOutNml::Float64 = 7.0 + 273.15,
    TcondInNml:: Float64 = 33.0 + 273.15,
    QclgNml:: Float64 = 10.2e6,
    c:: Vector{Float64} = [-0.00635, 0.008518, -0.015, -0.0002, 0.0],
    TevapOutBound:: Vector{Float64} = [7.0 + 273.15, 12.0 + 273.15],
    TcondInBound:: Vector{Float64} = [33.0 + 273.15, 38.0 + 273.15],
    CAP_min:: Float64 = 0.5,
    eps_mFlow:: Float64 = 1e-6)

    c1, c2, c3, c4, c5 = c
    TevapOutMin, TevapOutMax = TevapOutBound
    TcondInMin, TcondInMax = TcondInBound
    xlo, xhi = TevapOutMin - TevapOutNml, TevapOutMax - TevapOutNml
    ylo, yhi = TcondInMin  - TcondInNml,  TcondInMax  - TcondInNml
    fx_min = c1 + min(2c2*xlo, 2c2*xhi) + min(c5*ylo, c5*yhi)
    fy_max = c3 + max(2c4*ylo, 2c4*yhi) + max(c5*xlo, c5*xhi)
    @assert fx_min > 0 "Not monotonous within the TevapOut range"
    @assert fy_max < 0 "Not monotonous within the TcondIn range"

    vars = @variables begin
        TevapOut(t), [guess = 7.0 + 273.15]
        TcondIn(t), [guess = 33.0 + 273.15]
        Q_cap(t), [output = true]
        Q(t), [output = true]
        LR(t), [output = true]
    end
    pars = @parameters begin
        c1 = c[1]
        c2 = c[2]
        c3 = c[3]
        c4 = c[4]
        c5 = c[5]
        TevapOutNml = TevapOutNml
        TcondInNml = TcondInNml
        QclgNml = QclgNml
        CAP_min = CAP_min
    end
    @named port_1 = WaterPort() # Chilled water In
    @named port_2 = WaterPort() # Chilled water Out
    @named port_3 = WaterPort() # Cooling water In
    @named port_4 = WaterPort() # Cooling water Out
    x = TevapOut - TevapOutNml      
    y = TcondIn - TcondInNml     
    CAPFT = 1 + c1*x + c2*x^2 + c3*y + c4*y^2 + c5*x*y
    h_cwIn = instream(port_1.h_outflow)
    h_clIn = instream(port_3.h_outflow)
    Fcw = port_1.mFlow
    Fcl = port_3.mFlow
    Q_req = max(0, (h_cwIn - temperature2Enthalpy(7.0 + 273.15)) * max(eps_mFlow, Fcw))
    eqs = [
        TevapOut ~ enthalpy2TemperatureK(port_2.h_outflow),
        TcondIn ~ enthalpy2TemperatureK(instream(port_3.h_outflow)),
        Q_cap ~ QclgNml * max(CAPFT, CAP_min),
        Q ~ min(Q_req, Q_cap),
        LR ~ Q/Q_cap,
        port_1.h_outflow ~ instream(port_2.h_outflow),
        port_2.h_outflow ~ h_cwIn - Q/max(eps_mFlow, Fcw),
        port_3.h_outflow ~ instream(port_4.h_outflow),
        port_4.h_outflow ~ h_clIn + Q/max(eps_mFlow, Fcl),
        port_1.mFlow + port_2.mFlow ~ 0,
        port_3.mFlow + port_4.mFlow ~ 0,
    ]
    System(eqs, t, vars, pars; name, systems = [port_1, port_2, port_3, port_4])
end