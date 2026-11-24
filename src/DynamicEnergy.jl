module DynamicEnergy

module HeatTransfer
export HeatPin, Resistor, Capacitor
include("HeatTransfer.jl")
end

module ChilledWater
export WaterPort, Boundary_pT, MixingVolume, PressureDrop, Pump_dp
include("ChilledWater.jl")
end


export deindex, @deindex

function deindex(expr)
    expr isa Expr || return expr
    expr.head == :quote && return expr

    if expr.head == :ref && expr.args[1] isa Expr
        property = expr.args[1]
        if property.head == :. && length(property.args) == 2 &&
                property.args[2] isa QuoteNode && property.args[2].value isa Symbol
            length(expr.args) == 2 || throw(ArgumentError("@deindex expects one index per property"))
            owner = deindex(property.args[1])
            field = property.args[2].value
            index = deindex(expr.args[2])
            if index isa Integer
                return Expr(:., owner, QuoteNode(Symbol(field, "_", string(index))))
            end
            return :($(GlobalRef(Base, :getproperty))($owner,
                $(GlobalRef(Base, :Symbol))($(string(field, "_")),
                    $(GlobalRef(Base, :string))($index))))
        end
    end

    result = Expr(expr.head)
    append!(result.args, map(deindex, expr.args))
    return result
end

macro deindex(expr)
    return esc(deindex(expr))
end

end # module DynamicEnergy