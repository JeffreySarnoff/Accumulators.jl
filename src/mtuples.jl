

to_tuple(t::Tuple) = t
to_tuple(x::Integer) = (x,)
to_tuple(itr) = tuple(itr...)

function pad(a, b)
    N = max(length(a), length(b))
    Base.fill_to_length(a, 1, Val(N)), Base.fill_to_length(b, 1, Val(N))
end
function pad(a, b, c)
    N = max(max(length(a), length(b)), length(c))
    Base.fill_to_length(a, 1, Val(N)), Base.fill_to_length(b, 1, Val(N)), Base.fill_to_length(c, 1, Val(N))
end



# type-stable padding at end
fill_tolength(t::NTuple{N,Any}, val, ::Val{N}) where {N} = t
fill_tolength(t::Tuple{}, val, ::Val{1}) = (val,)
fill_tolength(t::Tuple{Any}, val, ::Val{2}) = (t..., val)
fill_tolength(t::Tuple{}, val, ::Val{2}) = (val, val)
#function fill_to_length(t::Tuple, val, ::Val{N}) where {N}
#    @inline
#    return (t..., ntuple(i -> val, N - length(t))...)
#end

@inline function fill_tolength(t::Tuple, val, ::Val{_N}) where {_N}
    M = length(t)
    N = _N::Int
    M > N && throw(ArgumentError("input tuple of length $M, requested $N"))
    if @generated
        quote
            (t..., $(fill(:val, (_N::Int) - length(t.parameters))...))
        end
    else
        (t..., fill(val, N - M)...)
    end
end

@inline function fill_tolength(t::Tuple, val, _N::Integer)
    fill_tolength(t::Tuple, val, Val{_N}())
end


# type-stable padding from start
fill_length(t::NTuple{N,Any}, val, ::Val{N}) where {N} = t
fill_length(t::Tuple{}, val, ::Val{1}) = (val,)
fill_length(t::Tuple{Any}, val, ::Val{2}) = (val, t...)
fill_length(t::Tuple{}, val, ::Val{2}) = (val, val)

@inline function fill_length(t::Tuple, val, ::Val{_N}) where {_N}
    M = length(t)
    N = _N::Int
    M > N && throw(ArgumentError("input tuple of length $M, requested $N"))
    if @generated
        quote
            ($(fill(:val, (_N::Int) - length(t.parameters))...), t...)
        end
    else
        (fill(val, N - M)..., t...)
    end
end

@inline function fill_length(t::Tuple, val, _N::Integer)
    fill_length(t::Tuple, val, Val{_N}())
end


# type stable fill from start and into end

fill_upto_length(t::NTuple{N,Any}, val, ::Val{N}, ::Val{B}, ::Val{E}) where {N,B,E} = t
fill_upto_length(t::Tuple{}, val, ::Val{1}, ::Val{B}, ::Val{E}) where {B,E} = (val,)
fill_upto_length(t::Tuple{Any}, val, ::Val{2}, ::Val{B}, ::Val{E}) where {B,E} = (val, t...)
fill_upto_length(t::Tuple{}, val, ::Val{2}, ::Val{B}, ::Val{E}) where {B,E} = (val, val)

@inline function fill_upto_length(t::Tuple, val, ::Val{_N}, ::Val{_B}, ::Val{_E} ) where {_N,_B,_E}
    M = length(t)
    N = _N::Int
    B = _B::Int
    E = _E::Int
    M > N && throw(ArgumentError("input tuple of length $M, requested $N"))
    (N - M) != (B + E) && 
        throw(ArgumentError("padding of (N-M, $(N-M)) != (B+E, $(B+E))"))
    if @generated
        quote
            #($(fill(:val, (B::Int) - length(t.parameters))...), t..., $(fill(:val, (E::Int) - length(t.parameters))...))
            ($(fill(:val, (B::Int) )...), t..., $(fill(:val, (E::Int))...))
        end
    else
        (fill(val, B)..., t...,  fill(val, E)...)
    end
end

@inline function fill_upto_length(t::Tuple, val, _N::Integer, _B::Integer, _E::Integer)
    fill_upto_length(t::Tuple, val, Val{_N}(), Val{_B}(), Val{_E}())
end


# type stable fill1 from start and fill2 into end

fill_upto_length(t::NTuple{N,Any}, valb, vale, ::Val{N}, ::Val{B}, ::Val{E}) where {N,B,E} = t
fill_upto_length(t::Tuple{}, valb, vale, ::Val{1}, ::Val{B}, ::Val{E}) where {B,E} = (val,)
fill_upto_length(t::Tuple{Any}, valb, vale, ::Val{2}, ::Val{B}, ::Val{E}) where {B,E} = (val, t...)
fill_upto_length(t::Tuple{}, valb, vale, ::Val{2}, ::Val{B}, ::Val{E}) where {B,E} = (val, val)

@inline function fill_upto_length(t::Tuple, valb, vale, ::Val{_N}, ::Val{_B}, ::Val{_E} ) where {_N,_B,_E}
    M = length(t)
    N = _N::Int
    B = _B::Int
    E = _E::Int
    M > N && throw(ArgumentError("input tuple of length $M, requested $N"))
    (N - M) != (B + E) && 
        throw(ArgumentError("padding of (N-M, $(N-M)) != (B+E, $(B+E))"))
    if @generated
        quote
            #($(fill(:val, (B::Int) - length(t.parameters))...), t..., $(fill(:val, (E::Int) - length(t.parameters))...))
            ($(fill(:valb, (B::Int) )...), t..., $(fill(:vale, (E::Int))...))
        end
    else
        (fill(valb, B)..., t...,  fill(vale, E)...)
    end
end

@inline function fill_upto_length(t::Tuple, valb, vale, _N::Integer, _B::Integer, _E::Integer)
    fill_upto_length(t::Tuple, valb, vale, Val{_N}(), Val{_B}(), Val{_E}())
end




# type-stable padding from start
fill_length(t::NTuple{N,Any}, val, ::Val{N}) where {N} = t
fill_length(t::Tuple{}, val, ::Val{1}) = (val,)
fill_length(t::Tuple{Any}, val, ::Val{2}) = (val, t...)
fill_length(t::Tuple{}, val, ::Val{2}) = (val, val)

@inline function fill_length(t::Tuple, val, ::Val{_N}) where {_N}
    M = length(t)
    N = _N::Int
    M > N && throw(ArgumentError("input tuple of length $M, requested $N"))
    if @generated
        quote
            ($(fill(:val, (_N::Int) - length(t.parameters))...), t...)
        end
    else
        (fill(val, N - M)..., t...)
    end
end

@inline function fill_length(t::Tuple, val, _N::Integer)
    fill_length(t::Tuple, val, Val{_N}())
end


# type-stable padding at end
fill_tolength(t::NTuple{N,Any}, val, ::Val{N}) where {N} = t
fill_tolength(t::Tuple{}, val, ::Val{1}) = (val,)
fill_tolength(t::Tuple{Any}, val, ::Val{2}) = (t..., val)
fill_tolength(t::Tuple{}, val, ::Val{2}) = (val, val)
#function fill_to_length(t::Tuple, val, ::Val{N}) where {N}
#    @inline
#    return (t..., ntuple(i -> val, N - length(t))...)
#end

@inline function fill_tolength(t::Tuple, val, ::Val{_N}) where {_N}
    M = length(t)
    N = _N::Int
    M > N && throw(ArgumentError("input tuple of length $M, requested $N"))
    if @generated
        quote
            (t..., $(fill(:val, (_N::Int) - length(t.parameters))...))
        end
    else
        (t..., fill(val, N - M)...)
    end
end

@inline function fill_tolength(t::Tuple, val, _N::Integer)
    fill_tolength(t::Tuple, val, Val{_N}())
end









# constructing from an iterator

# only define these in Base, to avoid overwriting the constructors
# NOTE: this means this constructor must be avoided in Core.Compiler!
if nameof(@__MODULE__) === :Base

    function tuple_type_tail(T::Type)
        @_foldable_meta # TODO: this method is wrong (and not :foldable)
        if isa(T, UnionAll)
            return UnionAll(T.var, tuple_type_tail(T.body))
        elseif isa(T, Union)
            return Union{tuple_type_tail(T.a),tuple_type_tail(T.b)}
        else
            T.name === Tuple.name || throw(MethodError(tuple_type_tail, (T,)))
            if isvatuple(T) && length(T.parameters) == 1
                va = unwrap_unionall(T.parameters[1])::Core.TypeofVararg
                (isdefined(va, :N) && isa(va.N, Int)) || return T
                return Tuple{Vararg{va.T,va.N - 1}}
            end
            return Tuple{argtail(T.parameters...)...}
        end
    end

    (::Type{T})(x::Tuple) where {T<:Tuple} = convert(T, x)  # still use `convert` for tuples

    Tuple(x::Ref) = tuple(getindex(x))  # faster than iterator for one element
    Tuple(x::Array{T,0}) where {T} = tuple(getindex(x))

    (::Type{T})(itr) where {T<:Tuple} = _totuple(T, itr)

    _totuple(::Type{Tuple{}}, itr, s...) = ()

    function _totuple_err(@nospecialize T)
        @noinline
        throw(ArgumentError("too few elements for tuple type $T"))
    end

    function _totuple(::Type{T}, itr, s::Vararg{Any,N}) where {T,N}
        @inline
        y = iterate(itr, s...)
        y === nothing && _totuple_err(T)
        t1 = convert(fieldtype(T, 1), y[1])
        # inference may give up in recursive calls, so annotate here to force accurate return type to be propagated
        rT = tuple_type_tail(T)
        ts = _totuple(rT, itr, y[2])::rT
        return (t1, ts...)::T
    end

    # use iterative algorithm for long tuples
    function _totuple(T::Type{All32{E,N}}, itr) where {E,N}
        len = N + 32
        elts = collect(E, Iterators.take(itr, len))
        if length(elts) != len
            _totuple_err(T)
        end
        (elts...,)
    end

    _totuple(::Type{Tuple{Vararg{E}}}, itr, s...) where {E} = (collect(E, Iterators.rest(itr, s...))...,)

    _totuple(::Type{Tuple}, itr, s...) = (collect(Iterators.rest(itr, s...))...,)

    # for types that `apply` knows about, just splatting is faster than collecting first
    _totuple(::Type{Tuple}, itr::Array) = (itr...,)
    _totuple(::Type{Tuple}, itr::SimpleVector) = (itr...,)
    _totuple(::Type{Tuple}, itr::NamedTuple) = (itr...,)
    _totuple(::Type{Tuple}, x::Number) = (x,) # to make Tuple(x) inferable

end


# `ntuple`, for constructing tuples of a given length

"""
    ntuple(f::Function, n::Integer)
Create a tuple of length `n`, computing each element as `f(i)`,
where `i` is the index of the element.
# Examples
```jldoctest
julia> ntuple(i -> 2*i, 4)
(2, 4, 6, 8)
```
"""
@inline function ntuple(f::F, n::Integer) where {F}
    # marked inline since this benefits from constant propagation of `n`
    t = n == 0 ? () :
        n == 1 ? (f(1),) :
        n == 2 ? (f(1), f(2)) :
        n == 3 ? (f(1), f(2), f(3)) :
        n == 4 ? (f(1), f(2), f(3), f(4)) :
        n == 5 ? (f(1), f(2), f(3), f(4), f(5)) :
        n == 6 ? (f(1), f(2), f(3), f(4), f(5), f(6)) :
        n == 7 ? (f(1), f(2), f(3), f(4), f(5), f(6), f(7)) :
        n == 8 ? (f(1), f(2), f(3), f(4), f(5), f(6), f(7), f(8)) :
        n == 9 ? (f(1), f(2), f(3), f(4), f(5), f(6), f(7), f(8), f(9)) :
        n == 10 ? (f(1), f(2), f(3), f(4), f(5), f(6), f(7), f(8), f(9), f(10)) :
        _ntuple(f, n)
    return t
end

function _ntuple(f::F, n) where {F}
    @noinline
    (n >= 0) || throw(ArgumentError(string("tuple length should be ≥ 0, got ", n)))
    ([f(i) for i = 1:n]...,)
end

function ntupleany(f, n)
    @noinline
    (n >= 0) || throw(ArgumentError(string("tuple length should be ≥ 0, got ", n)))
    (Any[f(i) for i = 1:n]...,)
end

# inferable ntuple (enough for bootstrapping)
ntuple(f, ::Val{0}) = ()
ntuple(f, ::Val{1}) = (@inline; (f(1),))
ntuple(f, ::Val{2}) = (@inline; (f(1), f(2)))
ntuple(f, ::Val{3}) = (@inline; (f(1), f(2), f(3)))

"""
    ntuple(f, ::Val{N})
Create a tuple of length `N`, computing each element as `f(i)`,
where `i` is the index of the element. By taking a `Val(N)`
argument, it is possible that this version of ntuple may
generate more efficient code than the version taking the
length as an integer. But `ntuple(f, N)` is preferable to
`ntuple(f, Val(N))` in cases where `N` cannot be determined
at compile time.
# Examples
```jldoctest
julia> ntuple(i -> 2*i, Val(4))
(2, 4, 6, 8)
```
"""
@inline function ntuple(f::F, ::Val{N}) where {F,N}
    N::Int
    (N >= 0) || throw(ArgumentError(string("tuple length should be ≥ 0, got ", N)))
    if @generated
        :(@ntuple $N i -> f(i))
    else
        Tuple(f(i) for i = 1:N)
    end
end

@inline function fill_to_length(t::Tuple, val, ::Val{_N}) where {_N}
    M = length(t)
    N = _N::Int
    M > N && throw(ArgumentError("input tuple of length $M, requested $N"))
    if @generated
        quote
            (t..., $(fill(:val, (_N::Int) - length(t.parameters))...))
        end
    else
        (t..., fill(val, N - M)...)
    end
end
