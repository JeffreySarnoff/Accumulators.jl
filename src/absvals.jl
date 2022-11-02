AccMinAbs(::Type{T}) where {T} = AccumMin(T, abs)
AccMaxAbs(::Type{T}) where {T} = AccumMax(T, abs)
AccMeanAbs(::Type{T}) where {T} = AccumMean(T, abs)
AccMeanAbs2(::Type{T}) where {T} = AccumMean(T, abs2)

AccSumAbs(::Type{T}) where {T} = AccumSum(T, abs)
AccProdAbs(::Type{T}) where {T} = AccumProd(T, abs)

AccMeanVarAbs(::Type{T}) where {T} = AccumMeanVar(T, abs)
AccMeanVarAbs2(::Type{T}) where {T} = AccumMeanVar(T, abs2)
