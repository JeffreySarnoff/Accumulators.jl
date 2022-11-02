AccMinAbs(::Type{T}) where {T} = AccumMin(T, abs)
AccMinAbs2(::Type{T}) where {T} = AccumMin(T, abs2)

AccMaxAbs(::Type{T}) where {T} = AccumMax(T, abs)
AccMaxAbs2(::Type{T}) where {T} = AccumMax(T, abs2)

AccMeanAbs(::Type{T}) where {T} = AccumMean(T, abs)
AccMeanAbs2(::Type{T}) where {T} = AccumMean(T, abs2)

AccSumAbs(::Type{T}) where {T} = AccumSum(T, abs)
AccSumAbs2(::Type{T}) where {T} = AccumSum(T, abs2)

AccProdAbs(::Type{T}) where {T} = AccumProd(T, abs)
AccProdAbs2(::Type{T}) where {T} = AccumProd(T, abs2)

AccMeanVarAbs(::Type{T}) where {T} = AccumMeanVar(T, abs)
AccMeanVarAbs2(::Type{T}) where {T} = AccumMeanVar(T, abs2)
