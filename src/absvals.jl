AccMinAbs(::Type{T}) = AccumMin(T, abs)
AccMaxAbs(::Type{T}) = AccumMax(T, abs)
AccMeanAbs(::Type{T}) = AccumMean(T, abs)
AccMeanAbs2(::Type{T}) = AccumMean(T, abs2)

AccSumAbs(::Type{T}) = AccumSum(T, abs)
AccProdAbs(::Type{T}) = AccumProd(T, abs)

AccMeanVarAbs(::Type{T}) = AccumMeanVar(T, abs)
AccMeanVarAbs2(::Type{T}) = AccumMeanVar(T, abs2)
