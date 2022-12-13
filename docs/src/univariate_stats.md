# Incremental Statistical Accumulators

### Univariate Statistics

- AccCount
- AccMinimum, AccMaximum, AccExtrema
- AccSum, AccProd
- AccMean, AccGeoMean, AccHarmMean, AccGenMean
- AccMeanVar, AccMeanStd, AccStats
- AccExpWtMean, AccExpWtMeanVar, AccExpWtMeanStd

#### How To Use

```
acc = AccMean(Float32)
for i in eachindex(data_sequence)
    acc( data_sequence[i] )
end
mean = acc()
```

```
acc = AccStats(Float32)
for i in eachindex(data_sequence)
    acc( data_sequence[i] )
end

stats = acc()
# (nobs = 200, mean = -0.09, var = 0.96, std = 0.98, skewness = -0.08, kurtosis = -0.01)
stats.std
# 0.98

### Two-variable Statistics

- AccCov
- AccCor

#### How To Use

- ?

### Multivariate Statistics

- ?

## Incremental Preprocessing

Each Incremental Accumulator is applied to each element of a data sequence `x₁ x₂ .. xₙ` in sequence order.
We support preapplying a function to each xᵢ as it is observed, before the observation is accumulated.
Two examples of use are to force observations to be non-negative, or to replace any `missings` with zeros.



