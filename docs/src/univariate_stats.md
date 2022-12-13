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
acc()
# (mean = 6.8f0, var = 42.7f0, std = 6.534524f0, skewness = 1.2196333f0, kurtosis = -0.1658349f0)

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



