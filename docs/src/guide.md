
## Apply a function, then accumulate

```
accum = Mean()
acc = Base.Fix2(accum, abs)

accum = Maximum()
acc = Base.Fix2(accum, x->sqrt(abs(x)))

accum = Mean()
clean(x) = ifelse(ismissing(x), accum.mean, x)
acc = Base.Fix2(accum, clean)

```

