# Extract raw posterior draws for given Stan parameter names (internal)

Returns a numeric vector (single param) or matrix (multiple params),
flattening the iter x chains dimensions.

## Usage

``` r
.extract_draws_raw(object, params)
```

## Arguments

- object:

  A `bvarnet` object.

- params:

  Character vector of Stan parameter names.

## Value

Numeric vector or matrix.
