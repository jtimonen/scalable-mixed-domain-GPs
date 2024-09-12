# Scalable mixed-domain Gaussian process modeling and model reduction for longitudinal data

Code for the experiments in the paper "Scalable mixed-domain Gaussian process modeling and model reduction for longitudinal data" by J. Timonen and H. Lähdesmäki.

This repo consists of two parts.

## experiments-01

This part uses `lgpr` and custom Stan code. Contains experiments of *scalability*.
These experiments require `lgpr` version at least 1.1.4 but less than 1.2.0.

```R
install.packages("lgpr")
```

## experiments-02

This part uses [lgpr2](https://github.com/jtimonen/lgpr2) and contains the *model reduction* experiments.
