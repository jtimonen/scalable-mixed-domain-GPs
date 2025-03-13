# Scalable mixed-domain Gaussian process modeling and model reduction for longitudinal data

Code for the experiments in the paper "Scalable mixed-domain Gaussian process modeling and model reduction for longitudinal data" by J. Timonen and H. Lähdesmäki.

These experiments require `lgpr` version at least 1.1.4 but less than 1.2.0. This repo consists of two parts.

## experiments-01

Contains Experiments 1-3 of the manuscript. This part uses [lgpr](https://github.com/jtimonen/lgpr) for fitting the exact GP models and custom Stan code for the approximate models. 


## experiments-02

Contains Experiment 4 and additional supplementary experiments. This part uses [lgpr](https://github.com/jtimonen/lgpr) for simulating the fake data and [lgpr2](https://github.com/jtimonen/lgpr2) for model fitting and reduction.
