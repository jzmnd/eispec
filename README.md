# eispec

Electrochemical impedance spectroscopy library written in Rust.

Inspired by [LEVM](https://jrossmacdonald.com/levmlevmw/) and other impedance spectroscopy fitting libraries.
This crate provides a way to run complex nonlinear least squares (CNLS) regression on frequency
dependent electrical impedance data with the goal of applying various circuit models to experimental data.

## New types

The following definitions can be found in `eispec::newtypes`
- Representation of electrical current, voltage, impedance and power which are all wrappers around `num::complex::Complex`
- Scalar quantities of frequency and angular frequency

## Electrical circuit components

Found in `eispec::components` and `eispec::circuits`.
These are general purpose building blocks for formulating circuit models.
The current available models are:
- Resistor, Capacitor, Inductor
- Constant Phase Element
- Various dielectric models (Debye, Cole-Cole, Cole-Davidson, Havriliak-Negami)
- Warburg Elements
- Gerischer Element

These can be combined in any way using series and parallel circuits
(`eispec::circuits::SeriesCircuit` and `eispec::circuits::ParallelCircuit`).

Users can also create their own components by implementing the `eispec::components::Component<T>` trait.
This simply requires the definition of the complex impedance as a function of frequency.

## Fitting

CNLS is performed using the [Levenberg–Marquardt algorithm](https://en.wikipedia.org/wiki/Levenberg–Marquardt_algorithm)
applied to the real and imaginary parts of the impedance similar to the approach taken by LEVM.
The code uses a modified version of [rmpfit](https://docs.rs/rmpfit/latest/rmpfit/) which is itself
a direct Rust implementation of the [CMPFIT](https://pages.physics.wisc.edu/~craigm/idl/cmpfit.html) code.

The user can run a fit by simply implementing the `eispec::fit::ImpedanceModel<T>` trait on their
data structure and defining their circuit model.

The data to be fit should in a new type that wraps `eispec::data::ImpedanceData<T>`.
Required data access methods can be applied using the `#[impl_impedance_data_accessors]` macro.

Model parameters are defined using `eispec::fit::ModelParameter` which stores:
- The initial parameter guess
- Whether the parameter is fixed or free in the fitting algorithm
- Any upper or lower limits on the parameter values

For example:

```rust
#[impl_impedance_data_accessors]
struct ImpedanceDataF64Wrapper(ImpedanceData<f64>);

impl ImpedanceModel<f64> for ImpedanceDataF64Wrapper {
    fn model(&self, params: &[f64]) -> Box<(dyn Component<f64>)> {
        \\ circuit model goes here
    }

let mut data = ImpedanceDataF64Wrapper::from_csv("impedance-data-file.csv").unwrap();
data.set_parameters(
    \\ define circuit model parameters and constraints
);

let result = data.fit().unwrap();
```

## Examples

Example fits can be found in the `examples` directory and can be run with e.g.:

```bash
cargo run --example simple_example
```
