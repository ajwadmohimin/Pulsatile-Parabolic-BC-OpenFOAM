# Pulsatile-Parabolic-BC-OpenFOAM
A boundary condition for OpenFOAM to produce a parabolic velocity profile from a time-varying flow rate. It can be useful when the inlet flow rate is period (i.e., in Cardiovascular Flow Simulation).

## Installation / Compiling
* Clone the repository and place the myQFourierFunc folder in your desired location.
* Open `terminal` in the myQFourierFunc folder.
* Run `wmake`
* The boundary condition is compiled in the `$FOAM_USER_LIBBIN`.

## Usage
```
<patchName>
    {
        type       myQFourierFunc;
        Q       ( (3 0) (0 2) (2 4) ); // list of complex values
        omega   7.854;
        value      uniform (0 0 0); // optional dummy value
    }
```
* Don't forget to modify `controlDict` by adding the library.
* Here, Q is a list of Fourier Coefficients (a + ib form). a and b is found by decomposing the flow curve into a Fourier series.
* omega can be found from $2 \pi T$, where, $T$ is the time-period of the pulsatile flow curve.

N.B.: If you find it useful, your star will be appreciated!
