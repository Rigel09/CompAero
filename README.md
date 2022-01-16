# CompAero
A python repository for compressible aerodynamics

The goal with this project is to make it easy to determine flow states. This is done by taking compressible aerodynamic tables and computing them. Also more functionality than what is included in the tables is included here as well. This is expected to be a python version of the compressible aerodynamic tables by Virginia tech http://www.dept.aoe.vt.edu/~devenpor/aoe3114/calc.html

Documentation for the project can be found at https://rigel09.github.io/CompAero/



### Install

Currently this packet is <u>not</u> available on PYPI. Install can be accomplished locally from the cloned project by running the following command. 

```python setup.py install```


### Overview 

This repository allows users to easily simulate the Compressible Aerodynamic Tables found in the back of common text books for Compressible Aerodynamics. By knowing the state for any single value or combination of values the user can easily know the rest of them. The repository also prints out the state of the flow in a human readable format. 

Supported Flows

- Isentropic Flows
- Normal Shocks
- Oblique Shocks
- Prandtl Meyer Expansion
- Rayleigh Flow
- Fanno Flow

Upcoming Flows

- Conical Shocks

Examples.

```python
temp = IsentropicRelations(gamma=1.4, a_aStar=1.0235)
print(temp)
```

Produces

```
|====================================================================|
|               Isentropic Flow State at Mach: 1.1749                |
|====================================================================|
|                                                                    |
|γ:---------------------------------------------------------------1.4|
|p0/p:                                                         2.3473|
|T0/T:---------------------------------------------------------1.2761|
|ρ0/ρ:                                                         1.8395|
|A/A*:---------------------------------------------------------1.0235|
|====================================================================|
```

Or:

```python
flow = NormalShockRelations(gamma=1.4, m2=0.5130)
print(flow)
```

Produces

```
|====================================================================|
|               Normal Shock Relations at Mach: 2.4999               |
|====================================================================|
|                                                                    |
|P2/P1:                                                        7.1243|
|ρ2/ρ1:--------------------------------------------------------3.3332|
|T2/T1:                                                        2.1374|
|P02/P01:------------------------------------------------------0.4991|
|P02/P1:                                                       8.5254|
|Dowstream Mach:------------------------------------------------0.513|
|====================================================================|
```

The Rayleigh and Fanno Classes can also simulate Heat and Friction

```python
t1 = 300
p1 = 1
po1 = IsentropicRelations.calc_P0_P(3, 1.4) * p1
flow = FannoFlowRelations(1.4, po_poSt=4.23456790)
flow.apply_pipe_parameters(0.4, 11, 0.005)
print()
print(flow)
print()
print("T2: ", t1 * flow.t2_t1)
print("P2: ", p1 * flow.p2_p1)
print("Po2: ", po1 * flow.po2_po1)
    
flow = RayleighFlowRelations(1.4, mach=1.5)
flow.simulate_heat_addition(1000, 275.2, 287)
print(flow)
```

Produces

```
|====================================================================|
|                  Fanno Relations at Mach: 3.0000                   |
|====================================================================|
|                                                                    |
|γ:                                                               1.4|
|T/T*:---------------------------------------------------------0.4286|
|P/P*:                                                         0.2182|
|ρ/ρ*:---------------------------------------------------------0.5092|
|4FL*/D:                                                       0.5222|
|U/U*:----------------------------------------------------------1.964|
|Flow Type:                                               SUPER_SONIC|
|                                                                    |
|========================= Pipe Parameters ==========================|
|Length For Chocked Flow:                                     10.4432|
|Is Flow Choked? :-----------------------------------------------True|
|Pipe Length:                                                      11|
|Pipe Diameter:---------------------------------------------------0.4|
|Friction Coefficient:                                          0.005|
|                                                                    |
|====================== Down Stream Conditions ======================|
|Mach:                                                              1|
|T/T*:---------------------------------------------------------0.4286|
|P/P*:                                                            1.0|
|P0/P0*:----------------------------------------------------------1.0|
|ρ/ρ*:                                                            1.0|
|4FL*/D:----------------------------------------------------------0.0|
|U/U*:                                                            1.0|
|                                                                    |
|================= Conditions Across Friction Area ==================|
|p2/p1:                                                        4.5826|
|ρ2/ρ1:---------------------------------------------------------1.964|
|T2/T1:                                                        2.3333|
|P02/P01:------------------------------------------------------0.2362|
|4FL*/D2 / 4FL*/D 1:                                              0.0|
|U2/U1:--------------------------------------------------------0.5092|
|====================================================================|


T2:  699.9999999081632
P2:  4.582575694187624
Po2:  8.674491157634302

|====================================================================|
|                 Rayleigh Relations at Mach: 1.5000                 |
|====================================================================|
|                                                                    |
|γ:                                                               1.4|
|T/T*:---------------------------------------------------------0.7525|
|P/P*:                                                         0.5783|
|ρ/ρ*:---------------------------------------------------------0.7685|
|P0/P0*:                                                       1.1215|
|U/U*:---------------------------------------------------------1.3012|
|T0/T0*:                                                       0.9093|
|Flow Type:-----------------------------------------------SUPER_SONIC|
|                                                                    |
|========================= Pipe Parameters ==========================|
|Heat Req. For Chocked Flow:                               27582.0562|
|Is Flow Choked? :----------------------------------------------False|
|Added Heat:                                                     1000|
|Gas Constant R:--------------------------------------------------287|
|Cp:                                                           1004.5|
|T01:-----------------------------------------------------------275.2|
|T02:                                                        276.1955|
|                                                                    |
|====================== Down Stream Conditions ======================|
|Mach:                                                         1.4869|
|T/T*:---------------------------------------------------------0.7525|
|P/P*:                                                          0.586|
|P0/P0*:-------------------------------------------------------1.1152|
|ρ/ρ*:                                                         0.7718|
|T0/T0*:-------------------------------------------------------0.9093|
|U/U*:                                                         1.2957|
|                                                                    |
|================= Conditions Across Heat Addition ==================|
|P2/P1:                                                        1.0133|
|ρ2/ρ1:--------------------------------------------------------1.0043|
|T2/T1:                                                         1.009|
|P02/P01:------------------------------------------------------0.9944|
|T02/T01:                                                      1.0036|
|U2/U1:--------------------------------------------------------0.9958|
|====================================================================|
```



### Class Functions

In some cases all the variables may not be necessary. Each of the classes contain static methods for each of the calculations so that a user can use them without creating  a class instance. In this case the class acts as a collective namespace for functions that make calculations for the same type of flows.

Example.

```python
NormalShockRelations.calc_po2_po1(1.5, 1.4)
IsentropicRelations.calc_P0_P(1.5, 1.4)
```

### Validation

Approximately 170 tests are run on the repository to ensure validation of both class constructors and internal functions for calculations. These tests may not catch all scenarios so please open an issue for any more further bugs. 



### Contributing

Any contributions must include a test and an example. Any bug fixes need a test that fails and then after the bug is fixed the test must pass. 
