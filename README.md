# CompAero
A python repository for compressible aerodynamics

The goal with this project is to make it easy to determine flow states. This is done by taking compressible aerodynamic tables and computing them. Also more functionality than what is included in the tables is included here as well. This is expected to be a python version of the compressible aerodynamic tables by Virginia tech http://www.dept.aoe.vt.edu/~devenpor/aoe3114/calc.html



This repository allows users to easily simulate the Compressible Aerodynamic Tables found in the back of common text books for Compressible Aerodynamics. By knowing the state for any single value or combination of values the user can easily know the rest of them. The repository also prints out the state of the flow in a human readable format. 

Examples.

```python
temp = IsentropicRelations(gamma=1.4, a_aStar=1.0235)
print(temp)
```

Produces

```
Isentropic Flow State at Mach: 1.1749  
---------------------------------------
p0/p:  2.3473
T0/T:  1.2761
ρ_0/ρ: 1.8395
A/A*:  1.0235
```

Or:

```python
flow = NormalShockRelations(gamma=1.4, m2=0.5130)
print(flow)
```

Produces

```
|=====Normal Shock Relations at Mach: 2.4999=====|
|                                                |
|p2/p1:                                    7.1243|
|ρ2/ρ1:------------------------------------3.3332|
|T2/T1:                                    2.1374|
|po2/po1:----------------------------------0.4991|
|po2/p1:                                   8.5254|
|Downstream Mach:                           0.513|
|================================================|
```

### Class Functions

In some cases all the variables may not be necessary. Each of the classes contain static methods for each of the calculations so that a user can use them without creating  a class instance. 

Example.

```python
NormalShockRelations.calc_po2_po1(1.5, 1.4)
IsentropicRelations.calc_p0_p(1.5, 1.4)
```

### Validation

Approximately 170 tests are run on the repository to ensure validation of both class constructors and internal functions for calculations. These tests may not catch all scenarios so please open an issue for any more further bugs. 



### Contributing

Any contributions must include a test and an example. Any bug fixes need a test that fails and then after the bug is fixed the test must pass. 
