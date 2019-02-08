# README #

This repository contains a physical, stochastic time reaction propagation library in D. Based on the 
Gillespie algorithm of 1977 it contains various implementations and enhancements for fast 
reaction index and timing sampling.

The algorithm models the time and reaction evolution as stated by Gillespie, see [1] and [2]. 
After [3], all properly or rather "related to the intended master equation" formulated 
algorithms, modelling physical time propagation algorithms are equivalent. This concerns 
algorithms with rejection as well as the ones which are rejection-free, as the current 
algorithm. The rejection-free algorithms form a subset of the algorithms with rejection. 

The default Gillespie algorithm was enhanced by two features. 
- If the reaction propensities are known, i. e. they do not have to be estimated, one of the 
random numbers needed by 
the original algorithm can be saved. As consequence the longer logarithmic operation is 
cancelled. See [4]. 

- The search of next reaction is done over the cumulative sum of provided propensities. This 
search can be enhanced by using memory space. In this case, any of available search 
algorithms can be applied to the cumulative sum range, which is naturally ordered. In [5] the 
binary search algorithm was applied, whereas in the present case, the search policy is 
managed by the standard library.

Example usage: 

```
import gillespied;

void main()
{
    import std.random : uniform, uniform01, rndGen; 
    import std.math : log, isNaN; 
    import std.range : put; 

	real[] inputPropensities = new real[uniform(1, ubyte.max)];
    foreach(ref el; inputPropensities) 
        el = - uniform01!real.log; 
    auto algorithm = gillespieAlgorithm;
    put(algorithm, inputPropensities); 
    assert(!algorithm.tau.isNaN); 
    assert(algorithm.tau != real.infinity); 
    assert(algorithm.index != inputPropensities.length); 
}
```

```
import gillespied;

void main()
{
    import std.stdio : writeln;
    import std.range : put;
    import std.algorithm.iteration : each;
    import std.algorithm.searching : all;
    import std.math : pow, lrint, isInfinity;
    // stoichiometric matrix for Michaelis-Menten kinetics
    auto VTransposed = [[-1, -1, 1, 0], [1, 1, -1, 0], [0, 1, -1, 1]];
    const nA3root = 84_446_888UL; // see [6]
    // define input params, // see [7]
    const vol = real(5e-13);
    const sInitConc = real(5e-7);
    const eInitConc = real(5e-10);
    const nAV = nA3root * nA3root * vol * nA3root;
    const k1 = real(3000.0); 
    const k2 = real(200.0); 
    const k3 = real(1000.0); 
    const tfinal = 4.2L;
    // declare working stuff
    real[] tVals;
    int[][] xVals;
    auto algorithm = gillespieAlgorithm(); 
    auto a = new real[VTransposed.length];
    // initialize initial values
    tVals ~= 0; 
    xVals ~= [cast(int)((sInitConc * nAV).lrint), cast(int)((eInitConc * nAV).lrint), 0, 0]; 
    // define rates 
    const k = [k1 * nAV/xVals[0][0], k2, k3]; 
    const c = [k[0]/nAV, k[1], k[2]];
    // run the simulation loop
    while(tVals[$-1] < tfinal)
    {
        // calculate reaction propensities from the current system state
        a[0] = c[0]*xVals[$-1][0]*xVals[$-1][1];
        a[1] = c[1]*xVals[$-1][2];
        a[2] = c[2]*xVals[$-1][2];
        // plug in the reaction propensities into the algorithm
        put(algorithm, a); 
        // recieve the timing for the next reaction
        tVals ~= tVals[$-1] + algorithm.tau; 
        // advance the system state container, by the copy of the last system state
        xVals ~= xVals[$-1].dup;
        // if the end of simulation is reached prematurely, exit the simulation loop
        if(tVals[$-1].isInfinity)
        {
            // in this unphysical case, update the time to the final time
            tVals[$-1] = tfinal;
            // ... and exit the simulation immediately.
            break; 
        }
        // update the system state accordingly to the generated reaction index 
        xVals[$-1][] += VTransposed[algorithm.index][]; 
        // all system participants are physical entities, whose amount cannot be negative
        assert(xVals[$-1].all!(el => el >= 0)); 
        // algebraic constraints hold
        assert(xVals[$-1][0] + xVals[$-1][2] + xVals[$-1][3] == xVals[0][0] + xVals[0][2] + xVals[0][3]); 
        assert(xVals[$-1][1] + xVals[$-1][2] == xVals[0][1] + xVals[0][2]); 
    }
    // generate the output
    writeln("time,S,E,ES,P");
    foreach(i, t; tVals)
        writeln(t,",",xVals[i][0],",",xVals[i][1],",",xVals[i][2],",",xVals[i][3]);
}
```

The package provides optional dependency on the random generator of mir.random. This is done 
by the optional feature of dub, so that this random generator take over, if already present 
in the project.

[1] D. T. Gillespie, J. Comput. Phys. 434, 403 (1976).

[2] D. T. Gillespie, 93555, 2340 (1977).

[3] S. A. Serebrinsky, Phys. Rev. E - Stat. Nonlinear, Soft Matter Phys. 83, 2010 (2011).

[4] W. Sandmann, Comput. Biol. Chem. J. 32, 292 (2008).

[5] H. Li and L. R. Petzold, Tech. Rep. 1 (2006). (logarithmic direct method)

[6] Ronald F. Fox and Theodore P. Hill. “A Proposed Exact Integer Value for Avogadro’s Number”

[7] Desmond J. Higham "Modeling and Simulating Chemical Reactions", 2008

Copyright: Copyright (c) 2019- Alexander Orlov. All rights reserved.

License: https://opensource.org/licenses/BSL-1.0, BSL License

Author: Alexander Orlov, sascha.orlov@gmail.com