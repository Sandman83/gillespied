# README #

This repository contains a physical time reaction propagation library in D. Based on the 
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

[1] D. T. Gillespie, J. Comput. Phys. 434, 403 (1976).

[2] D. T. Gillespie, 93555, 2340 (1977).

[3] S. A. Serebrinsky, Phys. Rev. E - Stat. Nonlinear, Soft Matter Phys. 83, 2010 (2011).

[4] W. Sandmann, Comput. Biol. Chem. J. 32, 292 (2008).

[5] H. Li and L. R. Petzold, Tech. Rep. 1 (2006). (logarithmic direct method)

Example usage: 

```
import gillespied;
import std.random : uniform, uniform01, rndGen; 
import std.math : log, isNaN; 
import std.range; 

void main()
{
    import std.stdio; 
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

The package provides optional dependency on the random generator of mir.random. This is done 
by the optional feature of dub, so that this random generator take over, if already present 
in the project.

Copyright: Copyright (c) 2019- Alexander Orlov. All rights reserved.

License: https://opensource.org/licenses/BSL-1.0, BSL License

Author: Alexander Orlov, sascha.orlov@gmail.com