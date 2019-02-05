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