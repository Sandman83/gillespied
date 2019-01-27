/**
Copyright: Copyright (c) 2019- Alexander Orlov. All rights reserved.
License: $(LINK2 https://opensource.org/licenses/BSL-1.0, BSL-1.0 License).
Authors: Alexander Orlov, $(LINK2 mailto:sascha.orlov@gmail.com, sascha.orlov@gmail.com) 
*/

/**
The module contains various physical time reaction propagation algorithms based on the Gillespie algorithm from 1977. 
In the current version this module contains
- the default version of Gillespie
- The modification introduced by Sandmann 2008, which makes possible faster sample reaction arrival times,
    in case all propensities of the reactions are known exactly. I. e. this modification can be applied, in case 
    propensities were not estimated. 
- The modification introduced by Petzold et. al. 2006. At the expense of allocating an ancillary propensities array, 
    where their cumulative sum is stored, the velocity of sampling the reaction index increases, dependent on the 
    chosen search policy of assumedSorted ranges.
- The module provides also a common testing routine, containing tests for all combinations of modifications. 
*/
module gillespied; 
import std.math : fabs, isNaN; 
import std.range : enumerate, ElementType, assumeSorted, isForwardRange, isInputRange, empty, front, back, put, iota;
import std.meta : AliasSeq;
import std.typecons : Flag;
import std.traits : isFloatingPoint, isIterable;  
import std.algorithm.sorting : isSorted;

/**
Having mir library imported, this one will prefere it over the standard library random generator, assuming its 
performance is better.
*/
version(Have_mir_random)
{
    import mir.random : rand, randIndex; 
}
else
{
    import std.random : uniform01, uniform; 
}

alias possibleTypes = AliasSeq!(float, double, real, size_t);
enum minTestedSizes = 0; 
enum maxTestedSizes = 8; // tested up to 24

/**
This struct models the time and reaction evolution as stated by Gillespie.
An instance of the struct is initiated by exposing reaction propensities to the constructor. Once an instance is
created, the opCall, which has the same interface as the constructor take over. 
Passing the reaction propensities (which is in general a range) updates the intern array to contain a cumulative sum of 
the rates. 
The array is persistent and therefore must be allocated only once. Furthermore, due to persistence an optimization via
binary search of the original algorithm can be achieved, see [3]. 
[1] D. T. Gillespie, J. Comput. Phys. 434, 403 (1976).
[2] D. T. Gillespie, 93555, 2340 (1977).
[3] H. Li and L. R. Petzold, Tech. Rep. 1 (2006). (logarithmic direct method)

// not used any more for defining the type:
It is forbidden to instanciate a propagator without knowledge of all existent reaction rates. This is due 
[4] S. A. Serebrinsky, Phys. Rev. E - Stat. Nonlinear, Soft Matter Phys. 83, 2010 (2011).

Another optimization is given by [5]. I. e. if the upper bound of reaction rates is found (last member of the 
cumulative sum working cache), the time step can be taken as 1/a0, without another random number generation and 
logarithmic operation. 
[5] W. Sandmann, Comput. Biol. Chem. J. 32, 292 (2008).
*/

alias Sandmann = Flag!"Sandmann"; 
alias LDM = Flag!"LDM"; 

/**
The floating type in use is defined either by the user on algorithm declaration or is chosen to be real, in case 
propensities come in form of integer values. 
*/
template FloatingType(T)
{
    static if(isFloatingPoint!T)
    {
        alias FloatingType = T; 
    }
    else
    {
        alias FloatingType = real; 
    }
}

private struct GillespieAlgorithm(
    Sandmann sandmann, // whether to sample another uniform value for timings, see [5]
    LDM ldm, // whether to use an array for storing cumulative sum, see [3]
    T // type of propensities
)
{
    static if(ldm == LDM.yes)
    {
        /**
        The only action to perform on creation is to blow up the working cache. Available only in case LDM is chosen
        */
        this(size_t val)
        out
        {
            assert(propensities.length); 
        }
        do
        {
            propsArr_.length = val; 
        }

        /**
        In LDM case, it is naturally to make the length of the object available.
        */
        size_t length() @nogc
        {
            return propensities.length; 
        }
    }

    /**
	consumes any range of propensities, where cumulativeFold is applicable at. 
	*/
	void put(P)(P props) if(isInputRange!P && is(ElementType!P == T))
    in
    {
        assert(!props.empty);
        static if(isFloatingPoint!T)
        {
            assert(!props.any!isNaN); 
        }
    }
    do
	{
        import std.algorithm.iteration : cumulativeFold;

        auto resCumFold = props.cumulativeFold!((a, b) => [a, b].sum); 
 
        static if(ldm == LDM.yes)
        {
            foreach(i, el; resCumFold.enumerate) 
            {
                propensities[i] = el;
            }
            
            assert(propensities.isSorted);
        }
        else
        {
            static assert(isForwardRange!P);
            import std.algorithm.iteration : fold;
            import std.algorithm.comparison : max;

            a0_ = resCumFold.fold!max; 
            nextReaction_ = getNextReaction(resCumFold);
        }
    }

    /**
    This method provides the arrival time to the next reaction. 
    Returns: the according floating point type value, corresponding to the provided propensities. This value is the 
    exponential distributed time needed for the next reaction to arrive. 
    */
	FloatingType!T tau() @nogc
	{
        typeof(return) retVal = cast(typeof(return))1.0/cast(typeof(return))a0; 

		static if(sandmann == Sandmann.no)
		{
            import std.math : log; 
            
            FloatingType!T rndNum = 0.0; 

            version(Have_mir_random)
            {
                rndNum = rand!(typeof(return)).fabs;
            }
            else
            {
                rndNum = uniform01!(typeof(return));
            }
            retVal = - retVal * log(rndNum); 
        }

		return retVal; 
	}

	/**
	This method yields the next reaction index, as stated by Gillespie. The sampling is done according to the original 
    sampling algorithm, taking advantage of stored propensities array if applicable.
	*/
	size_t index()() //not possible to mark @nogc directly, due to phobos uniform is not @nogc
	{
        static if(ldm == LDM.yes)
        {
            return getNextReaction(propensities); 
        }
        else
        {
            return nextReaction_; 
        }
	}

    private: 
	static if(ldm == LDM.yes)
    {
        T[] propsArr_; 
        auto propensities() @nogc 
        {
            return propsArr_;
        }
    }
    else
    {
        T a0_; 
        size_t nextReaction_; 
    }

    auto a0() @nogc
    {
        static if(ldm == LDM.yes)
        {
            return propensities.back; 
        }
        else
        {
            return a0_;
        }
    }

    size_t getNextReaction(R)(R range) 
    {
        if(a0 > 0)
        {
            T rndNum;

            static if(isFloatingPoint!T)
            {
                version(Have_mir_random)
                {
                    rndNum = rand!T.fabs; 
                }
                else
                {
                    rndNum = uniform01!T; 
                }
                rndNum *= a0; 
            }
            else
            {
                version(Have_mir_random)
                {
                    rndNum = randIndex!size_t(a0); 
                }
                else
                {
                    rndNum = uniform(T(0), a0); 
                }
            }
            
            /**
            The usage of persistent working cache allows to accelerate the search of next reaction by binary search.
            See [3] for reference.
            */
            static if(ldm == LDM.yes)
            {
                assert(range.isSorted);
                return range.assumeSorted.lowerBound(rndNum).length;  // range is propsarr
            }
            else
            {
                /*
                The simple case. Do not assume anything about the range. Count until the value is equal or exceeds the 
                random generated one. 
                */
                import std.algorithm.searching : countUntil; 
                return range.countUntil!"a > b"(rndNum);
            }
        }
        else
        {
            return typeof(return).max; 
        }
    }
}

/**
The convenience function to release a gillespie algorithm structure with some preset properties. 
In this case, the preset handles the memory allocation due using an internal buffer for storing passed propensities. 
Other properties do not affect memory allocations. 
*/
auto gillespieAlgorithm(
    Sandmann sandmann = Sandmann.yes,
    LDM ldm = LDM.no,
    T = real)(size_t val = 0)
{
    static if(ldm == LDM.yes)
    {
        return GillespieAlgorithm!(sandmann, ldm, T)(val);
    }
    else
    {
        return GillespieAlgorithm!(sandmann, ldm, T)();
    }
}

version(unittest):
import std.algorithm.iteration : sum;
import std.algorithm.searching : any;
import std.math : isInfinity, abs, approxEqual; 
/**
Common testing function.
*/
void testTemplate(Sandmann sandmann, LDM ldm, T, size_t l)()
{
    // 1. generate the environmental propensities of needed length
    T[] inputProps = new T[l]; 

    // 2a. Declare the usage of the algorithm
    GillespieAlgorithm!(sandmann, ldm, T) s; 

    // 3. In case the environment does not provide propensities, it is an error to use the algorithm.
    static if(l)
    {
        // 2b. If using Petzold's enhancement, the algorithm has not only to be declared but also to be allocated
        static if(ldm == LDM.yes)
        {
            s = typeof(s)(inputProps.length); 
        }

        // 4. The algorithm is intrinsically stochastic. Choosing a large amount of runs to be able to take an average. 
        enum runs = 1_000_000; 

        // 5. Declaring and initializing the checked result. 
        FloatingType!T res = 0.0; 
        FloatingType!T[] resArray = new FloatingType!T[runs]; 
        resArray[] = 0.0; 
        
        foreach(i; 0 .. runs)
        {
            // 7a. Fill the mocking propensities with random values. Sometimes, there will be a zero propensity. The 
            // index of the zero propensity has to be tracked to be checked later for not being choosen in all 
            // circumstances. 
            immutable indexWithZeroPropensity = inputProps.fill; 
            
            // 7b. Sometimes, there are cases, when propensities become zero altogether
            if((indexWithZeroPropensity >= inputProps.length) && !i)
            {
                inputProps[] = 0; 
            }
            
            // 8. Feed the propensities to the algorithm
            put(s, inputProps); 
            
            // 9. Calculate the physical arrival time for the next reaction. This is simulated by exp distribution with 
            // the single parameter of the intensity taken as the sum of all propensities feeded in the step 8. 
            immutable tau = s.tau; 

            if(tau.isInfinity)
            {
                // 10a. In case all propensites were zero, arriving time is infinit. In this case no valid reaction 
                // index can be generated. 
                immutable ind = s.index; 
                assert(ind >= inputProps.length);
            }
            else
            {
                // 10b. Otherwise, different assertions on the arrival timing can be done: 
                static if(sandmann == Sandmann.yes)
                {
                    // 11a. If the Sandmann modification is taken into account, the timing was not sampled, but taken 
                    // exactly. Assert on its value directly.
                    assert(approxEqual(tau, 1.0/inputProps.sum));
                }
                else
                {
                    // 11b. Otherwise, the mean of the arrival timings have to be inspected. Form a sum of differences 
                    // of the sampled timings and the expected ones. 
                    res = [res, (tau - 1.0/inputProps.sum)].sum; 
                    resArray[i] = tau - 1.0/inputProps.sum; 
                }

                // 12. In any case, the sampled reaction index cannot be the excluded one in step 7a, as the propensity
                // of the according reaction was set to zero. 
                assert(s.index != indexWithZeroPropensity);
            }
        }
        // 11c. The mean of the differences has to be approximately zero, as taking the average on all runs has to obey
        // the expectation of the exponential distribution.
        import std.conv : to;  
        assert(!resArray.any!isNaN);
        assert(approxEqual(resArray.sum/runs, res/runs));
        assert(approxEqual(abs(resArray.sum/runs), 0), abs(resArray.sum/runs).to!string);
        assert(approxEqual(abs(res/runs), 0), abs(res/runs).to!string);
    }
    else
    {
        // See 3.
        import std.exception : assertThrown; 
        static if(ldm == LDM.yes)
        {
            assertThrown!Error(typeof(s)(inputProps.length));
        }
        else
        {
            assertThrown!Error(put(s, inputProps));
        }
    }
}

// Auxillary function for testing. Mocks the environmental propensities. 
size_t fill(R)(R inputProps)
{
    alias T = ElementType!R;  
    
    foreach(ref el; inputProps)
    {
        T rndNum = 0; 
        
        auto maxEl = T.max/inputProps.length; 

        static if(isFloatingPoint!T)
        {
            version(Have_mir_random)
            {
                rndNum = rand!T.fabs * maxEl;
            }
            else
            {
                rndNum = uniform01!T * maxEl;
            }
        }
        else
        {
            version(Have_mir_random)
            {
                rndNum = randIndex!T(maxEl);
            }
            else
            {
                rndNum = uniform(T(0), maxEl);
            }
        }
        el = rndNum;
    }

    size_t retVal = inputProps.length;
    
    if(retVal > 1)
    {
        bool modify;

        version(Have_mir_random)
        {
            modify = rand!bool;
        }
        else
        {
            modify = uniform01!real <= 0.5;
        }

        if(modify)
        {
            version(Have_mir_random)
            {
                retVal = randIndex!size_t(retVal);
            }
            else
            {
                retVal = uniform(size_t(0), retVal);
            }
            inputProps[retVal] = 0;
        }
    }
    return retVal;
}