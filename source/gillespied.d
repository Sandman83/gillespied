module gillespied; 
import std.math : fabs; 
import std.range : enumerate, empty, front, back, put, ElementType, isInputRange, isForwardRange, assumeSorted;
import std.meta : AliasSeq;
import std.typecons : Flag;
import std.traits : isFloatingPoint, isIterable;  
import std.algorithm.sorting : isSorted;

version(Have_mir_random)
{
    import mir.random : rand, randIndex; 
}
else
{
    import std.random : uniform01, randomSample; 
}

alias possibleTypes = AliasSeq!(float, double, real, size_t);

/**
This auxiliary struct models the time and reaction evolution as stated by Gillespie.
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

private struct GillespieAlgorithm(
    Sandmann sandmann, // whether to sample another uniform value for timings, see [5]
    LDM ldm, // whether to use an array for storing cumulative sum, see [3]
    T // type of propensities
)
{
    static if(isFloatingPoint!T)
    {
        alias FloatingType = T; 
    }
    else
    {
        alias FloatingType = real; 
    }

    static if(ldm == LDM.yes)
    {
        /**
        The only action to perform on creation is to blow up the working cache. Available only in case LDM is chosen
        */
        this(size_t val)
        out
        {
            assert(propsArr_.length); 
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
            return propsArr_.length; 
        }
    }

    /**
	consumes any range of propensities, where cumulativeFold is applicable at. 
	*/
	void put(P)(P props) if(isInputRange!P && is(ElementType!P == T))
    in
    {
        assert(!props.empty);
    }
    do
	{
        import std.algorithm.iteration : cumulativeFold;

        auto resCumFold = props.cumulativeFold!((a, b) => (a + b)); 
 
        static if(ldm == LDM.yes)
        {
            foreach(i, el; resCumFold.enumerate) { propsArr_[i] = el; } 
            
            assert(propsArr_.isSorted);
        }
        else
        {
            static if(isForwardRange!P)
            {
                import std.algorithm.iteration : fold;
                import std.algorithm.comparison : max;

                a0_ = resCumFold.fold!max; 
                nextReaction_ = getNextReaction(resCumFold);
            }
            else
            {
                assert(0); 
            }
        }
    }

    /*
    real tau(T...)(T input) 
        if((sandmann == Sandmann.yes && T.length == 0) || 
           (sandmann == Sandmann.no  && T.length == 1  && is(T[0] == real)))
    in
    {
        import std.math : isNaN; 
        static foreach(i; input) 
        {
            assert(!i.isNaN); 
            assert(i > 0.0); 
            assert(i < 1.0); 
        }
    }
    do
    {
        auto retVal = cast(real)1.0/cast(real)a0; 

        static if(sandmann == Sandmann.no)
        {
            import std.math : log; 
            retVal = - retVal * log(input);
        }

        return retVal; 
    }
    */

	FloatingType tau() //@nogc
	{
        typeof(return) retVal = cast(typeof(return))1.0/cast(typeof(return))a0; 

		static if(sandmann == Sandmann.no)
		{
            import std.math : log; 

            version(Have_mir_random)
            {
                retVal = - retVal * log(rand!(typeof(return)).fabs);
            }
            else
            {
                retVal = - retVal * log(uniform01!(typeof(return)));
            }
        }

        static if(ldm == LDM.yes)
        {
            static if(isFloatingPoint!T)
            {
                propsArr_[] = propsArr_[]/propsArr_.back; 
            }
        }

		return retVal; 
	}

	/**
	This method yields the next reaction index, as stated by Gillespie. This method expects a normed state of the 
	working cache. 
	*/
	size_t index()
	{
        static if(ldm == LDM.yes)
        {
            return getNextReaction(propsArr_); 
        }
        else
        {
            return nextReaction_; 
        }
	}

    private : 
	static if(ldm == LDM.yes)
    {
        T[] propsArr_; 
    }
    else
    {
        T a0_; 
        size_t nextReaction_; 
    }

    auto a0()
    {
        static if(ldm == LDM.yes)
        {
            return propsArr_.back; 
        }
        else
        {
            return a0_;
        }
    }

    size_t getNextReaction(R)(R range)
    {
        static if(isFloatingPoint!T)
        {
            version(Have_mir_random)
            {
                const rndNum = rand!T.fabs; 
            }
            else
            {
                const rndNum = uniform01!T; 
            }
        }
        else
        {
            version(Have_mir_random)
            {
                const rndNum = randIndex!size_t(a0); 
            }
            else
            {
                const rndNum = a0.iota.randomSample(1).front; 
            }
        }
		
		/**
		The usage of persistent working cache allows to accelerate the search of next reaction by binary search.
		See [3] for reference.
		*/
		static if(ldm == LDM.yes)
		{   
            import std.conv : to; 
            assert(range.isSorted, to!string(range));
            return range.assumeSorted.lowerBound(rndNum).length;  // range is propsarr
		}
		else
		{
			/*
			The simple case. Do not assume anything about the range. Count until the value is equal or exceeds the 
			random generated one. 
			*/
            static if(isFloatingPoint!T)
            {
                import std.algorithm.iteration : map; 
                auto workRange = range.map!(el => el/a0()); // range is input propensities
            }
            else
            {
                auto workRange = range; // range is input propensities
            }

            import std.algorithm.searching : countUntil; 
            return workRange.countUntil!"a > b"(rndNum);
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
import std.math : isNaN, abs, approxEqual; 
void testTemplate(Sandmann sandmann, LDM ldm, T, size_t l)()
{
    T[] inputProps; 
    inputProps.length = l; 

    GillespieAlgorithm!(sandmann, ldm, T) s; 

    static if(l)
    {                        
        static if(ldm == LDM.yes)
        {
            s = typeof(s)(inputProps.length); 
        }
        
        static if(isFloatingPoint!T)
        {
            T res;
        }
        else
        {
            real res; 
        }

        res = 0.0; 

        enum runs = 1_000_000; 

        foreach(i; 0 .. runs)
        {
            immutable indexWithZeroPropensity = inputProps.fill; 

            assert(inputProps.sum > 0, ElementType!(typeof(inputProps)).stringof); 
            
            put(s, inputProps); 
            
            immutable tau = s.tau; 
            
            static if(sandmann == Sandmann.yes)
            {
                assert(approxEqual(tau, 1.0/inputProps.sum));
            }
            else
            {
                res += (tau - 1.0/inputProps.sum); 
            }
            assert(s.index != indexWithZeroPropensity); 
        }
        assert(approxEqual(abs(res/runs), 0));
    }
    else
    {
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

size_t fill(R)(R inputProps) 
{
    alias T = ElementType!R;  
    
    import core.stdc.limits : CHAR_BIT;
    
    enum boundExp = CHAR_BIT * size_t.sizeof/2;

    static if(isFloatingPoint!T)
    {
        version(Have_mir_random)
        {
            foreach(ref el; inputProps) { el = rand!T(boundExp).fabs; }
        }
        else
        {
            foreach(ref el; inputProps) { el = uniform01 * (1UL << boundExp); }
        }
        assert(!inputProps.any!isNaN); 
    }
    else
    {
        enum maxEl = T.max/CHAR_BIT; 

        version(Have_mir_random)
        {
            foreach(ref el; inputProps){ el = randIndex!T(maxEl); }
        }
        else
        {
            foreach(ref el; inputProps){ el = maxEl.iota.randomSample(1).front; }
        } 
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
                retVal = retVal.iota.randomSample(1).front;
            }
            inputProps[retVal] = 0;
        }         
    }
    return retVal; 
}