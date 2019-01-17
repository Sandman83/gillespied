module gillespied; 
import std.math : isNaN, fabs, log, abs, approxEqual; 
import std.algorithm.comparison : max;
import std.range;
import std.meta : AliasSeq;
import std.typecons : Flag;
import std.traits : EnumMembers, isFloatingPoint;  

debug import std.stdio; 

debug = bug; 

debug(bug)
{
    private alias possibleTypes = AliasSeq!(real);
    version = mirRandom;
}
else
{
    private alias possibleTypes = AliasSeq!(float, double, real, size_t);
    version = stdRandom;
}

version(mirRandom)
{
    import mir.random : rne, rand, randIndex; 
}
version(stdRandom)
{
    import std.random; 
}

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

struct GillespieAlgorithm(
    Sandmann sandmann, // whether to sample another uniform value for timings
    LDM ldm, // whether to use an array for storing cumulative sum
    //size_t max=(size_t.max/2 + 1UL), 
    T
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
        The only action to perform on creation is to blow up the working cache
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

        size_t length() @nogc
        {
            return propsArr_.length; 
        }
    }

    /**
	consumes any range of propensities, where cumulativeFold is applicable at. 
	*/
	void put(P)(P props) if(isForwardRange!P)
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
        }
        else
        {
            import std.algorithm.iteration : fold;

            a0_ = resCumFold.fold!max; 
            nextReaction_ = getNextReaction(resCumFold);
        }
    }

	FloatingType tau() //@nogc
	{
        typeof(return) retVal = cast(typeof(return))1.0/cast(typeof(return))a0; 

		static if(sandmann == Sandmann.no)
		{
            version(mirRandom)
            {
                retVal = - retVal * log(rne.rand!(typeof(return)).fabs);
            }
			version(stdRandom)
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
            version(mirRandom)
            {
                const rndNum = rne.rand!T.fabs; 
            }
            version(stdRandom)
            {
                const rndNum = uniform01!T; 
            }
        }
        else
        {
            version(mirRandom)
            {
                const rndNum = rne.randIndex!size_t(a0 + 1); 
            }
            version(stdRandom)
            {
                const rndNum = (a0 + 1).iota.randomSample(1).front; 
            }
        }
		
		/**
		The usage of persistent working cache allows to accelerate the search of next reaction by binary search.
		See [3] for reference.
		*/
		static if(ldm == LDM.yes)
		{
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
            return workRange.countUntil!"a >= b"(rndNum);
		}
    }
}

import core.stdc.limits : CHAR_BIT; 
import std.algorithm.iteration : sum;
import std.algorithm.searching : any;

static foreach(e1; EnumMembers!Sandmann)
{
    static foreach(e2; EnumMembers!LDM)
    {
        static foreach(T; possibleTypes)
        {
            static foreach(l; 0 .. CHAR_BIT)
            {
                unittest
                {
                    T[] inputProps; 
                    inputProps.length = l; 

                    GillespieAlgorithm!(e1, e2, T) s; 

                    static if(l)
                    {                        
                        static if(e2 == LDM.yes)
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
                            inputProps.fill; 
                            assert(inputProps.sum > 0, ElementType!(typeof(inputProps)).stringof); 
                            put(s, inputProps); 
                            auto tau = s.tau; 
                            
                            static if(e1 == Sandmann.yes)
                            {
                                assert(approxEqual(tau, 1.0/inputProps.sum));
                            }
                            else
                            {
                                res += (tau - 1.0/inputProps.sum); 
                            }
                        }
                        assert(approxEqual(abs(res/runs), 0));
                    }
                    else
                    {
                        import std.exception : assertThrown; 
                        static if(e2 == LDM.yes)
                        {
                            assertThrown!Error(typeof(s)(inputProps.length));
                        }
                        else
                        {
                            assertThrown!Error(put(s, inputProps));
                        }
                    }
                }
            }
        }
    }
}

version(unittest)
{
    void fill(R)(R inputProps) 
    {
        alias T = ElementType!R; 
        enum boundExp = CHAR_BIT * int.sizeof; 
        static if(isFloatingPoint!T)
        {
            version(mirRandom)
            {
                foreach(ref el; inputProps) { el = rne.rand!T(boundExp).fabs; }
            }
            version(stdRandom)
            {
                foreach(ref el; inputProps) { el = uniform01 * (1UL << boundExp); }
            }
            
            assert(!inputProps.any!isNaN); 
        }
        else
        {
            version(mirRandom)
            {
                foreach(ref el; inputProps){ el = rne.randIndex!T(T.max); }
            }
            version(stdRandom)
            {
                foreach(ref el; inputProps){ el = T.max.iota.randomSample(1).front; }
            }
            
        }
    }
    
}