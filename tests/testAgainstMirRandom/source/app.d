import mir.random; 
import gillespied;
import core.stdc.limits : CHAR_BIT;  
import std.traits : EnumMembers;
import std.range; 

import std.meta : AliasSeq;
//alias typesToCheck = possibleTypes; 

debug = bug; // <-- uncomment to turn bug on

debug(bug)
{
    alias typesToCheck = AliasSeq!(real);
}
else
{
    alias typesToCheck = AliasSeq!(float, double, size_t);
}

void main()
{
    static foreach(e1; EnumMembers!Sandmann)
    {
        static foreach(e2; EnumMembers!LDM)
        {
            static foreach(T; typesToCheck)
            {
                static foreach(l; 0 .. CHAR_BIT)
                {
                    testTemplate!(e1, e2, T, l); 
                }
            }
        }
    }
}