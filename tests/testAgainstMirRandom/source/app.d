import mir.random; 
import gillespied;
import core.stdc.limits : CHAR_BIT;  
import std.traits : EnumMembers;
import std.range; 

//debug = bug; // <-- uncomment to turn bug on
import std.meta : AliasSeq;
debug(bug)
{
    alias possibleTypes = AliasSeq!(real);
}
else
{
    import std.meta : AliasSeq;
    alias possibleTypes = AliasSeq!(float, double, size_t);
}

void main()
{
    static foreach(e1; EnumMembers!Sandmann)
    {
        static foreach(e2; EnumMembers!LDM)
        {
            static foreach(T; possibleTypes)
            {
                static foreach(l; 0 .. CHAR_BIT)
                {
                    testTemplate!(e1, e2, T, l); 
                }
            }
        }
    }
}