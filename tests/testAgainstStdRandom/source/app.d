import gillespied;
import core.stdc.limits : CHAR_BIT;  
import std.traits : EnumMembers;
import std.range; 

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