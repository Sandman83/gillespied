import gillespied;
import std.traits : EnumMembers;
import std.range; 
import mir.random; 

void main()
{
    static foreach(e1; EnumMembers!Sandmann)
    {
        static foreach(e2; EnumMembers!LDM)
        {
            static foreach(T; possibleTypes)
            {
                static foreach(l; minTestedSizes .. maxTestedSizes) 
                {
                    testTemplate!(e1, e2, T, l); 
                }
            }
        }
    }
}