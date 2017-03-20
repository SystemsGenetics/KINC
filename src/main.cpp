#include "linker.cpp"
#include <iostream>



#ifndef UNIT_TEST



int main(int argc, char* argv[])
{
   KINCFactory factory;
   return Ace::run("KINC",factory,argc,argv);
}



#else



void unit_test(Ace::Console& console)
{
   Spearman::runUnitTests(console);
}



int main(int argc, char *argv[])
{
   KINCFactory factory;
   return Ace::run("KINC",factory,argc,argv,unit_test);
}



#endif
