#include "linker.cpp"
#include <iostream>



int main(int argc, char* argv[])
{
   KINCFactory factory;
   return Ace::run("KINC",factory,argc,argv);
}
