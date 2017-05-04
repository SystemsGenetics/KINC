#include "linker.cpp"
#include <iostream>



#ifndef UNIT_TEST



int main(int argc, char* argv[])
{
   HelpFactory::getInstance().addItem(std::unique_ptr<AbstractHelpItem>(new EMatrixHelpItem));
   HelpFactory::getInstance().addItem(std::unique_ptr<AbstractHelpItem>(new CMatrixHelpItem));
   HelpFactory::getInstance().addItem(std::unique_ptr<AbstractHelpItem>(new SpearmanHelpItem));
   HelpFactory::getInstance().addItem(std::unique_ptr<AbstractHelpItem>(new RMTHelpItem));
   HelpFactory::getInstance().addItem(std::unique_ptr<AbstractHelpItem>(new PearsonHelpItem));
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
