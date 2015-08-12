#ifndef _INDEXER_
#define _INDEXER_

class Indexer {
  protected:
    char * indexdir;

  public:
    Indexer(char * indexdir);
    ~Indexer();

    void run();
};

#endif
