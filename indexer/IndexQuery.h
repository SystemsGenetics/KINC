#ifndef _INDEXQUERY_
#define _INDEXQUERY_

class IndexQuery {
  protected:
    char * indexdir;

  public:
    IndexQuery(char * indexdir);
    ~IndexQuery();

    void run();
};
#endif
