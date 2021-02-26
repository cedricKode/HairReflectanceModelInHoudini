#ifndef HAITCONVER_H
#define HAITCONVER_H
#include <vector>
#include "cyHairFile.h"

class HairConvert {
public:
    

    HairConvert(const char *filenameIn, const char *filenameOut);
    virtual ~HairConvert();
    
private:
    
    cyHairFile hair;
    float *dirs;

    const char *m_filenameIn;
    const char *m_filenameOut;
    
    void LoadHairModel( const char *filename, cyHairFile &hairfile, float *&dirs );

};

#endif
