#include <iostream>
#include "Convert.h"
#include "cyHairFile.h"


HairConvert::HairConvert(const char *filenameIn, const char *filenameOut)
{
    // assign file name to read
    m_filenameIn = filenameIn;
    m_filenameOut = filenameOut;
    LoadHairModel(filenameIn, hair, dirs);
}

HairConvert::~HairConvert()
{
    delete[]dirs;
}

void HairConvert::LoadHairModel( const char *filename, cyHairFile &hairfile, float *&dirs )
{

    // --------------------------------------------------------------------
    // -- The following section is based on the format by Cem Yuksel     --
    //          -- www.cemyuksel.com/research/hairmodels --
    // --------------------------------------------------------------------

    int result = hairfile.LoadFromFile( filename );

    switch( result ) {
        case CY_HAIR_FILE_ERROR_CANT_OPEN_FILE:
            printf("Error: Cannot open hair file!\n");
            return;
        case CY_HAIR_FILE_ERROR_CANT_READ_HEADER:
            printf("Error: Cannot read hair file header!\n");
            return;
        case CY_HAIR_FILE_ERROR_WRONG_SIGNATURE:
            printf("Error: File has wrong signature!\n");
            return;
        case CY_HAIR_FILE_ERROR_READING_SEGMENTS:
            printf("Error: Cannot read hair segments!\n");
            return;
        case CY_HAIR_FILE_ERROR_READING_POINTS:
            printf("Error: Cannot read hair points!\n");
            return;
        case CY_HAIR_FILE_ERROR_READING_COLORS:
            printf("Error: Cannot read hair colors!\n");
            return;
        case CY_HAIR_FILE_ERROR_READING_THICKNESS:
            printf("Error: Cannot read hair thickness!\n");
            return;
        case CY_HAIR_FILE_ERROR_READING_TRANSPARENCY:
            printf("Error: Cannot read hair transparency!\n");
            return;
        default:
            printf("Hair file \"%s\" loaded.\n", filename);
    }
    size_t hairCount = hairfile.GetHeader().hair_count;
    size_t pointCount = hairfile.GetHeader().point_count;
    printf("Number of hair strands = %d\n", int(hairCount) );
    printf("Number of hair points = %d\n", int(pointCount) );
    
    dirs = new float[3 * pointCount];

    if( hairfile.FillDirectionArray( dirs ) == 0 )
    {
        printf("Error: Cannot compute hair directions!\n");
    }
    /// end of Citation

    FILE *fp;
    fp = fopen (m_filenameOut, "w");
    const unsigned short *segments = hairfile.GetSegmentsArray();
    const float* points = hairfile.GetPointsArray();

    size_t index = 0;
    if ( segments )
    {
        // If segments array exists
        for ( size_t hairIndex=0; hairIndex <= hairCount/3 + 1; hairIndex++ )
        {
            size_t numSegments = segments[ hairIndex ];
            for (size_t j = 0; j <= 3*numSegments; j+=3)
            {
                for (size_t k = 0; k < 3; k++)
                {
                    std::string s = std::to_string(points[index + j + k]);
                    char const *pchar = s.c_str();
                    fputs(pchar, fp);
                    fputs(" ", fp);
                }
            }
            fputs("\n", fp);
            index += 3*numSegments + 3;
        }
    }
    else
    {
        // If segments array does not exist, use default segment count
        size_t dsegs = hairfile.GetHeader().d_segments;
        for ( size_t hairIndex=0; hairIndex <= hairCount/3 + 1; hairIndex++ )
        {
            size_t numSegments = dsegs;
            for (size_t j = 0; j <= 3*numSegments; j+=3)
            {
                for (size_t k = 0; k < 3; k++)
                {
                    std::string s = std::to_string(points[index + j + k]);
                    char const *pchar = s.c_str();
                    fputs(pchar, fp);
                    fputs(" ", fp);
                }
            }
            fputs("\n", fp);
            index += 3*numSegments + 3;
        }
    }

    fclose (fp);
}
