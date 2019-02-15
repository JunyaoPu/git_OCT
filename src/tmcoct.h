#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <stdbool.h>

#define PI 3.14159265359
#define WEIGHT 1E-4        /* Critical weight for roulette. */
#define CHANCE 0.1        /* Chance of roulette survival. */

#define STR_LEN 256        /* String length. */

#define Boolean char

#define SIGN(x) ((x)>=0 ? 1:-1)

#define sq(x) ((x)*(x) )
#define cube(x) ( (x)*(x)*(x) )
#define NUM_SUBSTEPS_RESOLUTION 6


// The max number of regions supported (MAX_REGIONS including 1 ambient region)
#define MAX_REGIONS 100

/****************** Stuctures *****************************/

// Vertices coordinates
typedef struct {
    double x, y, z;
} Vertex;

// Faces formed with triangles
typedef struct {
    double Nx, Ny, Nz;
    double d;
    Vertex V[3];
} TriangleFaces;

typedef struct TetrahedronStruct TetrahedronStruct;

struct TetrahedronStruct {
    int region;
    TriangleFaces Faces[4];
    TetrahedronStruct *adjTetrahedrons[4];
    Vertex Vertices[4];
};

// Probe Position Data
typedef struct {
    double current_scanning_x;
    double start_x;
    double end_x;
    double distance_Ascans;
} ProbeData;

typedef struct {
    double x, y, z;                        /* Cartesian coordinates.[cm] */
    double ux, uy, uz;                    /* directional cosines of a photon. */
    double w;                                /* weight. */

    Boolean dead;                            /* 1 if photon is terminated. */
    TetrahedronStruct *tetrahedron;        /* index to layer where the photon */
    int NextTetrahedron;
    /* packet resides. */
    double s;                        /* current step size. [cm]. */
    double sleft;                    /* step size left. dimensionless [-]. */
    double MinCos;
    
    double OpticalPath;
    double MaxDepth;
    double LikelihoodRatio;
    double LikelihoodRatioAfterFstBias;
    double LikelihoodRatioIncreaseFstBias;

    short FstBackReflectionFlag;
    double LocationFstBias;
    unsigned int NumBackwardsSpecularReflections;

} TetrahedronPhotonStruct;

typedef struct RegionStruct RegionStruct;
struct RegionStruct {
    double n;                  // refractive index of a region
    double muas;               // mua + mus
    double rmuas;              // 1/(mua+mus) = mutr = 1 / mua+mus
    double mua_muas;           // mua/(mua+mus)
    double g;                  // anisotropy
};

typedef struct {

    double probe_x, probe_y, probe_z;

    char outp_filename[STR_LEN];
    char inp_filename[STR_LEN];

    // the starting and ending offset (in the input file) for this simulation
    long begin, end;
    // ASCII or binary output
    char AorB;

    uint number_of_photons;

    double Rspecular;

    // Tetrahedrons variables
    short num_tetrahedrons;
    ProbeData probe;

    // Bias variables
    short int TypeSimulation;
    short int TypeBias;
    double BackwardBiasCoefficient;
    double TargetOpticalDepth;
    double CoherenceLengthSource;

    double TargetDepthMin;
    double TargetDepthMax;

    short RecordPhotonsFlag;
    double MaxCollectingRadius;
    double MaxCollectingAngleDeg;

    double ProbabilityAdditionalBias;
    double MaxRelativeContributionToBinPerPostProcessedSample;

    FILE *output_file_ptr; // check what's used in simulation from this file
    long int NumFilteredPhotons;
    long int NumFilteredPhotonsClassI;
    long int NumFilteredPhotonsClassII;

    double MaxFilteredOpticalDepth;
    short NumOpticalDepthLengthSteps;
    double OpticalDepthShift;
    double rndStepSizeInTissue;

    uint n_regions; // Number of regions including 1 region (with index 0) for ambient medium
    RegionStruct *regions;
} SimulationStruct;

typedef struct {

    double Rsp;    /* specular reflectance. [-] */

    double *ReflectanceClassI;
    double *ReflectanceClassI_Max;
    double *ReflectanceClassI_Sum;
    double *ReflectanceClassI_SumSq;
    unsigned long int *NumClassI_PhotonsFilteredInRange;

    double *ReflectanceClassII;
    double *ReflectanceClassII_Max;
    double *ReflectanceClassII_Sum;
    double *ReflectanceClassII_SumSq;
    unsigned long int *NumClassII_PhotonsFilteredInRange;

    double *MeanReflectanceClassI_Sum;
    double *MeanReflectanceClassI_SumSq;
    double *MeanReflectanceClassII_Sum;
    double *MeanReflectanceClassII_SumSq;

} OutStruct;

double *AllocVector(short, short);

double **AllocMatrix(short, short, short, short);

double ***AllocHyperMatrix(short, short, short, short, short, short);

void FreeVector(double *, short, short);

void FreeMatrix(double **, short, short, short, short);

void FreeHyperMatrix(double ***, short, short, short, short, short, short);

void nrerror(char *);

void Roulette(TetrahedronPhotonStruct *Photon_Ptr);

void ReadBiasParm(FILE *File_Ptr,
                  SimulationStruct *In_Ptr);

void WriteBiasParm(FILE *output_bias_file_ptr,
                   SimulationStruct *In_Ptr);

void CopyPhotonStruct(TetrahedronPhotonStruct *OrigPhoton_Ptr, TetrahedronPhotonStruct *DestPhoton_Ptr);

void Spin(double g,
          TetrahedronPhotonStruct *Photon_Ptr);

// IO functions to be used with mesh and bias files
extern int readints(int n_ints, int *temp, FILE *pFile);

extern int readfloats(int n_floats, float *temp, FILE *pFile);

extern int read_simulation_data(char *filename, SimulationStruct **simulations);

// Meshing
extern int read_mesh_param(char *filename, int *number_of_vertices,
                           int *number_of_tetrahedrons);

extern int interpret_arg(int argc, char *argv[], char **fpath_o, char **fpath_m,
                         char **fpath_b);

extern void usage(const char *prog_name);

extern void ShowVersion(char const *version);

extern TetrahedronStruct *MeshGraph(char *filename, SimulationStruct **simulations);

extern TetrahedronStruct *FindRootofGraph(TetrahedronStruct *tetrahedrons, SimulationStruct **simulations);

extern RegionStruct regionspecs[];

extern void WriteReflectanceClassI_and_ClassII_Filtered(SimulationStruct *In_Ptr,
                                                        OutStruct *Out_Ptr,
                                                        FILE *OutFileClassI_Scatt,
                                                        FILE *OutFileClassII_Scatt);

extern void WriteReflectanceClassI_and_ClassII(SimulationStruct *In_Ptr,
                                               OutStruct *Out_Ptr,
                                               FILE *OutFileClassI_Scatt,
                                               FILE *OutFileClassII_Scatt);

extern int Bias_Read(char *filename, SimulationStruct **simulations);

extern void WriteResult(SimulationStruct *In_Parm, OutStruct *Out_Parm);
