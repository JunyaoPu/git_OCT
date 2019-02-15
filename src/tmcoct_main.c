#define THINKCPROFILER 0

/* GNU cc does not support difftime() and CLOCKS_PER_SEC.*/
#define GNUCC 0

#if THINKCPROFILER
#include <profile.h>
#include <console.h>
#endif

#include "tmcoct.h"
#include "time.h"

RegionStruct regionspecs[MAX_REGIONS];

void InitOutputData(SimulationStruct, OutStruct *);

void FreeData(SimulationStruct, OutStruct *);

double Rspecular(double ni, double nt);

void LaunchPhoton(SimulationStruct *in_parm, double Rspecular, TetrahedronPhotonStruct *Photon_Ptr,
                  TetrahedronPhotonStruct *PhotonCont_Ptr, TetrahedronStruct *TetrahedronRoot);

void HopDropSpin(SimulationStruct *, TetrahedronPhotonStruct *, TetrahedronPhotonStruct *, OutStruct *);

Boolean HitBoundary(TetrahedronPhotonStruct *Photon_Ptr);


time_t PunchTime(char F, char *Msg) {
#if GNUCC
    return(0);
#else
    static clock_t ut0;    /* user time reference. */
    static time_t rt0;    /* real time reference. */
    double secs;
    char s[STR_LEN];

    if (F == 0) {
        ut0 = clock();
        rt0 = time(NULL);
        return (0);
    } else if (F == 1) {
        secs = (clock() - ut0) / (double) CLOCKS_PER_SEC;
        if (secs < 0) secs = 0;    /* clock() can overflow. */
        sprintf(s, "User time: %8.0lf sec = %8.2lf hr.  %s\n",
                secs, secs / 3600.0, Msg);
        puts(s);
        strcpy(Msg, s);
        return (difftime(time(NULL), rt0));
    } else if (F == 2) return (difftime(time(NULL), rt0));
    else return (0);
#endif
}

/***********************************************************
 *	Print the current time and the estimated finishing time.
 *
 *	P1 is the number of computed photon packets.
 *	Pt is the total number of photon packets.
 **********************************************************/
void PredictDoneTime(long P1, long Pt) {
    time_t now, done_time;
    struct tm *date;
    char s[80];

    now = time(NULL);
    date = localtime(&now);
    strftime(s, 80, "%H:%M %x", date);
    printf("Now %s, ", s);

    done_time = now +
                (time_t) (PunchTime(2, "") / (double) P1 * (Pt - P1));
    date = localtime(&done_time);
    strftime(s, 80, "%H:%M %x", date);
    printf("End %s\n", s);
}

/***********************************************************
 *	Execute Monte Carlo simulation for one independent run.
 **********************************************************/
void DoOneRun(short NumRuns, SimulationStruct *In_Ptr, TetrahedronStruct *RootTetrahedron) {
    register long i_photon;                                    /* index to photon. register for speed.*/
    OutStruct out_parm;                                        /* distribution of photons.*/
    TetrahedronPhotonStruct photon;                            /* The photon that is being traced */
    TetrahedronPhotonStruct photon_cont;                       /* The secondary photon due to the splitting in biased scattering which will be trace */

    long num_photons = In_Ptr->number_of_photons;

#if THINKCPROFILER
    InitProfile(200,200); cecho2file("prof.rpt",0, stdout);
#endif

    InitOutputData(*In_Ptr, &out_parm);                        /* initiating input and output structure */

    /* The specular reflection at the beginning */
    out_parm.Rsp = Rspecular(regionspecs[0].n, regionspecs[RootTetrahedron->region].n);

    i_photon = num_photons;

    PunchTime(0, "");

    photon.FstBackReflectionFlag = 0;  //No backreflection

    do {
        if (((num_photons - i_photon) % 10000) == 0 && i_photon < num_photons) {
            printf("%ld photons & %hd runs left, ", i_photon, NumRuns);
            PredictDoneTime(num_photons - i_photon, num_photons);
        }
        LaunchPhoton(In_Ptr, out_parm.Rsp, &photon, &photon_cont, RootTetrahedron);
        do {
            HopDropSpin(In_Ptr, &photon, &photon_cont, &out_parm);
        } while (!photon.dead);

        if (photon.FstBackReflectionFlag && In_Ptr->TypeBias != 3)
            i_photon++;

    } while (--i_photon);

    WriteResult(In_Ptr, &out_parm);
#if THINKCPROFILER
    exit(0);
#endif

}

int TMCOCT_run(int argc, char **argv) {

    char *filename_opt = NULL;
    char *filename_mesh = NULL;
    char *filename_bias = NULL;

    int n_simulations;
    SimulationStruct *simulations;

    int number_of_vertices = 0;
    int number_of_faces = 0;
    int number_of_tetrahedrons = 0;
    int mesh_param; // check vertices and no. of tetrahedrons
    int bias_flag;  // to check for errors in bias reading

    // Parse command-line arguments.
    if (interpret_arg(argc, argv, &filename_opt, &filename_mesh, &filename_bias)) {
        usage(argv[0]);
        return 1;
    }

    /*********************************
    ** Read the simulation inputs.
    *********************************/
    n_simulations = read_simulation_data(filename_opt, &simulations);
    if (n_simulations == 0) {
        printf("Something wrong with read_simulation_data!\n");
        return 1;
    }
    printf("Read %d simulations\n", n_simulations);

    mesh_param = read_mesh_param(filename_mesh, &number_of_vertices,
                                 &number_of_tetrahedrons);
    if (mesh_param == 0) {
        printf("Something wrong with the meshing parameters!\n");
        return 1;
    }

    printf("\n====================================\n");
    printf("Mesh Parameters:\n");
    printf("No. of nodes=%d\n", number_of_vertices);
    printf("No. of tetrahedrons=%d\n", number_of_tetrahedrons);
    number_of_faces = 4 * number_of_tetrahedrons;
    printf("No. of faces=%d\n", number_of_faces);
    printf("====================================\n");

    struct TetrahedronStruct *Graph;
    struct TetrahedronStruct *Root;

    printf("Constructing mesh graph. Please wait...\n");
    clock_t start = clock(), diff;
    Graph = MeshGraph(filename_mesh, &simulations);
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Mesh graph is created in %d seconds %d milliseconds.\n", msec / 1000, msec % 1000);

    /******************
    **  Regions
    ******************/

    for (uint i = 0; i < simulations->n_regions; i++) {
        regionspecs[i].n = simulations->regions[i].n;
        regionspecs[i].g = simulations->regions[i].g;
        regionspecs[i].muas = simulations->regions[i].muas;
        regionspecs[i].rmuas = simulations->regions[i].rmuas;
        regionspecs[i].mua_muas = simulations->regions[i].mua_muas;
    }

    for (simulations->probe.current_scanning_x = simulations->probe.start_x;
         simulations->probe.current_scanning_x < simulations->probe.end_x;
         simulations->probe.current_scanning_x += simulations->probe.distance_Ascans) {

        // probe locations
        simulations->probe_x = simulations->probe.current_scanning_x;
        simulations->probe_y = 0;
        simulations->probe_z = 0;

        printf("Finding the root of the graph...\n");

        Root = FindRootofGraph(Graph, &simulations);

        printf("Reading bias parameters.\n\n");
        bias_flag = Bias_Read(filename_bias, &simulations);
        if (bias_flag == 0) {
            printf("Something wrong with the bias reading, check file!\n");
            return 1;
        }

        clock_t start = clock(), diff;
        DoOneRun(n_simulations, simulations, Root);
        diff = clock() - start;
        float elapsedTime = diff * 1000 / CLOCKS_PER_SEC;
        printf("\n\n>>>>>>Simulation time: %f (ms)\n", elapsedTime);
    }

    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {

    ShowVersion("Version 1, 2016");

    return (TMCOCT_run(argc, argv));

}
