#include "tmcoct.h"
#include <gsl/gsl_rng.h>

const double G_COS_0_D = 1.0 - 1.0E-14;
const double G_COS_90_D = 1.0E-7;


#define STANDARDTEST 0
/* testing program using fixed rnd seed. */

#define PARTIALREFLECTION 1
/* 1=split photon, 0=statistical reflection. */

#define COSZERO (1.0-1.0E-14)
/* cosine of about 1e-6 rad. */

#define COS90D  1.0E-6
/* cosine of about 1.57 - 1e-6 rad. */

const gsl_rng_type *T_RNG_GSL;
gsl_rng *r_RNG_GSL;

/***********************************************************
 *	A random number generator from Numerical Recipes in C.
 ****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


/***********************************************************
 *	Generate a random number between 0 and 1.  Take a 
 *	number as seed the first time entering the function.  
 *	The seed is limited to 1<<15.  
 *	We found that when idum is too large, ran3 may return 
 *	numbers beyond 0 and 1.
 ****/
double RandomNum(void) {
    static Boolean first_time = 1;
    static unsigned long int idum;    /* seed for ran3. */

    if (first_time) {
        T_RNG_GSL = gsl_rng_default;
        r_RNG_GSL = gsl_rng_alloc(T_RNG_GSL);

#if STANDARDTEST /* Use fixed seed to test the program. */
        printf("*** Using seed to test. ***\n");
    printf("*** Remove the fixed seed after the test. ***\n");
    idum =  7923;

#else
        idum = (unsigned long int) time(NULL) % (1 << 15);
        /* use 16-bit integer as the seed. */
#endif
        printf("Random seed = %lu \n", idum);

        //ran3(&idum);
        gsl_rng_set(r_RNG_GSL, idum);

        first_time = 0;
        //idum = 1;
    }
    // TODO FIXME remove 1.0 -
    return (1.0 - gsl_rng_uniform(r_RNG_GSL));
}

/***********************************************************
 *	Compute the specular reflection. 
 *
 *	If the first layer is a turbid medium, use the Fresnel
 *	reflection from the boundary of the first layer as the 
 *	specular reflectance.
 *
 *	If the first layer is glass, multiple reflections in
 *	the first layer is considered to get the specular
 *	reflectance.
 *
 *	The subroutine assumes the Layerspecs array is correctly 
 *	initialized.
 ****/
double Rspecular(double ni, double nt) {
    double r1;
    /* direct reflections from the 1st and 2nd layers. */
    double temp;
    temp = (ni - nt) / (ni + nt);
    r1 = temp * temp;
    return (r1);
}

/***********************************************************
 *	Initialize a photon packet.
 ****/
void LaunchPhoton(SimulationStruct *in_parm, double Rspecular, TetrahedronPhotonStruct *Photon_Ptr,
                  TetrahedronPhotonStruct *PhotonCont_Ptr, TetrahedronStruct *TetrahedronRoot) {
    if (Photon_Ptr->FstBackReflectionFlag) {
        //Most often case for biased simulations
        //The photon will be re-initialized from where the bias was applied
        long double LikelihoodRatioTmp = Photon_Ptr->LikelihoodRatioAfterFstBias;
        CopyPhotonStruct(PhotonCont_Ptr, Photon_Ptr);

        //Apply one unbiased random scatter to the split photon
        Spin(regionspecs[Photon_Ptr->tetrahedron->region].g, Photon_Ptr);

        if (LikelihoodRatioTmp < 1)
            // Adjust the likelihoodRation to make the simulation unbiased
            Photon_Ptr->LikelihoodRatio = 1 - LikelihoodRatioTmp;
        else
            Photon_Ptr->LikelihoodRatio = 1;

    } else {

        Photon_Ptr->w = 1.0 - Rspecular;

        Photon_Ptr->x = in_parm->probe_x;
        Photon_Ptr->y = in_parm->probe_y;
        Photon_Ptr->z = in_parm->probe_z;

        Photon_Ptr->s = 0.0;
        Photon_Ptr->sleft = 0.0;

        Photon_Ptr->ux = Photon_Ptr->uy = 0.0;
        Photon_Ptr->uz = 1.0;
        Photon_Ptr->dead = 0.0;

        Photon_Ptr->OpticalPath = 0.0;
        Photon_Ptr->MaxDepth = 0.0;
        Photon_Ptr->LikelihoodRatio = 1.0;
        Photon_Ptr->LikelihoodRatioAfterFstBias = 1.0;

        Photon_Ptr->FstBackReflectionFlag = 0;
        Photon_Ptr->LocationFstBias = -1.0;
        Photon_Ptr->NumBackwardsSpecularReflections = 0;

        Photon_Ptr->LikelihoodRatioIncreaseFstBias = 0;

        Photon_Ptr->LikelihoodRatioIncreaseFstBias = 0;

        Photon_Ptr->NextTetrahedron = -1;

        Photon_Ptr->tetrahedron = TetrahedronRoot;
    }
}

/***********************************************************
 *	Choose (sample) a new theta angle for photon propagation
 *	according to the anisotropy.
 *
 *	If anisotropy g is 0, then
 *		cos(theta) = 2*rand-1.
 *	otherwise
 *		sample according to the Henyey-Greenstein function.
 *
 *	Returns the cosine of the polar deflection angle theta.
 ****/
double SpinTheta(double g) {
    double cost;

    if (g == 0.0)
        cost = 2.0 * RandomNum() - 1.0;
    else {
        double temp = (1.0 - g * g) / (1.0 - g + 2.0 * g * RandomNum());
        cost = (1.0 + g * g - temp * temp) / (2.0 * g);
        if (cost < -1.0) cost = -1.0;
        else if (cost > 1.0) cost = 1.0;
    }
    return (cost);
}

double SpinThetaForwardFstBias(double g) {
    double cost;

    if (g == 0.0)
        cost = RandomNum();
    else {
        double RandomNumTmp = RandomNum();
        double temp = RandomNumTmp / (1 - g) + (1 - RandomNumTmp) / sqrt(g * g + 1);
        cost = (g * g + 1 - 1. / (temp * temp)) / (2 * g);
        if (cost < -1) cost = -1;
        else if (cost > 1) cost = 1;
    }
    return (cost);
}

/***********************************************************
 *	Choose a new direction for photon propagation by 
 *	sampling the polar deflection angle theta and the 
 *	azimuthal angle psi.
 *
 *	Note:
 *  	theta: 0 - pi so sin(theta) is always positive 
 *  	feel free to use sqrt() for cos(theta).
 * 
 *  	psi:   0 - 2pi 
 *  	for 0-pi  sin(psi) is + 
 *  	for pi-2pi sin(psi) is - 
 ****/

void Spin(double g,
          TetrahedronPhotonStruct *Photon_Ptr) {
    double cost, sint;    /* cosine and sine of the */
    /* polar deflection angle theta. */
    double cosp, sinp;    /* cosine and sine of the */
    /* azimuthal angle psi. */
    double ux = Photon_Ptr->ux;
    double uy = Photon_Ptr->uy;
    double uz = Photon_Ptr->uz;
    double psi;

    cost = SpinTheta(g);
    sint = sqrt(1.0 - cost * cost);
    /* sqrt() is faster than sin(). */

    psi = 2.0 * PI * RandomNum(); /* spin psi 0-2pi. */
    cosp = cos(psi);
    if (psi < PI)
        sinp = sqrt(1.0 - cosp * cosp);
        /* sqrt() is faster than sin(). */
    else
        sinp = -sqrt(1.0 - cosp * cosp);

    if (fabs(uz) > COSZERO) {    /* normal incident. */
        Photon_Ptr->ux = sint * cosp;
        Photon_Ptr->uy = sint * sinp;
        Photon_Ptr->uz = cost * SIGN(uz);
        /* SIGN() is faster than division. */
    } else {
        /* regular incident. */
        double temp = sqrt(1.0 - uz * uz);
        Photon_Ptr->ux = sint * (ux * uz * cosp - uy * sinp)
                         / temp + ux * cost;
        Photon_Ptr->uy = sint * (uy * uz * cosp + ux * sinp)
                         / temp + uy * cost;
        Photon_Ptr->uz = -sint * cosp * temp + uz * cost;
    }

}

void SpinBias37(double g,
                TetrahedronPhotonStruct *Photon_Ptr,
                TetrahedronPhotonStruct *PhotonCont_Ptr,
                SimulationStruct *In_Ptr)
{
    double g_squared = sq(g);

    double cost, sint;    /* cosine and sine of the */
    /* polar deflection angle theta. */
    double cosp, sinp;    /* cosine and sine of the */
    /* azimuthal angle psi. */
    double ux = Photon_Ptr->ux;
    double uy = Photon_Ptr->uy;
    double uz = Photon_Ptr->uz;
    double ux_Orig = ux;
    double uy_Orig = uy;
    double uz_Orig = uz;
    double psi;


    double BackwardBiasCoefficient = In_Ptr->BackwardBiasCoefficient;
    double costg;
    double costg1, costg2;

    double BiasCoefficientTmp = 0;

    short int ReachedTargetOpticalDepthFlag = 0;
    short int ThisIsFirstBackwardBiasFlag = 0;

    if (Photon_Ptr->z > In_Ptr->TargetDepthMin
        && Photon_Ptr->z < In_Ptr->TargetDepthMax
        && Photon_Ptr->uz > 0
        && !Photon_Ptr->FstBackReflectionFlag) {
        // The bias backwards will be applied only if the photon is going forward
        // The status of the photon prior to the bias will be saved

        CopyPhotonStruct(Photon_Ptr, PhotonCont_Ptr);
        ReachedTargetOpticalDepthFlag = 1;
        Photon_Ptr->FstBackReflectionFlag = 1; // Bias backwards only once
        Photon_Ptr->LocationFstBias = Photon_Ptr->z;
        BiasCoefficientTmp = BackwardBiasCoefficient;
        ThisIsFirstBackwardBiasFlag = 1;
    }

    /**********************************
	 ** Biased Direction towards probe
	 **********************************/
    double vx = In_Ptr->probe_x-Photon_Ptr->x;
    double vy = In_Ptr->probe_y-Photon_Ptr->y;
    double vz = In_Ptr->probe_z-Photon_Ptr->z;
    double LengthVector = sqrt(sq(vx) + sq(vy) + sq(vz));
    vx /= LengthVector;
    vy /= LengthVector;
    vz /= LengthVector;
    /*********************************/

    if ((Photon_Ptr->FstBackReflectionFlag
         || Photon_Ptr->NumBackwardsSpecularReflections > 0)
        && !ReachedTargetOpticalDepthFlag)
    {
        // It was biased at least once before, and is moving backwards
        ReachedTargetOpticalDepthFlag = 2;

        double NextStepSize = -log(In_Ptr->rndStepSizeInTissue) / regionspecs[Photon_Ptr->tetrahedron->region].muas;
        double CurrentDistanceToOrigin = sqrt(sq(Photon_Ptr->x-In_Ptr->probe_x) + sq(Photon_Ptr->y-In_Ptr->probe_y) + sq(Photon_Ptr->z-In_Ptr->probe_z));

        if (NextStepSize >= CurrentDistanceToOrigin
            && acos(-vz) <= In_Ptr->MaxCollectingAngleDeg * PI / 180) {
            ReachedTargetOpticalDepthFlag = 1;
        }

        BiasCoefficientTmp = BackwardBiasCoefficient;

    }

    int BiasFunctionRandomlySelected = 0;
    if (ReachedTargetOpticalDepthFlag) {
        // Photon reached target optical layer or may undergo an additional biased
        // scattering or unbiased scattering

        // BiasFunctionRandomlySelected=1 means use biased scattering and 2 means unbiased scattering
        if (RandomNum() <= In_Ptr->ProbabilityAdditionalBias)
            BiasFunctionRandomlySelected = 1;
        else
            BiasFunctionRandomlySelected = 2;

        if (ReachedTargetOpticalDepthFlag == 1
            || BiasFunctionRandomlySelected == 1) {
            /*************************************************************************
			** The photon is within the target depth and going forward
			** The additional biased scattering is randomly chosen
			** So the scattering is biased Henyey-Greenstein scattering
			*************************************************************************/
            cost = SpinThetaForwardFstBias(BiasCoefficientTmp);
            ux = vx;
            uy = vy;
            uz = vz;
        } else {
            /**************************************************************************************
			** The photon is within the target depth but the scattering is randomly selected is
			** unbiased scattering
			** or the photon is already going backward or it is out of target depth
			**************************************************************************************/
            cost = SpinTheta(g);
        }
    } else {
        /**************************************************************************
		**  The photon is not within the target depth or it is not going forward
		**  so do unbiased scattering
		**************************************************************************/
        cost = SpinTheta(g);
    }

    if (cost < -1)
        cost = -1;
    else if (cost > 1)
        cost = 1;

    sint = sqrt(1.0 - cost * cost);
    /* sqrt() is faster than sin(). */

    psi = 2.0 * PI * RandomNum(); /* spin psi 0-2pi. */
    cosp = cos(psi);
    if (psi < PI)
        sinp = sqrt(1.0 - cosp * cosp);
        /* sqrt() is faster than sin(). */
    else
        sinp = -sqrt(1.0 - cosp * cosp);

    if (fabs(uz) > COSZERO) {    /* normal incident. */
        Photon_Ptr->ux = sint * cosp;
        Photon_Ptr->uy = sint * sinp;
        Photon_Ptr->uz = cost * SIGN(uz);
        /* SIGN() is faster than division. */
    } else {        /* regular incident. */
        double temp = sqrt(1.0 - uz * uz);
        Photon_Ptr->ux = sint * (ux * uz * cosp - uy * sinp)
                         / temp + ux * cost;
        Photon_Ptr->uy = sint * (uy * uz * cosp + ux * sinp)
                         / temp + uy * cost;
        Photon_Ptr->uz = -sint * cosp * temp + uz * cost;
    }


    costg = ux_Orig * Photon_Ptr->ux
            + uy_Orig * Photon_Ptr->uy
            + uz_Orig * Photon_Ptr->uz;
    costg2 = costg;
    costg1 = vx * Photon_Ptr->ux
             + vy * Photon_Ptr->uy
             + vz * Photon_Ptr->uz;

    if (BiasCoefficientTmp) {
        double one_plus_a_squared = 1 + sq(BiasCoefficientTmp);
        double sqrt_one_plus_a_squared = sqrt(one_plus_a_squared);
        double LikelihoodRatioIncreaseFactor;
        if (ReachedTargetOpticalDepthFlag == 1)
            /****************************************************************************************
			 ** Likelihood for the first Bias scattering. Equation (8) of the paper:
			 ** Malektaji, Siavash, Ivan T. Lima, and Sherif S. Sherif. "Monte Carlo simulation of
			 ** optical coherence tomography for turbid media with arbitrary spatial distributions."
			 ** Journal of biomedical optics 19.4 (2014): 046001-046001.
			 ****************************************************************************************/
            LikelihoodRatioIncreaseFactor = ((1.0 - g_squared) * (sqrt_one_plus_a_squared - 1.0 + BiasCoefficientTmp) *
                                             sqrt(cube(one_plus_a_squared - 2.0 * BiasCoefficientTmp * cost))) /
                                            (2.0 * BiasCoefficientTmp * (1.0 - BiasCoefficientTmp) *
                                             sqrt_one_plus_a_squared * sqrt(cube(1.0 + g_squared - 2.0 * g * costg)));
        else {
            double cost1, cost2;
            if (BiasFunctionRandomlySelected == 1) {
                cost1 = cost;
                cost2 = costg2;
            } else {
                cost1 = costg1;
                cost2 = cost;
            }
            /*******************************************************************************************************
			 **  The likelihood ratio of additional biased scatterings, whether the biased or the unbiased
			 **  probability density function is randomly selected, is calculated according to the equation (9)
			 **  of the paper:
			 **  Malektaji, Siavash, Ivan T. Lima, and Sherif S. Sherif. "Monte Carlo simulation of
			 **  optical coherence tomography for turbid media with arbitrary spatial distributions."
			 **  Journal of biomedical optics 19.4 (2014): 046001-046001.
			 *******************************************************************************************************/
            double pdf1 = (sqrt_one_plus_a_squared * BiasCoefficientTmp * (1.0 - BiasCoefficientTmp)) /
                          ((sqrt_one_plus_a_squared - 1.0 + BiasCoefficientTmp) *
                           sqrt(cube(one_plus_a_squared - 2.0 * BiasCoefficientTmp * cost1)));

            double pdf2 = (1.0 - g_squared) / (2.0 * cube(sqrt(1.0 + g_squared - 2.0 * g * cost2)));

            LikelihoodRatioIncreaseFactor = pdf2 /
                                            (In_Ptr->ProbabilityAdditionalBias * pdf1 +
                                             (1.0 - In_Ptr->ProbabilityAdditionalBias) * pdf2);

        }

        Photon_Ptr->LikelihoodRatio *= LikelihoodRatioIncreaseFactor;
        if (ThisIsFirstBackwardBiasFlag == 1) {
            // In case there was a sure backward bias and that was the very first one
            Photon_Ptr->LikelihoodRatioAfterFstBias = Photon_Ptr->LikelihoodRatio;
            Photon_Ptr->LikelihoodRatioIncreaseFstBias = LikelihoodRatioIncreaseFactor;
        }
    }
}

/***********************************************************
 *	Move the photon s away in the current layer of medium.  
 ****/
void Hop(TetrahedronPhotonStruct *Photon_Ptr) {
    double s = Photon_Ptr->s;
    Photon_Ptr->x += s * Photon_Ptr->ux;
    Photon_Ptr->y += s * Photon_Ptr->uy;
    Photon_Ptr->z += s * Photon_Ptr->uz;

    Photon_Ptr->OpticalPath += s;
    if (Photon_Ptr->MaxDepth < Photon_Ptr->z)
        Photon_Ptr->MaxDepth = Photon_Ptr->z;
}
/***********************************************************
 *	If uz != 0, return the photon step size in glass, 
 *	Otherwise, return 0.
 *
 *	The step size is the distance between the current 
 *	position and the boundary in the photon direction.
 *
 *	Make sure uz !=0 before calling this function.
 ****/

/***********************************************************
 *	Pick a step size for a photon packet when it is in 
 *	tissue.
 *	If the member sleft is zero, make a new step size 
 *	with: -log(rnd)/(mua+mus).
 *	Otherwise, pick up the leftover in sleft.
 *
 *	Layer is the index to layer.
 *	In_Ptr is the input parameters.
 ****/
void StepSizeInTissue(TetrahedronPhotonStruct *Photon_Ptr) {

    if (Photon_Ptr->sleft == 0.0) {  /* make a new step. */
        double rnd;

        do rnd = RandomNum();
        while (rnd <= 0.0);    /* avoid zero. */

        Photon_Ptr->s = -log(rnd) * regionspecs[Photon_Ptr->tetrahedron->region].rmuas;
    } else {    /* take the leftover. */
        Photon_Ptr->s = Photon_Ptr->sleft * regionspecs[Photon_Ptr->tetrahedron->region].rmuas;
        Photon_Ptr->sleft = 0.0;
    }
}

/***********************************************************
 *	Check if the step will hit the boundary.
 *	Return 1 if hit boundary.
 *	Return 0 otherwise.
 *
 * 	If the projected step hits the boundary, the members
 *	s and sleft of Photon_Ptr are updated.
 ****/
Boolean HitBoundary(TetrahedronPhotonStruct *Photon_Ptr) {
    TetrahedronStruct *tetrahedron;
    tetrahedron = Photon_Ptr->tetrahedron;

    double min_distance = 1e10;
    int index_of_tetrahedron_with_min_distance = -1;

    double cos_normals_and_photon_direction[4];
    double distance_from_face_in_photon_direction;
    double perpendicular_distance;

    for (int i = 0; i < 4; i++) {
        cos_normals_and_photon_direction[i] = (Photon_Ptr->tetrahedron->Faces[i].Nx) * Photon_Ptr->ux +
                      (Photon_Ptr->tetrahedron->Faces[i].Ny) * Photon_Ptr->uy +
                      (Photon_Ptr->tetrahedron->Faces[i].Nz) * Photon_Ptr->uz;
    }

    for (int i = 0; i < 4; i++) {
        if (cos_normals_and_photon_direction[i] < 0) {
            perpendicular_distance = ((Photon_Ptr->tetrahedron->Faces[i].Nx) * Photon_Ptr->x +
                      (Photon_Ptr->tetrahedron->Faces[i].Ny) * Photon_Ptr->y +
                      (Photon_Ptr->tetrahedron->Faces[i].Nz) * Photon_Ptr->z + Photon_Ptr->tetrahedron->Faces[i].d);
            distance_from_face_in_photon_direction = -perpendicular_distance / cos_normals_and_photon_direction[i];

            if (distance_from_face_in_photon_direction < min_distance) {
                min_distance = distance_from_face_in_photon_direction;
                index_of_tetrahedron_with_min_distance = i;
                Photon_Ptr->MinCos = cos_normals_and_photon_direction[i];
            }
        }
    }

    Boolean hit_boundary;
    if (Photon_Ptr->s > min_distance) {
        /* not horizontal & crossing. */

        Photon_Ptr->sleft = (Photon_Ptr->s - min_distance) * regionspecs[tetrahedron->region].muas;
        Photon_Ptr->s = min_distance;
        hit_boundary = 1;
        Photon_Ptr->NextTetrahedron = index_of_tetrahedron_with_min_distance;

    } else {

        hit_boundary = 0;
        Photon_Ptr->NextTetrahedron = -1;

    }
    return hit_boundary;
}

/***********************************************************
 *	Drop photon weight inside the tissue (not glass).
 *
 *  The photon is assumed not dead. 
 *
 *	The weight drop is dw = w*mua/(mua+mus).
 *
 *	The dropped weight is assigned to the absorption array 
 *	elements.
 ****/
void Drop(TetrahedronPhotonStruct *Photon_Ptr) {
    double dwa;    /* absorbed weight.*/
    TetrahedronStruct *tetrahedron = Photon_Ptr->tetrahedron;

    dwa = Photon_Ptr->w * regionspecs[tetrahedron->region].mua_muas;
    Photon_Ptr->w -= dwa;
}

/***********************************************************
 *	The photon weight is small, and the photon packet tries 
 *	to survive a roulette.
 ****/
void Roulette(TetrahedronPhotonStruct *Photon_Ptr) {
    if (Photon_Ptr->w == 0.0)
        Photon_Ptr->dead = 1;
    else if (RandomNum() < CHANCE) /* survived the roulette.*/
        Photon_Ptr->w /= CHANCE;
    else
        Photon_Ptr->dead = 1;
}

/***********************************************************
 *	Compute the Fresnel reflectance.
 *
 *	Make sure that the cosine of the incident angle a1
 *	is positive, and the case when the angle is greater 
 *	than the critical angle is ruled out.
 *
 * 	Avoid trigonometric function operations as much as
 *	possible, because they are computation-intensive.
 ****/
double
RFresnel(SimulationStruct *In_Ptr, TetrahedronPhotonStruct *Photon_Ptr, double *ux_reflected, double *uy_reflected,
         double *uz_reflected, double *ux_refracted, double *uy_refracted, double *uz_refracted) {

    double r;

    double cosa, sina, sinb, cosb;

    cosa = -(Photon_Ptr->MinCos);

    double ni = regionspecs[Photon_Ptr->tetrahedron->region].n;
    double nt;
    if (Photon_Ptr->tetrahedron->adjTetrahedrons[Photon_Ptr->NextTetrahedron] == NULL)
        nt = regionspecs[0].n;
    else {

        int next_region = Photon_Ptr->tetrahedron->adjTetrahedrons[Photon_Ptr->NextTetrahedron]->region;
        nt = regionspecs[next_region].n;
    }

    sina = sqrt(1.0 - cosa * cosa);
    sinb = ni * sina / nt;
    cosb = sqrt(1.0 - sinb * sinb);

    if (cosa > G_COS_0_D) {
        r = (ni - nt)
            / (ni + nt);
        r *= r;
        cosb = 1;
    } else if (cosa < G_COS_90_D) {
        r = 1.0;
    } else {
        sina = sqrt(1 - cosa * cosa);
        sinb = ni * sina / nt;
        if (sinb >= 1) {
            r = 1.0;
        } else {
            cosb = sqrt(1.0 - sinb * sinb);
            double cosacosb = cosa * cosb;
            double sinasinb = sina * sinb;
            double sinacosb = sina * cosb;
            double cosasinb = cosa * sinb;
            double cap = cosacosb - sinasinb;
            double cam = cosacosb + sinasinb;
            double sap = sinacosb + cosasinb;
            double sam = sinacosb - cosasinb;
            r = 0.5 * sam * sam * (cam * cam + cap * cap) / (sap * sap * cam * cam);
        }
    }

    double sbdivsa = ni / nt;
    *ux_refracted = -cosb * Photon_Ptr->tetrahedron->Faces[Photon_Ptr->NextTetrahedron].Nx +
                    sbdivsa * (cosa * Photon_Ptr->tetrahedron->Faces[Photon_Ptr->NextTetrahedron].Nx + Photon_Ptr->ux);
    *uy_refracted = -cosb * Photon_Ptr->tetrahedron->Faces[Photon_Ptr->NextTetrahedron].Ny +
                    sbdivsa * (cosa * Photon_Ptr->tetrahedron->Faces[Photon_Ptr->NextTetrahedron].Ny + Photon_Ptr->uy);
    *uz_refracted = -cosb * Photon_Ptr->tetrahedron->Faces[Photon_Ptr->NextTetrahedron].Nz +
                    sbdivsa * (cosa * Photon_Ptr->tetrahedron->Faces[Photon_Ptr->NextTetrahedron].Nz + Photon_Ptr->uz);

    *ux_reflected = 2 * cosa * Photon_Ptr->tetrahedron->Faces[Photon_Ptr->NextTetrahedron].Nx + Photon_Ptr->ux;
    *uy_reflected = 2 * cosa * Photon_Ptr->tetrahedron->Faces[Photon_Ptr->NextTetrahedron].Ny + Photon_Ptr->uy;
    *uz_reflected = 2 * cosa * Photon_Ptr->tetrahedron->Faces[Photon_Ptr->NextTetrahedron].Nz + Photon_Ptr->uz;

    return r;

}

/***********************************************************
 *	Record the photon weight exiting the first layer(uz<0), 
 *	no matter whether the layer is glass or not, to the 
 *	reflection array.
 *
 *	Update the photon weight as well.
 ***********************************************************/
void RecordValidationR(double Refl, SimulationStruct *In_Ptr, TetrahedronPhotonStruct *Photon_Ptr, OutStruct *Out_Ptr) {
    double x = Photon_Ptr->x;
    double y = Photon_Ptr->y;

    short SpatialFilterFlag = 0;
    if (sqrt(sq(x-In_Ptr->probe_x) + sq(y-In_Ptr->probe_y)) <
        In_Ptr->MaxCollectingRadius        /* The position of photon should be within the radius of the probe */
        && acos(-Photon_Ptr->uz) < In_Ptr->MaxCollectingAngleDeg * PI / 180) {
        if ((In_Ptr->TypeBias == 0 || In_Ptr->TypeBias == 3))
            SpatialFilterFlag = 1;

        if ((Photon_Ptr->FstBackReflectionFlag && In_Ptr->TypeBias == 37 ||
             Photon_Ptr->LocationFstBias == Photon_Ptr->MaxDepth
             && In_Ptr->TypeBias != 37 ||
             Photon_Ptr->NumBackwardsSpecularReflections > 0))
            SpatialFilterFlag = 1;
    }

    if (SpatialFilterFlag) {

        In_Ptr->NumFilteredPhotons++;
        double FilterOpticalPath = Photon_Ptr->OpticalPath;
        double FilterOpticalDepth = FilterOpticalPath / 2.0;

        unsigned int VectorPosition = (unsigned int) (
                (FilterOpticalPath / 2.0 - In_Ptr->OpticalDepthShift)
                / (In_Ptr->CoherenceLengthSource / NUM_SUBSTEPS_RESOLUTION));

        int ii;
        for (ii = VectorPosition - NUM_SUBSTEPS_RESOLUTION / 2;
             ii < VectorPosition + NUM_SUBSTEPS_RESOLUTION / 2;
             ii++) {
            if (ii >= 0 && ii < In_Ptr->NumOpticalDepthLengthSteps) {
                short FilterClassI_Flag =
                        Photon_Ptr->MaxDepth > (FilterOpticalDepth
                                                - In_Ptr->CoherenceLengthSource / 2.0);

                double tmpPhotonContribution = Photon_Ptr->w * (1.0 - Refl) * Photon_Ptr->LikelihoodRatio;
                if (FilterClassI_Flag) {
                    if (ii == VectorPosition)
                        In_Ptr->NumFilteredPhotonsClassI++;
                    Out_Ptr->ReflectanceClassI_Sum[ii] += tmpPhotonContribution;
                    if (Out_Ptr->ReflectanceClassI_Max[ii] < tmpPhotonContribution)
                        Out_Ptr->ReflectanceClassI_Max[ii] = tmpPhotonContribution;
                    Out_Ptr->ReflectanceClassI_SumSq[ii] += sq(tmpPhotonContribution);
                    Out_Ptr->NumClassI_PhotonsFilteredInRange[ii]++;

                } else {
                    if (ii == VectorPosition)
                        In_Ptr->NumFilteredPhotonsClassII++;
                    if (Out_Ptr->ReflectanceClassII_Max[ii] < tmpPhotonContribution)
                        Out_Ptr->ReflectanceClassII_Max[ii] = tmpPhotonContribution;
                    Out_Ptr->ReflectanceClassII_Sum[ii] += tmpPhotonContribution;
                    Out_Ptr->ReflectanceClassII_SumSq[ii] += sq(tmpPhotonContribution);
                    Out_Ptr->NumClassII_PhotonsFilteredInRange[ii]++;

                }
            }
        }
    }

    Photon_Ptr->w *= Refl;

}

/***********************************************************
 *	Decide whether the photon will be transmitted  or be 
 *	reflected on the bottom boundary (uz>0) of the current 
 *	layer.
 *
 *	If the photon is transmitted, move the photon to 
 *	"layer+1". If "layer" is the last layer, record the 
 *	transmitted weight as transmittance. See comments for 
 *	CrossUpOrNot.
 *
 *	Update the photon parmameters.
 ****/
void CrossOrNot(SimulationStruct *In_Ptr, TetrahedronPhotonStruct *Photon_Ptr, OutStruct *Out_Ptr) {

    double uz_reflected = 0, ux_reflected = 0, uy_reflected = 0;    /* cosines of transmission alpha. */
    double uz_refracted = 0, ux_refracted = 0, uy_refracted = 0;

    double r = 0.0;                        /* reflectance */

    r = RFresnel(In_Ptr, Photon_Ptr, &ux_reflected, &uy_reflected, &uz_reflected, &ux_refracted, &uy_refracted,
                 &uz_refracted);

#if PARTIALREFLECTION
    if (Photon_Ptr->tetrahedron->adjTetrahedrons[Photon_Ptr->NextTetrahedron] == NULL && r < 1.0) { /* We use actual photon splitting if the photon is exiting the media
																							 ** We store the part the is exiting and the rest will be traced */
        Photon_Ptr->ux = ux_refracted;
        Photon_Ptr->uy = uy_refracted;
        Photon_Ptr->uz = uz_refracted;

        if (Photon_Ptr->uz < 0)
            RecordValidationR(r, In_Ptr, Photon_Ptr, Out_Ptr);
        else {
            Photon_Ptr->NumBackwardsSpecularReflections++;
            Photon_Ptr->w *= r;
        }
        Photon_Ptr->ux = ux_reflected;
        Photon_Ptr->uy = uy_reflected;
        Photon_Ptr->uz = uz_reflected;

    }
        /* If photon is not exiting the medium we use the statistical splitting */
    else if (RandomNum() > r) {/* transmitted */

        Photon_Ptr->tetrahedron = Photon_Ptr->tetrahedron->adjTetrahedrons[Photon_Ptr->NextTetrahedron];

        Photon_Ptr->ux = ux_refracted;
        Photon_Ptr->uy = uy_refracted;
        Photon_Ptr->uz = uz_refracted;

        Photon_Ptr->NextTetrahedron = -1;

    } else                        /* reflected. */
    {
        Photon_Ptr->NextTetrahedron = -1;

        if (Photon_Ptr->uz > 0)
            Photon_Ptr->NumBackwardsSpecularReflections++;

        Photon_Ptr->ux = ux_reflected;
        Photon_Ptr->uy = uy_reflected;
        Photon_Ptr->uz = uz_reflected;

    }
#else	/* use just statistical splitting */
    if(RandomNum() > r) {		/* transmitted. */
    
    if(Photon_Ptr->tetrahedron->adjTetrahedrons[Photon_Ptr->NextTetrahedron]==NULL) {
      
      

      Photon_Ptr->ux = ux_refracted;
      Photon_Ptr->uy = uy_refracted;
      Photon_Ptr->uz = uz_refracted;

      if (Photon_Ptr->uz<0)
            RecordValidationR(0.0, In_Ptr, Photon_Ptr, Out_Ptr);

      Photon_Ptr->dead = 1;
    }
    else {

      /* Siavash Refraction Direction Cosines */
      Photon_Ptr->ux = ux_refracted;
      Photon_Ptr->uy = uy_refracted;
      Photon_Ptr->uz = uz_refracted;
      Photon_Ptr->tetrahedron=Photon_Ptr->tetrahedron->adjTetrahedrons[Photon_Ptr->NextTetrahedron];
      Photon_Ptr->NextTetrahedron=-1;
    }
  }
  else {
    /* reflected. */

    Photon_Ptr->ux = ux_reflected;
    Photon_Ptr->uy = uy_reflected;
    Photon_Ptr->uz = uz_reflected;
    Photon_Ptr->NextTetrahedron=-1;
  }
#endif

}

/***********************************************************
 ****/
/***********************************************************
 *	Move the photon packet in glass layer.
 *	Horizontal photons are killed because they will
 *	never interact with tissue again.
 ****/
/***********************************************************
 *	Set a step size, move the photon, drop some weight, 
 *	choose a new photon direction for propagation.  
 *
 *	When a step size is long enough for the photon to 
 *	hit an interface, this step is divided into two steps. 
 *	First, move the photon to the boundary free of 
 *	absorption or scattering, then decide whether the 
 *	photon is reflected or transmitted.
 *	Then move the photon in the current or transmission 
 *	medium with the unfinished stepsize to interaction 
 *	site.  If the unfinished stepsize is still too long, 
 *	repeat the above process.  
 ****/
void HopDropSpinInTissue(SimulationStruct *In_Ptr,
                         TetrahedronPhotonStruct *Photon_Ptr,
                         TetrahedronPhotonStruct *PhotonCont_Ptr,
                         OutStruct *Out_Ptr) {
    StepSizeInTissue(Photon_Ptr);
    if (HitBoundary(Photon_Ptr)) {
        Hop(Photon_Ptr);    /* move to boundary plane. */
        CrossOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
    } else {
        Hop(Photon_Ptr);
        Drop(Photon_Ptr);

        switch (In_Ptr->TypeBias) {
            case 0:                        /* No importance sampling */
                Spin(regionspecs[Photon_Ptr->tetrahedron->region].g, Photon_Ptr);
                break;
            case 37:                    /* BOE2 Biasing by Lima et al. */
                SpinBias37(regionspecs[Photon_Ptr->tetrahedron->region].g, Photon_Ptr, PhotonCont_Ptr, In_Ptr);
        }
    }
}

void HopDropSpin(SimulationStruct *In_Ptr, TetrahedronPhotonStruct *Photon_Ptr, TetrahedronPhotonStruct *PhotonCont_Ptr,
                 OutStruct *Out_Ptr) {
    HopDropSpinInTissue(In_Ptr, Photon_Ptr, PhotonCont_Ptr, Out_Ptr);

    if (Photon_Ptr->w < WEIGHT && !Photon_Ptr->dead)
        Roulette(Photon_Ptr);

}

void CopyPhotonStruct(TetrahedronPhotonStruct *OrigPhoton_Ptr, TetrahedronPhotonStruct *DestPhoton_Ptr) {

    DestPhoton_Ptr->x = OrigPhoton_Ptr->x;
    DestPhoton_Ptr->y = OrigPhoton_Ptr->y;
    DestPhoton_Ptr->z = OrigPhoton_Ptr->z;    /* Cartesian coordinates.[cm] */
    DestPhoton_Ptr->ux = OrigPhoton_Ptr->ux;
    DestPhoton_Ptr->uy = OrigPhoton_Ptr->uy;
    DestPhoton_Ptr->uz = OrigPhoton_Ptr->uz;/* directional cosines of a photon. */
    DestPhoton_Ptr->w = OrigPhoton_Ptr->w;            /* weight. */
    DestPhoton_Ptr->dead = OrigPhoton_Ptr->dead;        /* 1 if photon is terminated. */

    DestPhoton_Ptr->MinCos = OrigPhoton_Ptr->MinCos;
    DestPhoton_Ptr->tetrahedron = OrigPhoton_Ptr->tetrahedron;
    DestPhoton_Ptr->NextTetrahedron = OrigPhoton_Ptr->NextTetrahedron;

    DestPhoton_Ptr->s = OrigPhoton_Ptr->s;            /* current step size. [cm]. */
    DestPhoton_Ptr->sleft = OrigPhoton_Ptr->sleft;        /* step size left. dimensionless [-]. */

    DestPhoton_Ptr->OpticalPath = OrigPhoton_Ptr->OpticalPath;
    DestPhoton_Ptr->MaxDepth = OrigPhoton_Ptr->MaxDepth;
    DestPhoton_Ptr->LikelihoodRatio = OrigPhoton_Ptr->LikelihoodRatio;
    DestPhoton_Ptr->FstBackReflectionFlag = OrigPhoton_Ptr->FstBackReflectionFlag;
    DestPhoton_Ptr->LocationFstBias = OrigPhoton_Ptr->LocationFstBias;
    DestPhoton_Ptr->NumBackwardsSpecularReflections = OrigPhoton_Ptr->NumBackwardsSpecularReflections;

}
