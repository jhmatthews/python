#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "atomic.h"
#include "python.h"


// define dot product function todo these will need to be new subroutines - some already exist in python!!

/*double dotProduct (double v1[], double v2[]) {
    double dotResult;
    dotResult = ((v1[0] * v2[0]) + (v1[1] * v2[1]) + (v1[2] * v2[2]));
    printf("Dot product from function call of the components is \"%f\".\n", dotResult);
    return (dotResult);
}
*/
// define normalize function
// todo split the magnitude part out of this function

int normalize (double v1[], double v2[]) {
    double fMag;

    fMag = sqrt( pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
    if (fMag == 0) {
        return 1;
    }
    else {
        v2[0] = v1[0] / fMag;
        v2[1] = v1[1] / fMag;
        v2[2] = v1[2] / fMag;
        return 0;
    }
}
/*
void crossProduct(double v1[], double v2[], double vR[]) {
//    double verbosity;
//    verbosity = 1;                                      // todo pull in the global verbosity flag

    if (verbosity == 10) {
        printf("Cross product function - first 4 components are \" %f %f %f %f \".\n", v1[1], v2[2], v1[2], v2[1]);
    }

    vR[0] =   ( (v1[1] * v2[2]) - (v1[2] * v2[1]) );

    if (verbosity == 10) {
        printf("Cross product function - first result is %f \".\n", vR[0]);
        printf("Cross product function - second 4 components are \" %f %f %f %f \".\n", v1[0], v2[2], v1[2], v2[0]);
    }

    vR[1] = - ( (v1[0] * v2[2]) - (v1[2] * v2[0]) );

    if (verbosity == 10) {
        printf("Cross product function - second result is %f \".\n", vR[1]);
        printf("Cross product function - third 4 components are \" %f %f %f %f \".\n", v1[0], v2[1], v1[1], v2[0]);
    }

    vR[2] =   ( (v1[0] * v2[1]) - (v1[1] * v2[0]) );

    if (verbosity == 10) {
        printf("Cross product function - third result is %f \".\n", vR[2]);
    }
}
*/
double magnitude (double v1[]) {
    double Mag;
    Mag = sqrt( pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
//    printf("Magnitude is %f \".\n", Mag);
    return Mag;
}


int
poltest (p, pp)
    PhotPtr p;                              // the original wind photon structure
    PhotPtr pp;                             // the new extracted photon structure
{
    double *vect1;
    double *vect2;
    double adotb3;
//    double *new_lmn;
    vect1 = p->orig;                         // vect1 now points to the address of the inbound photon
    vect2 = pp->lmn;                        // vect2 now points to the address of the outbound photon
//    new_lmn = p->lmn;                         // orig maybe points to the original photon direction
//    Log ("POLTEST p %p L %f M %f N %f \n", p, p->lmn[0], p->lmn[1], p->lmn[2]);
//    Log ("POLTEST pp %p, L %f M %f N %f  \n", pp, pp->lmn[0], pp->lmn[1], pp->lmn[2]);
//    Log ("POLTEST P %p Lorig %f Morig %f Norig %f L %f M %f N %f \n", p, vect1[0], vect1[1], vect1[2], vect2[0], vect2[1], vect2[2]);
    int i, j, mainloop, r1, c1, r2, c2, k, errno;
//    int count, n,
//    int verbosity;                                        // there is a global verbosity flag in python.h
//    int rand(void);                                       // not needed here
    double stokes_parameter_q_dimensionless, stokes_parameter_u_dimensionless, stokes_parameter_v_dimensionless;
    double thomson_scattering_angle, polarization_degree, intensity, diff;
    //double polarization_angle_chi;
    //double array_double[6];
    int loop;
    int ref_axis;                           // reference axis = 0 if z axis, 1 if x axis
//    double adotb2, adotb3;                   // dotproducts
    double prod1, prod2;
//    double prod3;
    double costheta, Iout;
//    double dotResult;
    double magnitude_a,magnitude_b;
    double cosphi, sinphi;
    double phi_dot_result;
    double magnitude_r_obs;                  // magnitude of r axis in observer Frame of Reference
    double magnitude_r_out;                  // magnitude of r axis in scattering Frame of Reference
    double mag_mag;                          //
    double cross_mag;                        //
    double StokesVector[4][4] = {};          // Stokes vector
    double RotationVector[3][3] = {};        // Rotation vector
    double StokesOut[4][1] = {};             // Stokes matrix scattered photon
    double StokesObs[3][1] = {};             // Stokes matrix scattered photon observer frame of reference
    double StokesIn[4][1] = {};              // Stokes Matrix for incoming photon
//    double vect1[3];                         // n_in scattering plane vector
//    double vect2[3];                         // n_out scattering plane vector
    double vectl[3];                         // l_scattering plane un-normalized vector
    double vectResult[3];                    // generic cross product result vector
    double vectrObsResult[3];                // unnormalized r_observer frame vector
    double vectlObsResult[3];                // unnormalized l_observer frame vector
    double vectNormal[3];                    // normalized r_scattering plane vector
    double vectlNormal[3];                   // normalized l_scattering plane vector
    double vectrObsNormal[3];                // normalized r_observer frame vector
    double vectlObsNormal[3];                // normalized l_observer frame vector*/
    double zaxis[3] = {0.0, 0.0, 1.0};       // z axis is the reference axis for our observer
    double xaxis[3] = {1.0, 0.0, 0.0};       // x axis is the alternative reference axis for our observer
    double vectxprod[3];                     // cross product r_obs and r_out to calculate sinphi
    double sinphi_numerator;                 // (r_obs x r_s).n_obs

    /* linearly polarized light is completely defined by the Stokes parameters Q and U (V=0) with the polarization
     * angle chi given by tan 2chi = U/Q and the angle gamma given by sin2gamma = V/I - but gamma is zero for
     * linearly polarized light */

    stokes_parameter_q_dimensionless = 0;    // Initially the photon is unpolarized
    stokes_parameter_u_dimensionless = 0;    // Initially the photon is unpolarized
    stokes_parameter_v_dimensionless = 0;    // Initially the photon is unpolarized
    thomson_scattering_angle = 0;            // Initialize Thomson scattering angle
    //polarization_angle_chi = 0;            // Initialize polarization angle chi  - it should lie between 0 and pi
    polarization_degree = 0;                 // Initially unpolarized degree
//    adotb = 0.0;                             // Initialize dot product of input position and output position vectors
//    adotb2 = 0.0;                            // DEBUG: Second way to calculate dot product
    adotb3 = 0.0;                            // DEBUG: Third way to calculate dot product
    prod1 = 0;                               // Initialize dot product of x input position and x output position
    prod2 = 0;                               // Initialize dot product of y input position and y output position
    magnitude_a = 0;                         // Initialize magnitude of input vector
    magnitude_b = 0;                         // Initialize magnitude of output vector
    costheta = 0;                            // Initialize angle between the vectors
    cosphi = 0;                              // Initialize angle between scattered frame and observer frame of reference
    sinphi = 0;                              // Initilize
//   n = 6;                                   // Initialize random number counter
//    intensity = 0;                           // Initialize the wave intensity
    intensity = p->w;                        // set intensity to photon weight
    loop = 100;                              // Loop counter for how many times to repeat the calculations
    // todo loop could be an input parameter
    mainloop = 0;                            // main program loop counter
    diff = 0;                                // difference between cos theta and 1
    errno = 1;                               // error detector
    Iout = 0;                                // Intensity of outgoing photon
    phi_dot_result = 0.0;


    /* Input position parameters will be in 3D cartesian coordinates on an x,y,z axis.
     * The first 3 coordinates entered represent the initial photon direction vector components.
     * The last 3 coordinates entered represent the final photon direction vector components. */


//     verbosity = 10;                          // verbosity indicator. 10 means print every debug line to screen
//
//    /* Output will be written to file called write.txt
//
//
//    FILE *fp;
//    fp = fopen ("write.txt","w");
//    if (fp == NULL) {
//        printf ("File not created properly, errno = %d\n", errno);
//        return 1;
//    }
//
//    // Create header of output file
//
//    fprintf (fp, "n_i_in          n_j_in     n_k_in       n_i_out      n_j_out      n_k_out       "
//            "INTENSITY  r_i_s   r_j_s   r_k_s   q        u       v      "
//            "ctheta angle     P      l_i_s    l_j_s   l_k_s  r_i_o   r_j_o    r_k_o   l_i_o   l_j_o  l_k_o"
//            "\n\n");
//    //todo add observer q and u and decide what we want to log and where- we don't need all this info
//    //todo decide what we want for diagnostics
//
//
//    // It is possible to enter the 6 parameters on the command line  - or they can be randomly generated
//
//    //if (verbosity == 10) { printf("This program was called with \"%s\".\n", argv[0]); }
//    // Initialize random number generator
//
//    //     srand((signed) time(NULL));
//
//    // Generate and Print 6 random numbers
//
//    /* NOT NEEDED RANDOM NUMBER GENERATION ROUTINE
//     *
//    for (mainloop = 0; mainloop < loop; mainloop++) {
//        if (argc <= 1) {
//
//            /*for (i = 0; i < n; i++) {
//                array_double[i] = (RAND_MAX - (2 * rand()));  //This generates positive and negative random numbers
//                printf("array_double[%d] is initialized with %f.\n", i, array_double[i]);
//            }
//
//
//            // First calculate the dot product
//
//            prod1 = (array_double[0] * array_double[3]);
//            prod2 = (array_double[1] * array_double[4]);
//            prod3 = (array_double[2] * array_double[5]);
//
//            vect1[0] = array_double[0];
//            vect1[1] = array_double[1];
//            vect1[2] = array_double[2];
//
//            vect2[0] = array_double[3];
//            vect2[1] = array_double[4];
//            vect2[2] = array_double[5];
//
//            if (verbosity == 10)
//            {
//                printf("Product of the x components is \"%f\".\n", prod1);
//                printf("Product of the y components is \"%f\".\n", prod2);
//                printf("Product of the z components is \"%f\".\n", prod3);
//            }
//
//            adotb = (prod1 + prod2 + prod3);                    // this calculates the dot product from products
//
//            if (verbosity == 10) {
//
//                printf("Dot product from prod1 plus prod2 plus prod3 of the components is \"%f\".\n", adotb);
//            }
//
//            //Use function call for dot product
//
//            dotResult = 0;                                      // initialize dot product result as we are in a loop
//
//            adotb3 = (dotProduct(vect1, vect2));                // dotproduct result from input and output vectors
//
//            if (verbosity == 10) {
//                printf("Dot product assigned as a result of calling dotProduct function call is \"%f\".\n", adotb3);
//            }
//*/
//            // Now calculate the magnitude of input and output vectors
//
//            /*if (magnitude(vect1) == 0) {
//                printf("ERROR in magnitude a function\n");      // magnitude is zero - quit loop unable to continue
//                return 1;
//            } else {
//                printf("Magnitude a is %f\n", Mag);
//                magnitude_a = Mag;                              // magnitude of input vector
//                printf("Magnitude a is %f\n", magnitude_a);
//            };
//             */
//
// /*           // magnitude of input vector
//            magnitude_a = magnitude(vect1);                     // magnitude of input vector
//
//            if (magnitude_a == 0) {
//               printf("ERROR in magnitude a function\n");      // magnitude is zero - quit loop unable to continue
//               return 1;                                       // TODO what error codes should I be using?check python
//            }
//
//            if (verbosity == 10) {
//                printf("Magnitude of input vector is \"%f\".\n", magnitude_a);
//            }
//
//
//            magnitude_b = magnitude(vect2);                     //  magnitude of output vector
//
//            if (magnitude_b == 0) {
//                printf("ERROR in magnitude b function\n");
//                return 1;                                   // TODO what error codes should I be using?check python
//            }
//
//            if (verbosity == 10) { printf("Magnitude of output vector is \"%f\".\n", magnitude_b); }
//        }
//
//        // end of no command line input parameters loop i.e. random number generator loop
//
//        */
//        //*******************************COMMAND LINE INPUT****************************************************
//  /*      if (argc > 1) {
//            // Read the input parameters from the COMMAND LINE and store in an array called argv with a pointer
//            // of the "count" of number of elements  - we really only need 6 elements
//
//            for (count = 1; count < argc; count++) {
//                printf("argv[%d] = %s\n", count, argv[count]);
//            }
//
//            // todo this could be an input parameter
//
//            mainloop = 100;                             // no. of times to create a random set of input parameters
//
//            // Now we have captured our initial and final coordinates.
//            // To calculate the cosine of the angle between these 2 points. We use the following formulae:
//            // a · b = ax × bx + ay × by + az × bz
//            // a · b = |a| × |b| × cos(θ)
//
//            // First calculate the dot product
//            //TODO: This is basically a repetition of the block above with the counter number increased
//            //TODO - this could be better coded
//
//            if (verbosity == 10) {
//                printf("Initial x component is \"%s\".\n", argv[1]);
//                printf("Final x component is \"%s\".\n", argv[4]);
//            }
//
//            prod1 = (atof(argv[1]) * atof(argv[4]));
//
//            if (verbosity == 10) {
//                printf("Product of the x components is \"%f\".\n", prod1);
//                printf("Initial y component is \"%s\".\n", argv[2]);
//                printf("Final y component is \"%s\".\n", argv[5]);
//            }
//
//            prod2 = (atof(argv[2]) * atof(argv[5]));
//
//            if (verbosity == 10) {
//                printf("Product of the y components is \"%f\".\n", prod2);
//
//                printf("Initial z component is \"%s\".\n", argv[3]);
//                printf("Final z component is \"%s\".\n", argv[6]);
//            }
//
//            prod3 = (atof(argv[3]) * atof(argv[6]));
//
//            if (verbosity == 10) {
//                printf("Product of the z components is \"%f\".\n", prod3);
//            }
//
//            adotb = (prod1 + prod2 + prod3);
//
//            if (verbosity == 10) {
//                printf("Dot product1 of the components is \"%f\".\n", adotb);
//            }
//            // put terminal input into 2 vectors for later calculations in functions
//            // todo if this was moved up to the top then the rest of the code could be equivalent with the random loop
//
//            vect1[0] = atof(argv[1]);                        //TODO: this needs tested - star or not?
//            vect1[1] = atof(argv[2]);
//            vect1[2] = atof(argv[3]);
//            vect2[0] = atof(argv[4]);
//            vect2[1] = atof(argv[5]);
//            vect2[2] = atof(argv[6]);
//
//            if (verbosity == 10) {
//
//                printf("component 1 of vector 1: \"%f\".\n", vect1[0]);         //TODO: test
//                printf("component 2 of vector 1: \"%f\".\n", vect1[1]);
//                printf("component 3 of vector 1: \"%f\".\n", vect1[2]);
//                printf("component 1 of vector 2: \"%f\".\n", vect2[0]);
//                printf("component 2 of vector 2: \"%f\".\n", vect2[1]);
//                printf("component 3 of vector 2: \"%f\".\n", vect2[2]);
//            }
//
//            // Use function call to get dotproduct

    adotb3 = (dot(vect1, vect2));

// magnitude of input vector
    magnitude_a = magnitude(vect1);                     // magnitude of input vector
    if (magnitude_a == 0)
    {
      Error ("poltest: Input vector magnitude is zero \n");
      return (1);
    }
//    printf("Magnitude of input vector is \"%f\".\n", magnitude_a);


//magnitude of output vector
    magnitude_b = magnitude(vect2);                     //  magnitude of output vector
    if (magnitude_b == 0)
    {
        Error ("poltest: Error in magnitude of output vector \n");
        return (1);
    }
//    printf("Magnitude of output vector is \"%f\".\n", magnitude_b);

//
//            if (verbosity == 10) { printf("Dot product3 of the components is \"%f\".\n", adotb3); }
//
//
// Now calculate the magnitude of each vector TODO: I think we can use the function to do this too
//
//            magnitude_a = sqrt((atof(argv[1]) * atof(argv[1])) + (atof(argv[2]) * atof(argv[2])) +
//                               (atof(argv[3]) * atof(argv[3])));
//
//            if (verbosity == 10) { printf("Magnitude of input vector is \"%f\".\n", magnitude_a); }
//
//            magnitude_b = sqrt((atof(argv[4]) * atof(argv[4])) + (atof(argv[5]) * atof(argv[5])) +
//                               (atof(argv[6]) * atof(argv[6])));
//
//            if (verbosity == 10) { printf("Magnitude of output vector is \"%f\".\n", magnitude_b);}
//        }
//        // end of CLI input loop
//
        //  *****************Code from here onwards is common to both CLI and random generation of input parameters

        if (magnitude_a && magnitude_b) {

            //************************* LOCAL FRAME OF REFERENCE OF PHOTON **************************************
            // We have 2 non zero magnitudes so we can proceed to calculate our triad of axes and scattering angle
            // First calculate r scattering axis: normalized cross product between incoming vector and outgoing
            // vector
            //**************************** r axis local ********************************************************

            int m,n;
            ref_axis = 0;                                                   // initialize reference axis as 0 = z-axis

            n = memcmp( vect1, vect2, sizeof(*vect1) );

            if (n == 0) {

                Log ("Poltest: Input vector is the same as output vector.\n");    // input and output vectors are the same

                m = memcmp(vect1, zaxis, sizeof(*vect1));                   // check if they are the z axis

                if (m == 0) {

                    Log ("Poltest: Both vectors are z axis.\n");                  // input and output vectors are the z axis

                    vectResult[0] = xaxis[0];                              // set outgoing scattering vector as x axis
                    vectResult[1] = xaxis[1];                              // set outgoing scattering vector as x axis
                    vectResult[2] = xaxis[2];                              // set outgoing scattering vector as x axis
                    ref_axis = 1;                                          // reference axis flag is set to x

                    Log ("Poltest: Set reference axis as x axis.\n");              // input and output vectors are the x axis

                } else {

                    vectResult[0] = zaxis[0];                               // set outgoing scattering vector as z axis
                    vectResult[1] = zaxis[1];                               // set outgoing scattering vector as z axis
                    vectResult[2] = zaxis[2];                               // set outgoing scattering vector as z axis

                    Log ("Poltest: Set reference axis as z axis.\n");              // input and output vectors are the z axis
                }

            } else {

                cross(vect1, vect2, vectResult);

//            if (verbosity == 10) {
//
//                printf("r scattering axis (unnormalized): Cross product of vector in and vector out x:"
//                               " %.4f, y: %.4f, z: %.4f\n", vectResult[0], vectResult[1], vectResult[2]);
//            }
            }
            // now normalize the r axis in the scattering frame - new vector is called vectNormal

            if (normalize(vectResult,vectNormal) != 0)
            {
              // TODO: This is a very rare case and unlikely to happen in practice in the main code
              // but for completeness, this means that the input vector and the output vector are in the same
              // scattering plane.
              // - if r_s = r_obs i.e. the dotproduct of vect1 and vect2 = 0
              // then pick any vector perpendicular to vect2 e.g. parallel to or equal to r_obs.
              // Then print a warning to screen

              Error ("poltest: Error in r scattering normalize function\n");
              Log ("Poltest ERROR Input and output vectors are anti-parallel\n");
              Log ("Poltest ERROR choose any perpendicular vector to output vector to be the r scattering axis\n");
              // TODO add this into the code or just skip this particular photon**************************!!!!!!!
              return (1);
            }

            // print out this r axis information

//            printf("r scattering axis (normalized): x: %.4f, y: %.4f, z: %.4f\n",
//                   vectNormal[0], vectNormal[1], vectNormal[2]);

            //**************************** l axis local ********************************************************

            // Now we have r we can calculate l in the scattering plane
            // Calculate l scattering axis: normalized cross product between outgoing unit vector and r axis vector

            // todo is this the correct way round for this cross product??
            // todo use python cross product routine

            cross(vect2, vectResult, vectl);

//            printf("l scattering axis (unnormalized): Cross product of r vector and vector out x: "
//                           "%.4f, y: %.4f, z: %.4f\n", vectl[0], vectl[1], vectl[2]);

            // now normalize the l axis in the scattering frame - new vector is called vectlNormal
            //todo the vector names aren't consistent - better to keep a theme going here

            if (normalize(vectl,vectlNormal) != 0)
            {
                Error ("poltest: Error in l scattering normalize function\n");
                return (1);
            }

            // print out this l axis information

//            printf("l scattering axis (normalized): x: %.4f, y: %.4f, z: %.4f\n",
//                   vectlNormal[0], vectlNormal[1], vectlNormal[2]);


            //***********************OBSERVER FRAME OF REFERENCE***********************************************
            // Calculate r observer frame: normalized cross product between outgoing n vector and (0,0,1) ie z axis
            //**************************** r axis observer ********************************************************
            //todo add error checking

            if (ref_axis == 0) {
                cross(vect2, zaxis, vectrObsResult);
            } else {
                cross(vect2, xaxis, vectrObsResult);
            }

            // todo use log command
//            if (verbosity == 10) {
//                printf("r observer frame (unnormalized): Cross product of photon outgoing vector and z axis x:"
//                               " %.4f, y: %.4f, z: %.4f\n", vectrObsResult[0], vectrObsResult[1], vectrObsResult[2]);
//            }

            // todo use python normalize command
            if (normalize(vectrObsResult,vectrObsNormal) != 0) {
                Log ("ERROR in r observer frame normalize function\n");
                Error ("poltest:ERROR in r observer frame normalize function \n");
              return (1);                                   // TODO add error check on poltest return
            }

            //todo use log command in python

//            printf("r observer frame (normalized): x: %.4f, y: %.4f, z: %.4f\n",
//                   vectrObsNormal[0], vectrObsNormal[1], vectrObsNormal[2]);

            //**************************** l axis observer ********************************************************
            // Calculate l observer axis: normalized cross product between outgoing unit vector and r observer axis

            //todo add error checking

            cross(vect2, vectrObsNormal, vectlObsResult);

            // todo use python log command

//            printf("l observer frame (unnormalized): Cross product of vector out and r observer vector x: "
//                         "%.4f, y: %.4f, z: %.4f\n", vectlObsResult[0], vectlObsResult[1], vectlObsResult[2]);

            if (normalize(vectlObsResult,vectlObsNormal) != 0)
            {
                Error ("poltest: in l observer frame  normalize function\n");
                return (1);
            }

//            printf("l observer frame (normalized): x: %.4f, y: %.4f, z: %.4f\n",
//                   vectlObsNormal[0], vectlObsNormal[1], vectlObsNormal[2]);

            //**************************SCATTERING ANGLE*****************************************************
            /* The cosine of the scattering angle between the input and output vectors
             * is equivalent to the dot product divided by the
             * multiplication of the input and output vector magnitudes: a · b = |a| × |b| × cos(θ) */

            // todo add error checking

            costheta = adotb3 / (magnitude_a * magnitude_b);

//            printf("Cosine of angle between the input and output n scattering vectors: \"%f\".\n", costheta);

            // If the beam is passing straight through then angle = zero hence no point converting to degrees

            // TODO: check that this case is correct else fix the code so that costheta != 1.0

            diff = fabs(costheta) - 1.0000;

//            if (verbosity == 10) {
//                printf("the difference is \"%f\".\n", diff);
//            }

            if (fabs(diff) > 10e-8) {
                thomson_scattering_angle = (acos(costheta) * (180 / 3.14159265));  // convert angle to degrees
            }

//            if (verbosity == 10) { printf("Thomson scattering angle: \"%f\" degrees.\n", thomson_scattering_angle);}

            //**************************POLARIZATION DEGREE*****************************************************
            /* Now we have calculated the direction of propagation of the photon
             * and 2 vectors: inwards and outwards
             * Polarization in the plane can be calculated */

            // todo add error checking

            polarization_degree = (1 - (costheta * costheta)) / (1 + (costheta * costheta));


            // TODO: make this a for loop or a switch?


/*            if (polarization_degree <= 0.00001) {
                printf("Stream of photons is unpolarized.\n");
            }

            if (polarization_degree == 1.0) {
                printf("Stream of photons is completely polarized.\n");
            } else if (polarization_degree > 0.00001 && polarization_degree < 1.0) {
                printf("Stream of photons is partially polarized.\n");
            }
*/
//            printf("Polarization degree: \"%f\".\n", polarization_degree);

            //**************************INTENSITY*****************************************************

            // todo add a positive limit to the intensity random number generator

            // intensity = rand();                // Random intensity is chosen. It will be irrelevant anyway!
//            printf("POLTEST p weight or intensity: \"%f\".\n", intensity);



            //**************************STOKES VECTOR SETUP*****************************************************
            /* set up stokes matrix for unpolarized light*/
            /* Stokes Vector = [I, Q, U, V] nb we are not using v as it is for circularly polarized light only
             * Stokes Vector (Dimensionless) = Stokes Vector / Intensity
             * [1,q,u,v] = [1,0,0,0] for unpolarized light */
            //**************************INPUT STOKES VECTOR*****************************************************

            StokesIn[0][0] = intensity;

            // set up stokes vector

            StokesVector[0][0] = (0.75 * ((costheta * costheta) + 1));

            StokesVector[0][1] = (0.75 * ((costheta * costheta) - 1));

            StokesVector[1][0] = (0.75 * ((costheta * costheta) - 1));

            StokesVector[1][1] = (0.75 * ((costheta * costheta) + 1));

            StokesVector[2][2] = (1.5 * costheta);

            StokesVector[3][3] = (1.5 * costheta);


            // Displaying the stokes vector

            // Set up some counters for the loops

            // element counters i j and k to go through the matrices

            i = 0;
            j = 0;
            k = 0;

            // row and column indicators for input and output matrices
            // todo calculate this rather than hard code it?

            r1 = 4;
            c1 = 4;
            r2 = 4;
            c2 = 1;

/*            printf("\nStokes vector:\n");

            for (i = 0; i < 4; ++i) {
                for (j = 0; j < 4; ++j) {
                    printf("%f  ", StokesVector[i][j]);
                    if (j == 3) {
                        printf("\n\n");
                    }
                }
            }
*/

            //**************************OUTPUT STOKES VECTOR*****************************************************
            // Initializing all elements of StokesOut matrix to 0
            // todo is this even needed?

            for (i = 0; i < 4; ++i) {
                    StokesOut[i][0] = 0;
            }

            // Multiplying matrices StokesArray and StokesIn and storing result in StokesOut matrix

            for (i = 0; i < r1; ++i) {
                for (j = 0; j < c2; ++j) {
                    for (k = 0; k < c1; ++k) {
                        StokesOut[i][j] += StokesVector[i][k] * StokesIn[k][j];
                    }
                }
            }

            // Displaying the result

 //           printf("\nStokes Output Matrix:\n");

//            for (i = 0; i < 4; ++i)
//            {
//                    printf("%f  \n\n", StokesOut[i][0]);
//            }

            //**************************DIMENSIONLESS OUTPUT STOKES VECTOR*****************************************
            // Displaying the dimensionless stokes vector result

//            printf("\nDimensionless Stokes Output Matrix:\n");
            Iout = StokesOut[0][0];
            for (i = 0; i < 4; ++i) {
                    StokesOut[i][0] = (StokesOut[i][0] / Iout);
                    // printf("%f  \n\n", StokesOut[i][0]);
            }

            //**************************GET q,u,v *****************************************************************
            stokes_parameter_q_dimensionless = StokesOut[1][0];
            stokes_parameter_u_dimensionless = StokesOut[2][0];
            stokes_parameter_v_dimensionless = StokesOut[3][0];     // we don't use v here

//            printf("Final photon q stoke: \"%f\".\n", stokes_parameter_q_dimensionless);
//            printf("Final photon u stoke: \"%f\".\n", stokes_parameter_u_dimensionless);
            // printf("Final photon v stoke: \"%f\".\n", stokes_parameter_v_dimensionless);

            //**************************WRITE FILE OUTPUT LINE*****************************************
            // write all values calculated so far to file

/*            printf("%12.0f %12.0f %12.0f "
                            "%12.0f %12.0f %12.0f %12.0f"
                            "%7.3f %7.3f %7.3f "
                            "%7.3f %7.3f %7.3f "
                            "%7.3f %7.3f %7.3f "
                            "%7.3f %7.3f %7.3f "
                            "%7.3f %7.3f %7.3f "
                            "%7.3f %7.3f %7.3f \n",
                    vect1[0], vect1[1], vect1[2],
                    vect2[0], vect2[1], vect2[2], intensity,
                    vectNormal[0], vectNormal[1], vectNormal[2],
                    stokes_parameter_q_dimensionless, stokes_parameter_u_dimensionless,
                    stokes_parameter_v_dimensionless,
                    costheta, thomson_scattering_angle, polarization_degree,
                    vectlNormal[0], vectlNormal[1], vectlNormal[2],
                    vectrObsNormal[0], vectrObsNormal[1], vectrObsNormal[2],
                    vectlObsNormal[0], vectlObsNormal[1], vectlObsNormal[2]);
*/

            //**************************SANITY CHECK POLARIZATION DEGREE *****************************************
            // Calculate the polarization degree here - (a double check and should be equal to value calculated
            // above)

            polarization_degree = sqrt((stokes_parameter_q_dimensionless*stokes_parameter_q_dimensionless)
                                           + (stokes_parameter_u_dimensionless*stokes_parameter_u_dimensionless));
            //fprintf(fp, "%7.3f\n", polarization_degree);


            //**************************PREPARE STOKES ROTATION VECTOR *****************************************


            //*****************CALCULATE PHI - THE CLOCKWISE ROTATION ANGLE BETWEEN SCATTER *********************
            // ************************PLANE AND MERIDIAN AND THE OBSERVER PLANE ********************************

            //cos phi is the dot product of r_obs and r_out divided by the magnitude of the 2 vectors

            //calculate magnitude_r_obs
            // todo logically this should come up above with the r obs calculations
            magnitude_r_obs = magnitude(vectrObsResult);

            if (magnitude_r_obs == 0)
            {
              Error ("poltest: Null magnitude for r vector in observer frame of reference \n");
              return (1);
            }

 //           printf("Magnitude of r vector in observer frame of reference is %f \n",magnitude_r_obs);



            //calculate magnitude_r_out
            magnitude_r_out = magnitude(vectResult);

            if (magnitude_r_out == 0)
            {
                Error ("poltest: Null magnitude for r vector out in photon frame of reference \n");
                return (1);
            }

//            printf("Magnitude of r outgoing vector in scattering frame of reference is %f \n",magnitude_r_out);



            //calculate r_obs_dot_r_out
            // todo I think the normalised vectors would be fine here instead of dividing by the magnitudes

            phi_dot_result = dot(vectResult, vectrObsResult);

            //todo add verbosity flag on printfs

//            printf("Phi dot result is %f \n",phi_dot_result);

            mag_mag = magnitude_r_obs * magnitude_r_out;                    //the two magnitudes are already nonzero

//            printf("Product of magnitude of the observer by the outgoing magnitude is: %f \n",mag_mag);

            cosphi = (phi_dot_result / mag_mag);

//            printf("Cos phi is: %f \n",cosphi);


            // calculate r_obs_cross_r_out

            cross(vectrObsResult, vectResult, vectxprod);

            // take scalar prod of vectxprod and direction to observer as both should be perp to dir of obs. par
            // or anti parallel

            cross_mag = magnitude(vectxprod);

            sinphi_numerator = dot(vectxprod,vect2);

            //todo vect 2 is not a unit vector so it is too big (but the sign is correct) calculate it!!
/*
            if (cross_mag == 0) {

                printf("ERROR null magnitude in cross product between vector r observer FoR and "
                               "vector r scattering frame of reference\n");
                return (1);                             // TODO what error codes should I be using?check python
            }

            printf("Magnitude of the cross product between vector r observer frame of reference and "
                           "vector r scattering frame of reference is: %f\n",cross_mag);
            */
            // sinphi = (cross_mag / mag_mag);

            sinphi = (sinphi_numerator / (mag_mag*magnitude_b));

//            printf("Sin phi and  cos phi are: %f %f \n",sinphi, cosphi);


            // todo initialize rotation vector

            // set the Rotation vector

            RotationVector[0][0] = (1);

            RotationVector[1][1] = ((2 * (cosphi * cosphi)) - 1);

            RotationVector[1][2] = (2 * (sinphi * cosphi));

            RotationVector[2][1] = (-2 * sinphi * cosphi);

            RotationVector[2][2] = ((2 * cosphi * cosphi) - 1);

            // Displaying the Rotation vector

            i = 0;
            j = 0;
            k = 0;

            r1 = 3;
            c1 = 3;
            r2 = 3;
            c2 = 1;
/*
            printf("\nRotation vector:\n");

            for (i = 0; i < r1; ++i) {
                for (j = 0; j < c1; ++j) {
                    printf("%f  ", RotationVector[i][j]);
                    if (j == c1 - 1) {
                        printf("\n\n");
                    }
                }
            }
*/
            // Initializing all elements of StokesObs matrix to 0
            // todo is this even needed?

            for (i = 0; i < 4; ++i) {
                    StokesObs[i][0] = 0;
            }

            // Displaying the result
/*
            printf("\nInitial Stokes Observer Matrix:\n");


            for (i = 0; i < 3; ++i) {
                    printf("%f  \n\n", StokesObs[i][0]);
            }
*/
            // todo - the matrix should be I Q U not i q u v as it is just now!!! I needs to be 1 I think.
            // (i)(1)
            // (q)(-1)
            // (u)(0)

            i = 0;
            j = 0;
            k = 0;

            r1 = 3;
            c1 = 3;
            r2 = 3;
            c2 = 1;

            for (i = 0; i < r1; ++i) {
                for (j = 0; j < c2; ++j) {
                    for (k = 0; k < c1; ++k) {
                        StokesObs[i][j] += RotationVector[i][k] * StokesOut[k][j];
                    }
                }
            }
            // Displaying the result
/*
            printf("\nFinal Stokes Observer Matrix:\n");


            for (i = 0; i < 3; ++i) {
                printf("%f  \n\n", StokesObs[i][0]);
            }
            */
            pp->q = StokesObs[1][0];
            pp->u = StokesObs[2][0];
 //           printf("POLTEST: Observer q is %f  \n\n", pp->q);
 //           printf("POLTEST: Observer u is %f  \n\n", pp->u);
            if (pp->q < -1 || pp->q > 1) {
                Error("Poltest: Abnormal Stokes parameter q is %f \n", pp->q);
            }
            if (pp->u < -1 || pp->u > 1) {
                Error("Poltest: Abnormal Stokes parameter u is %f \n", pp->u);
            }
        } else {
            Error ("Poltest: A vector magnitude is zero. Unable to calculate Thomson scattering angle "
                           "or Polarization degree\n");
        }

    // end of main loop
//    fclose (fp);
//    if (verbosity == 10) {printf ("File created successfully\n");}
    return (0);
}