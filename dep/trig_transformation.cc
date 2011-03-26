
/*
  transform R^3 into the simplex defined by nucleotide composition space with
  A+C+G+T = 1 and all positive, using filipe ribeiro's transformation.
  
  c1 = cos(a)^2
  c2 = sin(a)^2 cos(b)^2
  c3 = sin(a)^2 sin(b)^2 cos(c)^2
  c4 = sin(a)^2 sin(b)^2 sin(c)^2   // 3 degrees of freedom: the angles 'a', 'b', and 'c'.  

*/

    void r3_to_composition_spherical(gsl_vector const* x, 
                                     double * comp, 
                                     double comp_gradient[4][3])
    {
        double sin[3], sin2[3], cos[3], cos2[3];

        double ran_offset = 0.0;
        //double ran_offset = M_PI / 4.0;

        for (size_t d = 0; d != 3; ++d)
        {
            sin[d] = gsl_sf_sin(gsl_vector_get(x, d) - ran_offset);
            cos[d] = gsl_sf_cos(gsl_vector_get(x, d) - ran_offset);
            sin2[d] = gsl_pow_2(sin[d]);
            cos2[d] = gsl_pow_2(cos[d]);
        }
    
        comp[0] = cos2[0];
        comp[1] = sin2[0] * cos2[1];
        comp[2] = sin2[0] * sin2[1] * cos2[2];
        comp[3] = sin2[0] * sin2[1] * sin2[2];

        if (! (normalized(comp, 4, 1e-10) &&
               all_positive(comp, 4)))
        {
            fprintf(stderr, "probability_transform: invalid composition\n");
            exit(2);
        }

        comp_gradient[0][0] = -2.0 * cos[0] * sin[0];
        comp_gradient[0][1] = 0.0;
        comp_gradient[0][2] = 0.0;
        
        comp_gradient[1][0] = 2.0 * sin[0] * cos[0] * cos2[1];
        comp_gradient[1][1] = -2.0 * sin[1] * cos[1] * sin2[0];
        comp_gradient[1][2] = 0.0;

        comp_gradient[2][0] = 2.0 * cos[0] * sin[0] * sin2[1] * cos2[2];
        comp_gradient[2][1] = 2.0 * cos[1] * sin[1] * sin2[0] * cos2[2];
        comp_gradient[2][2] = -2.0 * sin[2] * cos[2] * sin2[0] * sin2[1];

        comp_gradient[3][0] = 2.0 * cos[0] * sin[0] * sin2[1] * sin2[2];
        comp_gradient[3][1] = 2.0 * cos[1] * sin[1] * sin2[0] * sin2[2];
        comp_gradient[3][2] = 2.0 * cos[2] * sin[2] * sin2[0] * sin2[1];

    }


    //transform normalized 4D coordinates to spherical coordinates
    void composition_to_r3_spherical(double const* comp, double * sphere)
    {
        double sin2[3];
        
        sphere[0] = GSL_REAL(gsl_complex_arccos_real(sqrt(comp[0])));
        sin2[0] = gsl_pow_2(gsl_sf_sin(sphere[0]));

        sphere[1] = GSL_REAL(gsl_complex_arccos_real(sqrt(comp[1] / sin2[0])));
        sin2[1] = gsl_pow_2(gsl_sf_sin(sphere[1]));

        sphere[2] = GSL_REAL(gsl_complex_arccos_real(sqrt(comp[2] / (sin2[0] * sin2[1]))));

    }
