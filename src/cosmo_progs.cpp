#include "cosmo_progs.h"

namespace fast {
    
    /* dynamical time at z in seconds */
    float t_dynamical(float z){
        return sqrt(3.*M_PI/16.0) / sqrt(G * Deltac_nonlinear(z) * rho_critz(z) * Msun / (CMperMPC*CMperMPC*CMperMPC));
    }
    
    /* helper function defined in Barkana 2002 */
    double barkana_I (float x){
        return ( pow(x,4.5)/(1-x) ) + ( (9.0/7.0)*pow(x,3.5) ) + ( (9.0/5.0)*pow(x,2.5) ) + ( 3*pow(x,1.5) ) +
        ( 9*pow(x,0.5) ) - ( 4.5*log((1+sqrt(x))/(1-sqrt(x))) );
    }
    
    /*
     analytic Lya damping wing optical depth contributed by a neutral IGM from zstart to zend (zstart > zend) at observed redshift z
     
     from Barkana 2002
     */
    double anal_taudamp(float z, float zstart, float zend){
        float taudampcoeff;
        
        if (z <= zstart){
            fprintf(stderr, "cosmo_progs, anal_taudamp: Warning, using damping wing in resonance.  Returning arbitrary high tau, 1e5\n");
            return 1e5;
        }
        
        if (!(zstart > zend)){
            fprintf(stderr, "cosmo_progs, anal_taudamp: ERROR. zstart must be > zend!\nReturning 0\n");
            return 0;
        }
        
        taudampcoeff = 0.00232 * (OMb*hlittle/0.0325) * pow(OMm/0.35, -0.5);
        
        return taudampcoeff * pow((1+z)/7.0, 1.5) *
        ( barkana_I((1+zstart)/(1+z)) - barkana_I((1+zend)/(1+z)) );
    }
    
    /* returns the hubble "constant" (in 1/sec) at z */
    double hubble(float z){
        return Ho*sqrt(OMm*pow(1+z,3) + OMr*pow(1+z,4) + OMl);
    }
    
    
    /* returns hubble time (in sec), t_h = 1/H */
    double t_hubble(float z){
        return 1.0/hubble(z);
    }
    
    /* returns the proper mean baryonic density at z in g/cm^3 */
    double mean_rho(float z){
        return OMb*(3.0*Ho*Ho / (8.0*M_PI*G)) * pow(1+z, 3);
    }
    
    /* returns the gas temperature assuming adiabatic cooling after instantaneous CMB decoupling */
    double T_gas_adiabatic(float z){
        float zdec = z_gasCMBdecoupling();
        return T_cmb*(1+zdec)*pow( (1+z)/(1+zdec), 2);
    }
    
    
    /* returns the cosmological jeans mass (in solar masses) at z, density Del=rho/<rho>, temperature T, and \mu
     Taken from Barkana & Loeb (2001), eq. 52
     */
    double M_jeans_general(float z, float Del, float T, float mu){
        float rho = mean_rho(z)*Del;
        return 2.9/sqrt(rho)*pow(k_B*T/G/mu/m_p, 1.5)/Msun;
    }
    
    
    /* returns the cosmological jeans mass (in solar masses) at z, before the region was ionized */
    double M_jeans_neutralIGM(float z){
        if (z< z_gasCMBdecoupling())
            return 5730/sqrt(OMm*hlittle*hlittle/0.15)*pow(OMb*hlittle*hlittle/0.022, -0.6)*pow((1+z)/10.0, 1.5);
        return 1.35e5 /sqrt(OMm*hlittle*hlittle/0.15);
    }
    
    
    /* the redshift at which the gas temperature diverges from the CMB temperature */
    float z_gasCMBdecoupling(){
        return 137*pow(OMb*hlittle*hlittle/0.022, 0.4)-1;
    }
    
    
    /*
     Function NEUTRAL_FRACTION returns the hydrogen neutral fraction, chi, given:
     hydrogen density (cm^-3)
     ionization rate (s^-1)
     */
    double neutral_fraction(double density, double gamma_ion){
        double chi, b;
        
        // approximation chi << 1
        chi = density * alphaB_10k / gamma_ion;
        if (chi < TINY){ return 0;}
        if (chi < 1e-5)
            return chi;
        
        //  this code, while mathematically accurate, is numerically buggy for very small x_HI, so i will use valid approximation x_HI <<1 above when x_HI < 1e-5, and this otherwise... the two converge seemlessly
        //get solutions of quadratic of chi (neutral fraction)
        b = -2 - gamma_ion / (density*alphaB_10k);
        chi = ( -b - sqrt(b*b - 4) ) / 2.0; //correct root
        return chi;
    }
    
    
    
    /* ionization rate (1/s) from neutral fraction (assuming mean density at z) */
    double nftogamma(double nf, float z){
        double nfguess, gammamax, gammamin, gammaguess;
        float epsilon = 1e-7;
        if ( (nf<0) || (nf>1) ){
            fprintf(stderr, "Error in function nftogamma.  Neutral fraction must be between 0 and 1.\n");
            return 0;
        }
        
        if (nf < TINY){
            return 1e20;
        }
        
        gammamax = 1e-4;
        gammamin = 0;
        while (1){
            gammaguess = (gammamax+gammamin)/2.0;
            nfguess = neutral_fraction(No*pow(1+z, 3), gammaguess);
            if ( fabs((nfguess-nf)/nf) > epsilon){  /* guess isn't close enough */
                if (nfguess > nf){ gammamin = gammaguess;}
                else { gammamax = gammaguess;}
            }
            else { break;}  /* guess is acceptable */
            
        }
        
        return gammaguess;
    }
    
    
    /* returns the case A hydrogen recombination coefficient (Abel et al. 1997) in cm^3 s^-1*/
    double alpha_A(double T){
        double logT, ans;
        logT = log(T/(double)1.1604505e4);
        ans = pow(EULER, -28.6130338 - 0.72411256*logT - 2.02604473e-2*pow(logT, 2)
                  - 2.38086188e-3*pow(logT, 3) - 3.21260521e-4*pow(logT, 4)
                  - 1.42150291e-5*pow(logT, 5) + 4.98910892e-6*pow(logT, 6)
                  + 5.75561414e-7*pow(logT, 7) - 1.85676704e-8*pow(logT, 8)
                  - 3.07113524e-9 * pow(logT, 9));
        return ans;
    }
    
    /*
     T in K, M in Msun, mu is mean molecular weight
     from Barkana & Loeb 2001
     
     SUPRESS = 0 for no radiation field supression;
     SUPRESS = 1 for supression (step function at z=z_ss, at v=v_zz)
     */
    float TtoM(float z, float T, float mu){
        return 7030.97 / hlittle * sqrt( omega_mz(z) / (OMm*Deltac_nonlinear(z)) ) *
        pow( T/(mu * (1+z)), 1.5 );
        /*  if (!SUPRESS || (z >= z_re) ) // pre-reionization or don't worry about supression
         return 7030.97 / hlittle * sqrt( omega_mz(z) / (OMm*Deltac_nonlinear(z)) ) *
         pow( T/(mu * (1+z)), 1.5 );
         
         if (z >= z_ss) // self-shielding dominates, use T = 1e4 K
         return 7030.97 / hlittle * sqrt( omega_mz(z) / (OMm*Deltac_nonlinear(z)) ) *
         pow( 1.0e4 /(mu * (1+z)), 1.5 );
         
         // optically thin
         return 7030.97 / hlittle * sqrt( omega_mz(z) / (OMm*Deltac_nonlinear(z)) ) *
         pow( VcirtoT(v_ss, mu) /(mu * (1+z)), 1.5 );
         */
    }
    
    
    /* v in km/s, mu is mean molecular weight */
    float VcirtoT(float v, float mu){
        v *= 1.0e5; // converts to cm/s
        return 0.5 * mu * m_p * v * v / k_B;
    }
    
    /* returns V in km/s
     from Barkana & Loeb 2000 */
    float MtoVcir(float z, double M){
        return 23.4 * pow(M*hlittle/1.0e8, 1.0/3.0) * pow(OMm*Deltac_nonlinear(z)/(omega_mz(z)*18*M_PI*M_PI), 1.0/6.0) * sqrt((1+z)/10.0);
    }
    
    
    /* returns Rvir in comoving Mpc
     from Barkana & Loeb 2000 */
    float MtoRvir(float z, double M){
        return 0.00784/hlittle * pow(M*hlittle/1.0e8, 1.0/3.0) * pow(OMm*Deltac_nonlinear(z)/(omega_mz(z)*18*M_PI*M_PI), -1.0/3.0);
    }
    
    
    
    /* critical density at redshift z Msun Mpc^-3 */
    double rho_critz(float z){
        return RHOcrit * (OMm*pow(1+z,3) + OMl + OMr*pow(1+z,4) + OMk*pow(1+z, 2));
    }
    
    
    /* Physical (non-linear) overdensity at virialization (relative to critical density)
     i.e. answer is rho / rho_crit
     In Einstein de sitter model = 178
     (fitting formula from Bryan & Norman 1998) */
    double Deltac_nonlinear(float z){
        double d;
        d = omega_mz(z) - 1.0;
        return 18*M_PI*M_PI + 82*d - 39*d*d;
    }
    
    
    /* Omega matter at redshift z */
    double omega_mz(float z){
        return OMm*pow(1+z,3) / (OMm*pow(1+z,3) + OMl + OMr*pow(1+z,4) + OMk*pow(1+z, 2));
    }
    
    
    /* Comoving volume (in Mpc^3) probed by past light cone, per z interval */
    double dVdz(float z){
        double d;
        d = d_L(z);
        return 4*M_PI*d*d*C_LIGHT*fabs(dtdz(z)) / (1+z) / CMperMPC;
    }
    
    
    /* arcmin to comoving Mpc */
    double arcmintoMpc(float z, float arcmin){
        return (1+z)*d_A(z)*arcmin*M_PI/180.0/60.0;
    }
    
    /* angular diameter distance (proper Mpc) */
    double d_A(float z){
        return d_L(z)/pow(1.0+z,2);
    }
    
    
    /*  in proper Mpc */
    double d_L_derivs(double z, void * params){
        return (1+z) * fabs(dtdz(z));
    }
    double d_L(float z){
        double result, error;
        gsl_function F;
        double rel_tol  = FRACT_FLOAT_ERR; //<- relative tolerance
        gsl_integration_workspace * w
        = gsl_integration_workspace_alloc (1000);
        
        F.function = &d_L_derivs;
        
        gsl_integration_qag (&F, 0, z, 0, rel_tol,
                             1000, GSL_INTEG_GAUSS61, w, &result, &error);
        gsl_integration_workspace_free (w);
        
        return C_LIGHT * (1+z) * result / CMperMPC;
    }
    
    
    /* function INVSINH returns the inverse hyperbolic sine of parameter x */
    double invsinh (double x){
        return log( x + sqrt(pow(x,2) + 1) );
    }
    
    
    /* function GETTIME returns the age of the universe, at a given redshift parameter, z.
     This assumes zero curvature.  (from Weinberg 1989) */
    double gettime(double z){
        double term1, const1;
        term1 = invsinh( sqrt( OMl/OMm ) * pow(1+z, -3.0/2.0) );
        const1 = 2 * sqrt( 1 + OMm/OMl ) / (3 * Ho) ;
        return (term1 * const1);
    }
    
    
    /* function TIMESINCE returns the time elapsed from parameter redshift zsource to parameter redshift z */
    double timesince(double z, double zsource){
        if (z > zsource){printf("Invalid usage of TIMESINCE. Z must be <= Zsource."); return 0;}
        else if ((zsource-z) < TINY){ return 0;}
        return (gettime(z) - gettime(zsource));
    }
    
    
    /* function TTOZ returns the Z value corresponding to the time, t, elapsed since t(Zsource) */
    double ttoz (double t, double zsource){
        double tguess, zmax, zmin, zguess;
        float epsilon = 1e-6;
        
        if (t < TINY){
            fprintf(stderr, "ERROR: in ttoz: time, t=%.20e, must be greater than zero!\n", t);
            return zsource;
        }
        
        tguess = 0.0;
        zmax = zsource;
        zmin = 0.0;
        
        while (1){
            zguess = (zmax+zmin)/2.0;
            tguess = timesince(zguess, zsource);
            if ( fabs((tguess-t)/t) > epsilon){  /* guess isn't close enough */
                if (tguess > t){ zmin = zguess;}
                else { zmax = zguess;}
            }
            else { break;}  /* guess is acceptable */
        }
        
        return zguess;
    }
    
    
    
    /* function rtoz returns the Z value corresponding to the comoving distance, R, in MPC, traveled by a photon emitted at zsource */
    double rtoz (double R, double zsource){
        double Rguess, zmax, zmin, zguess;
        float epsilon = 1e-6;
        
        if (R < TINY){
            fprintf(stderr, "ERROR: in rtoz: distance, R=%.20e, must be greater than zero!\n", R);
            return zsource;
        }
        
        Rguess = 0.0;
        zmax = zsource;
        zmin = 0.0;
        
        while (1){
            zguess = (zmax+zmin)/2.0;
            Rguess = comovingdistance(zguess, zsource)/CMperMPC;
            if ( fabs((Rguess-R)/R) > epsilon){  /* guess isn't close enough */
                if (Rguess > R){ zmin = zguess;}
                else { zmax = zguess;}
            }
            else { break;}  /* guess is acceptable */
        }
        
        return zguess;
    }
    
    
    /* function DTDZ returns the value of dt/dz at the redshift parameter z. */
    double dtdz(float z){
        double x, dxdz, const1, denom, numer;
        x = sqrt( OMl/OMm ) * pow(1+z, -3.0/2.0);
        dxdz = sqrt( OMl/OMm ) * pow(1+z, -5.0/2.0) * (-3.0/2.0);
        const1 = 2 * sqrt( 1 + OMm/OMl ) / (3.0 * Ho) ;
        
        numer = dxdz * (1 + x*pow( pow(x,2) + 1, -0.5));
        denom = x + sqrt(pow(x,2) + 1);
        return (const1 * numer / denom);
    }
    
    
    
    /* comoving distance (in cm) per unit redshift */
    double drdz(float z){
        return (1.0+z)*C_LIGHT*dtdz(z);
    }
    
    /* Returns the proper distance (in cm) traved by photon (i.e. length of a null geodesic) from zsource to z */
    double properdistance(double z, double zsource){
        return C_LIGHT*timesince(z, zsource);
    }
    
    
    /* Return the comoving distance (in cm) traveled by photon (i.e. length of a null geodesic) from zsource to z */
    double dcom_dz(double z, void *params){
        return drdz(z);
    }
    double comovingdistance(double z, double zsource){
        double result, error;
        gsl_function F;
        double rel_tol  = FRACT_FLOAT_ERR; //<- relative tolerance
        gsl_integration_workspace * w
        = gsl_integration_workspace_alloc (1000);
        
        if (z >= zsource){
            fprintf(stderr, "Error in function comovingdistance! First parameter must be smaller than the second.\n");
            return 0;
        }
        
        F.function = &dcom_dz;
        
        gsl_integration_qag (&F, zsource, z, 0, rel_tol,
                             1000, GSL_INTEG_GAUSS61, w, &result, &error);
        gsl_integration_workspace_free (w);
        
        return result;
    }
    
    
    /* Return the thomspon scattering optical depth from zstart to zend through fully ionized IGM.
     The hydrogen reionization history is given by the zarry and xHarry parameters, in increasing
     redshift order of length len.*/
    typedef struct{
        float *z, *xH;
        int len;
    } tau_e_params;
    double dtau_e_dz(double z, void *params){
        float xH, xi;
        int i=1;
        tau_e_params p = *(tau_e_params *)params;
        
        if ((p.len == 0) || !(p.z))
            return (1+z)*(1+z)*drdz(z);
        else{
            // find where we are in the redshift array
            if (p.z[0]>z) // ionization fraction is 1 prior to start of array
                return (1+z)*(1+z)*drdz(z);
            while ( (i < p.len) && (p.z[i] < z) ) {i++;}
            if (i == p.len)
                return 0;
            
            // linearly interpolate in redshift
            xH = p.xH[i-1] + (p.xH[i] - p.xH[i-1])/(p.z[i] - p.z[i-1]) * (z - p.z[i-1]);
            /*
             fprintf(stderr, "in taue: Interpolating between xH(%f)=%f and xH(%f)=%f to obtain xH(%f)=%f\n",
             p.z[i-1], p.xH[i-1], p.z[i], p.xH[i], z, xH);
             */
            xi = 1.0-xH;
            if (xi<0){
                fprintf(stderr, "in taue: funny buisness xi=%e, changing to 0\n", xi);
                xi=0;
            }
            if (xi>1){
                fprintf(stderr, "in taue: funny buisness xi=%e, changing to 1\n", xi);
                xi=1;
            }
            
            return xi*(1+z)*(1+z)*drdz(z);
        }
    }
    double tau_e(float zstart, float zend, float *zarry, float *xHarry, int len){
        double prehelium, posthelium, error;
        gsl_function F;
        double rel_tol  = 1e-3; //<- relative tolerance
        gsl_integration_workspace * w
        = gsl_integration_workspace_alloc (1000);
        tau_e_params p;
        
        if (zstart >= zend){
            fprintf(stderr, "Error in function taue! First parameter must be smaller than the second.\n");
            return 0;
        }
        
        F.function = &dtau_e_dz;
        p.z = zarry;
        p.xH = xHarry;
        p.len = len;
        F.params = &p;
        if ((len > 0) && zarry)
            zend = zarry[len-1] - FRACT_FLOAT_ERR;
        
        
        
        if (zend > Zreion_HeII){// && (zstart < Zreion_HeII)){
            if (zstart < Zreion_HeII){
                gsl_integration_qag (&F, Zreion_HeII, zstart, 0, rel_tol,
                                     1000, GSL_INTEG_GAUSS61, w, &prehelium, &error);
                gsl_integration_qag (&F, zend, Zreion_HeII, 0, rel_tol,
                                     1000, GSL_INTEG_GAUSS61, w, &posthelium, &error);
            }
            else{
                prehelium = 0;
                gsl_integration_qag (&F, zend, zstart, 0, rel_tol,
                                     1000, GSL_INTEG_GAUSS61, w, &posthelium, &error);
            }
        }
        else{
            posthelium = 0;
            gsl_integration_qag (&F, zend, zstart, 0, rel_tol,
                                 1000, GSL_INTEG_GAUSS61, w, &prehelium, &error);
        }
        gsl_integration_workspace_free (w);
        
        return SIGMAT * ( (N_b0+He_No)*prehelium + N_b0*posthelium );
    }
    
    
    /* Function Lya_crosssec returns the differential Lya absorption crosssection (in cm) at frequency nu, and gas temperature T (in Kelvin) */
    double Lya_crosssec(double nu, double T){
        double dopwidth, avoigt, x;
        
        dopwidth = Lya_dopwidth(T);
        avoigt = 49854330.99 / dopwidth;
        x = ( nu - Ly_alpha_HZ ) / dopwidth;
        return (0.011046) * voigt1( x, avoigt ) / dopwidth;
    }
    
    
    double Lya_dopwidth(double T){
        return 1056924075.0 * sqrt(T);
    }
    
    /* function VOIGT1 computes the normalized voigt function for any value of x and any
     positive value of a.
     this function is the C port of the fortran function written by G. Rybicki */
    double voigt1(double x, double a){
        
        //double avoigt,
        double cvoigt[31];
        int VOIGT_INIT = -1;
        
        double a1, a2, b1, b2, e, q1, q2, s, t,zi, zr, vg1;
        int n, k;
        
        if (VOIGT_INIT != 10){ // first time in function, initialize cvoigt
            k = -16;
            for (n=0; n<31; n++){
                k++;
                cvoigt[n] = 0.0897935610625833 * exp(-k*k/9.0);
            }
            VOIGT_INIT = 10;
        }
        
        
        q1 = 9.42477796076938;
        q2 = 0.564189583547756;
        
        /* check if a = 0 (Doppler) */
        if ( fabs(a) < TINY ){
            return ( q2 * exp(-x*x) );
        }
        
        /* calculation for the general case */
        a1 = 3*a;
        a2 = a*a;
        e = exp(-q1*a);
        if (a < 0.1){
            zr = 0.5 * (e + 1.0/e) * cos(q1*x);
            zi = 0.5 * (e - 1.0/e) * sin(q1*x);
            vg1 = q2 * exp(a2 - x*x) * cos(2*a*x);
        }
        else{
            zr = e * cos(q1*x);
            zi = e * sin(q1*x);
            vg1 = 0;
        }
        
        b1 = (1 - zr) * a * 1.5;
        b2 = -zi;
        s = -8 - 1.5 * x;
        t = s*s + 2.25 * a2;
        for (n=0; n<31; n++){
            t += s + 0.25;
            s += 0.5;
            b1 = a1 - b1;
            b2 = -b2;
            if (t > 2.5e-12)
                vg1 += cvoigt[n] * (b1 + b2*s) / t;
            else
                vg1 -= cvoigt[n] * a * 29.608813203268;
        }
        
        return vg1;
    }
    
    
    /* function HeI_ion_crosssec returns the HI ionization cross section at parameter frequency
     (taken from Verner et al (1996) */
    double HeI_ion_crosssec(double nu){
        double x,y;//,Fy;
        
        if (nu < HeI_NUIONIZATION)
            return 0;
        
        x = nu/NU_over_EV/13.61 - 0.4434;
        y = sqrt(x*x + pow(2.136, 2));
        return  9.492e-16*((x-1)*(x-1) + 2.039*2.039) *
        pow(y, (0.5 * 3.188 - 5.5))
        * pow(1.0 + sqrt(y/1.469), -3.188);
    }
    
    
    /* function HeII_ion_crosssec returns the HeII ionization cross section at parameter frequency
     (taken from Osterbrock, pg. 14) */
    double HeII_ion_crosssec(double nu){
        double epsilon, Z = 2;
        
        if (nu < HeII_NUIONIZATION)
            return 0;
        
        if (nu == HeII_NUIONIZATION)
            nu+=TINY;
        
        epsilon = sqrt(nu/HeII_NUIONIZATION - 1);
        return (6.3e-18)/Z/Z * pow(HeII_NUIONIZATION/nu, 4)
        * pow(EULER, 4-(4*atan(epsilon)/epsilon)) / (1-pow(EULER, -2.*M_PI/epsilon));
    }
    
    
    /* function HI_ion_crosssec returns the HI ionization cross section at parameter frequency
     (taken from Osterbrock, pg. 14) */
    double HI_ion_crosssec(double nu){
        double epsilon, Z = 1;
        
        if (nu < NUIONIZATION)
            return 0;
        
        if (nu == NUIONIZATION)
            nu+=TINY;
        
        epsilon = sqrt( nu/NUIONIZATION - 1);
        return (6.3e-18)/Z/Z * pow(NUIONIZATION/nu, 4)
        * pow(EULER, 4-(4*atan(epsilon)/epsilon)) / (1-pow(EULER, -2.*M_PI/epsilon));
    }
    
} /* namespace */






