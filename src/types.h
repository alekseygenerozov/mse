#include <math.h>
#include <cstdlib>
#include <map>
#include <vector>


extern "C"
{
#define CVODE_EXTENDED_PRECISION
 
//#define VERBOSE
 
#ifndef __FOUND_ROOT
#define ___FOUND_ROOT
#define FOUND_ROOT ((roots_found[i_root] == 1) || (roots_found[i_root] == -1))
#endif 

#ifndef __CONSTANTS
#define __CONSTANTS

// Default constants //
extern double CONST_G;
extern double CONST_G_P2;
extern double CONST_G_P3;
extern double CONST_C_LIGHT;
extern double CONST_C_LIGHT_P2;
extern double CONST_C_LIGHT_P3;
extern double CONST_C_LIGHT_P4;
extern double CONST_C_LIGHT_P5;
extern double CONST_C_LIGHT_P6;
extern double CONST_MSUN;
extern double CONST_R_SUN;
extern double CONST_L_SUN;
extern double CONST_KM_PER_S;
extern double CONST_PER_PC3;

// Default parameters //
extern double relative_tolerance;
extern double absolute_tolerance_eccentricity_vectors;
extern bool include_quadrupole_order_terms;
extern bool include_octupole_order_binary_pair_terms;
extern bool include_octupole_order_binary_triplet_terms;
extern bool include_hexadecupole_order_binary_pair_terms;
extern bool include_dotriacontupole_order_binary_pair_terms;
extern bool include_double_averaging_corrections;
extern double epsilon;
extern int random_seed;

extern bool include_flybys;
extern bool flybys_correct_for_gravitational_focussing;
extern int flybys_velocity_distribution;
extern int flybys_mass_distribution;
extern int flybys_N_enc;
extern int flybys_N_not_impulsive;
extern int flybys_reference_binary;
extern double flybys_t_next_encounter;
extern double flybys_mass_distribution_lower_value;
extern double flybys_mass_distribution_upper_value;
extern double flybys_encounter_sphere_radius;
extern double flybys_stellar_density;
extern double flybys_stellar_relative_velocity_dispersion;
extern double flybys_W_max;
extern double flybys_total_encounter_rate_at_R_enc;
extern double flybys_stellar_density_at_R_enc;
extern double flybys_internal_mass;
extern double flybys_internal_semimajor_axis;

extern int binary_evolution_CE_energy_flag;
extern int binary_evolution_CE_spin_flag;

extern double chandrasekhar_mass;
extern double eddington_accretion_factor;
extern double nova_accretion_factor;
extern double alpha_wind_accretion;
extern double beta_wind_accretion;

extern double mstar_gbs_tolerance;
extern double mstar_collision_tolerance;

extern double triple_mass_transfer_primary_star_accretion_efficiency_no_disk;
extern double triple_mass_transfer_secondary_star_accretion_efficiency_no_disk;
extern double triple_mass_transfer_primary_star_accretion_efficiency_disk;
extern double triple_mass_transfer_secondary_star_accretion_efficiency_disk;
extern double triple_mass_transfer_inner_binary_alpha_times_lambda;

extern double nbody_analysis_fractional_semimajor_axis_change_parameter;
extern double nbody_analysis_fractional_integration_time;
extern double nbody_analysis_maximum_integration_time;

extern double nbody_dynamical_instability_direct_integration_time_multiplier;
extern double nbody_semisecular_direct_integration_time_multiplier;
extern double nbody_supernovae_direct_integration_time_multiplier;
extern double nbody_other_direct_integration_time_multiplier;

#ifdef IGNORE
#define CONST_C_LIGHT		(double)	63239726386.8
#define CONST_C_LIGHT_P2	(double)	CONST_C_LIGHT*CONST_C_LIGHT
#define CONST_C_LIGHT_P4	(double)	CONST_C_LIGHT_P2*CONST_C_LIGHT_P2
#define CONST_C_LIGHT_P5	(double)	CONST_C_LIGHT_P4*CONST_C_LIGHT
#define CONST_MSUN          (double)    1.0
#define CONST_R_SUN         (double)    0.00464913034382
#define CONST_L_SUN         (double)    2.71040410975e+14
#endif

#define TWOPI               (double)    2.0 * M_PI
#define yr_to_Myr           (double)    1.0e-6
#define Myr_to_yr           (double)    1.0e6
#define c_1div2             (double)    1.0/2.0
#define c_1div3             (double)    1.0/3.0
#define c_1div4             (double)    1.0/4.0
#define c_1div5             (double)    1.0/5.0
#define c_1div6             (double)    1.0/6.0
#define c_1div7             (double)    1.0/7.0
#define c_1div8             (double)    1.0/8.0
#define c_1div10            (double)    1.0/10.0
#define c_1div15            (double)    1.0/15.0
#define c_1div16            (double)    1.0/16.0
#define c_1div30            (double)    1.0/30.0
#define c_2div3             (double)    2.0/3.0
#define c_3div2             (double)    3.0/2.0
#define c_3div4             (double)    3.0/4.0
#define c_3div5             (double)    3.0/5.0
#define c_3div8             (double)    3.0/8.0
#define c_3div32            (double)    3.0/32.0
#define c_3div1024          (double)    3.0/1024.0
#define c_4div15            (double)    4.0/15.0
#define c_5div2             (double)    5.0/2.0
#define c_5div8             (double)    5.0/8.0
#define c_5div16            (double)    5.0/16.0
#define c_5div64            (double)    5.0/64.0
#define c_7div8             (double)    7.0/8.0
#define c_8div5             (double)    8.0/5.0
#define c_8div7             (double)    8.0/7.0
#define c_9div2             (double)    9.0/2.0
#define c_9div8             (double)    9.0/8.0
#define c_9div16            (double)    9.0/16.0
#define c_9div32            (double)    9.0/32.0
#define c_9div64            (double)    9.0/64.0
#define c_10div3            (double)    10.0/3.0
#define c_11div18           (double)    11.0/18.0
#define c_15div2            (double)    15.0/2.0
#define c_15div4            (double)    15.0/4.0
#define c_15div8            (double)    15.0/8.0
#define c_15div16           (double)    15.0/16.0
#define c_16div5            (double)    16.0/5.0
#define c_25div16           (double)    25.0/16.0
#define c_25div64           (double)    25.0/64.0
#define c_27div4            (double)    27.0/4.0
#define c_27div64           (double)    27.0/64.0
#define c_31div2            (double)    31.0/2.0
#define c_32div5            (double)    32.0/5.0
#define c_35div3            (double)    35.0/3.0
#define c_37div96           (double)    37.0/96.0
#define c_45div8            (double)    45.0/8.0
#define c_64div3            (double)    64.0/3.0
#define c_64div5            (double)    64.0/5.0
#define c_73div24           (double)    73.0/24.0
#define c_105div4096        (double)    105.0/4096.0
#define c_121div304         (double)    121.0/304.0
#define c_185div16          (double)    185.0/16.0
#define c_255div8           (double)    255.0/8.0
#define c_304div15          (double)    304.0/15.0
#endif

#ifndef __TABLES
#define __TABLES
#define MAX_ORDER (int) 5

#define TABLEWIDTH_A (int) 3
#define TABLELENGTH_A (int) 10
const double A_TABLE[TABLELENGTH_A][TABLEWIDTH_A] =
{{2, 0, -0.500000000000000}, {2, 2, 1.50000000000000}, {3, 
  1, -1.50000000000000}, {3, 3, 2.50000000000000}, {4, 0, 
  0.375000000000000}, {4, 2, -3.75000000000000}, {4, 4, 
  4.37500000000000}, {5, 1, 1.87500000000000}, {5, 
  3, -8.75000000000000}, {5, 5, 7.87500000000000}};

#define A_mn_table_max_order (int) 15
const double A_mn_table[15][15+1] = {{-0.500000000000000, 0, 1.50000000000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0}, {0, -1.50000000000000, 0, 2.50000000000000, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0}, {0.375000000000000, 0, -3.75000000000000, 
  0, 4.37500000000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 
  1.87500000000000, 0, -8.75000000000000, 0, 7.87500000000000, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0}, {-0.312500000000000, 0, 6.56250000000000, 
  0, -19.6875000000000, 0, 14.4375000000000, 0, 0, 0, 0, 0, 0, 0, 0, 
  0}, {0, -2.18750000000000, 0, 19.6875000000000, 
  0, -43.3125000000000, 0, 26.8125000000000, 0, 0, 0, 0, 0, 0, 0, 
  0}, {0.273437500000000, 0, -9.84375000000000, 0, 54.1406250000000, 
  0, -93.8437500000000, 0, 50.2734375000000, 0, 0, 0, 0, 0, 0, 0}, {0,
   2.46093750000000, 0, -36.0937500000000, 0, 140.765625000000, 
  0, -201.093750000000, 0, 94.9609375000000, 0, 0, 0, 0, 0, 
  0}, {-0.246093750000000, 0, 13.5351562500000, 0, -117.304687500000, 
  0, 351.914062500000, 0, -427.324218750000, 0, 180.425781250000, 0, 
  0, 0, 0, 0}, {0, -2.70703125000000, 0, 58.6523437500000, 
  0, -351.914062500000, 0, 854.648437500000, 0, -902.128906250000, 0, 
  344.449218750000, 0, 0, 0, 0}, {0.225585937500000, 
  0, -17.5957031250000, 0, 219.946289062500, 0, -997.089843750000, 0, 
  2029.79003906250, 0, -1894.47070312500, 0, 660.194335937500, 0, 0, 
  0}, {0, 2.93261718750000, 0, -87.9785156250000, 0, 747.817382812500,
   0, -2706.38671875000, 0, 4736.17675781250, 0, -3961.16601562500, 0,
   1269.60449218750, 0, 0}, {-0.209472656250000, 0, 21.9946289062500, 
  0, -373.908691406250, 0, 2368.08837890625, 0, -7104.26513671875, 0, 
  10893.2065429688, 0, -8252.42919921875, 0, 2448.52294921875, 
  0}, {0, -3.14208984375000, 0, 124.636230468750, 
  0, -1420.85302734375, 0, 7104.26513671875, 0, -18155.3442382813, 0, 
  24757.2875976563, 0, -17139.6606445313, 0, 4733.81103515625}};
  
#define TABLEWIDTH_B (int) 8
#define TABLELENGTH_B (int) 47
#define HIGHEST_POWER_ECCP2_IN_B_TABLE (int) 3
const double B_TABLE[TABLELENGTH_B][TABLEWIDTH_B] =
{{2, 0, 0, 0, 1.00000000000000, 1.50000000000000, 0, 0}, {2, 1, 1, 
  0, -2.00000000000000, -0.500000000000000, 0, 0}, {2, 2, 0, 0, 
  0.500000000000000, -0.500000000000000, 0, 0}, {2, 2, 0, 
  2, -0.500000000000000, 0, 0, 0}, {2, 2, 2, 0, 2.50000000000000, 0, 
  0, 0}, {3, 0, 0, 0, 1.00000000000000, 3.00000000000000, 
  0.375000000000000, 0}, {3, 1, 1, 
  0, -2.50000000000000, -1.87500000000000, 0, 0}, {3, 2, 0, 0, 
  0.500000000000000, -0.375000000000000, -0.125000000000000, 0}, {3, 
  2, 0, 2, -0.500000000000000, -0.125000000000000, 0, 0}, {3, 2, 2, 0,
   3.75000000000000, 0.625000000000000, 0, 0}, {3, 3, 1, 
  0, -1.87500000000000, 1.87500000000000, 0, 0}, {3, 3, 1, 2, 
  1.87500000000000, 0, 0, 0}, {3, 3, 3, 0, -4.37500000000000, 0, 0, 
  0}, {4, 0, 0, 0, 1.00000000000000, 5.00000000000000, 
  1.87500000000000, 0}, {4, 1, 1, 
  0, -3.00000000000000, -4.50000000000000, -0.375000000000000, 0}, {4,
   2, 0, 0, 0.500000000000000, -0.125000000000000, -0.375000000000000,
   0}, {4, 2, 0, 2, -0.500000000000000, -0.375000000000000, 0, 0}, {4,
   2, 2, 0, 5.25000000000000, 2.62500000000000, 0, 0}, {4, 3, 1, 
  0, -2.25000000000000, 1.87500000000000, 0.375000000000000, 0}, {4, 
  3, 1, 2, 2.25000000000000, 0.375000000000000, 0, 0}, {4, 3, 3, 
  0, -7.00000000000000, -0.875000000000000, 0, 0}, {4, 4, 0, 0, 
  0.375000000000000, -0.750000000000000, 0.375000000000000, 0}, {4, 4,
   0, 2, -0.750000000000000, 0.750000000000000, 0, 0}, {4, 4, 0, 4, 
  0.375000000000000, 0, 0, 0}, {4, 4, 2, 0, 
  5.25000000000000, -5.25000000000000, 0, 0}, {4, 4, 2, 
  2, -5.25000000000000, 0, 0, 0}, {4, 4, 4, 0, 7.87500000000000, 0, 0,
   0}, {5, 0, 0, 0, 1.00000000000000, 7.50000000000000, 
  5.62500000000000, 0.312500000000000}, {5, 1, 1, 
  0, -3.50000000000000, -8.75000000000000, -2.18750000000000, 0}, {5, 
  2, 0, 0, 0.500000000000000, 
  0.250000000000000, -0.687500000000000, -0.0625000000000000}, {5, 2, 
  0, 2, -0.500000000000000, -0.750000000000000, -0.0625000000000000, 
  0}, {5, 2, 2, 0, 7.00000000000000, 7.00000000000000, 
  0.437500000000000, 0}, {5, 3, 1, 0, -2.62500000000000, 
  1.31250000000000, 1.31250000000000, 0}, {5, 3, 1, 2, 
  2.62500000000000, 1.31250000000000, 0, 0}, {5, 3, 3, 
  0, -10.5000000000000, -3.93750000000000, 0, 0}, {5, 4, 0, 0, 
  0.375000000000000, -0.687500000000000, 0.250000000000000, 
  0.0625000000000000}, {5, 4, 0, 2, -0.750000000000000, 
  0.625000000000000, 0.125000000000000, 0}, {5, 4, 0, 4, 
  0.375000000000000, 0.0625000000000000, 0, 0}, {5, 4, 2, 0, 
  7.00000000000000, -6.12500000000000, -0.875000000000000, 0}, {5, 4, 
  2, 2, -7.00000000000000, -0.875000000000000, 0, 0}, {5, 4, 4, 0, 
  13.1250000000000, 1.31250000000000, 0, 0}, {5, 5, 1, 
  0, -2.18750000000000, 4.37500000000000, -2.18750000000000, 0}, {5, 
  5, 1, 2, 4.37500000000000, -4.37500000000000, 0, 0}, {5, 5, 1, 
  4, -2.18750000000000, 0, 0, 0}, {5, 5, 3, 0, -13.1250000000000, 
  13.1250000000000, 0, 0}, {5, 5, 3, 2, 13.1250000000000, 0, 0, 
  0}, {5, 5, 5, 0, -14.4375000000000, 0, 0, 0}};

#define TABLEWIDTH_D (int) 8
#define TABLELENGTH_D (int) 15
const double D_TABLE[TABLELENGTH_D][TABLEWIDTH_D] = 
{
    {2, 0, 0, 0, 0, 0, 0, 1}, \
    {2, 0, 2, 0, 0, 0, 0, 2}, \
    {2, 0, 2, 0, 2, 0, 0, 3}, \
    {2, 0, 2, 0, 0, 0, 2, 4}, \
    {2, 2, 0, 0, 0, 0, 0, 5}, \
    {2, 2, 0, 2, 0, 0, 0, 6}, \
    {2, 2, 0, 0, 0, 2, 0, 7}, \
    {3, 1, 0, 1, 0, 0, 0, 8}, \
    {3, 1, 2, 0, 1, 1, 1, 9}, \
    {3, 1, 2, 1, 0, 0, 0, 10}, \
    {3, 1, 2, 1, 0, 0, 2, 11}, \
    {3, 1, 2, 1, 2, 0, 0, 12}, \
    {3, 3, 0, 1, 0, 0, 0, 13}, \
    {3, 3, 0, 1, 0, 2, 0, 14}, \
    {3, 3, 0, 3, 0, 0, 0, 15}
};

/* the numbers in each last entry of the D table refer to the functions defined below */

inline double D_TABLE_FUNC1(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 2.0*(sqrt_ef_p2_minus_one + asec_minus_ef);
}
inline double D_TABLE_FUNC1_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

inline double D_TABLE_FUNC2(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return (1.0 - ep_p2)*( c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + asec_minus_ef);
}
inline double D_TABLE_FUNC2_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return -2.0*ep*( c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + asec_minus_ef);
}

inline double D_TABLE_FUNC3(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return -c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) - asec_minus_ef;
}
inline double D_TABLE_FUNC3_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

inline double D_TABLE_FUNC4(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_2div3*one_div_ef_p2*ef_p2_minus_one*sqrt_ef_p2_minus_one;
}
inline double D_TABLE_FUNC4_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

inline double D_TABLE_FUNC5(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return ep_p2*( c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + asec_minus_ef);
}
inline double D_TABLE_FUNC5_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 2.0*ep*( c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + asec_minus_ef);
}

inline double D_TABLE_FUNC6(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_2div3*one_div_ef_p2*ef_p2_minus_one*sqrt_ef_p2_minus_one;
}
inline double D_TABLE_FUNC6_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

inline double D_TABLE_FUNC7(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return -c_1div3*one_div_ef_p2*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) - asec_minus_ef;
}
inline double D_TABLE_FUNC7_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC8(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 2.0*(c_1div3*one_div_ef_p1*sqrt_ef_p2_minus_one*(1.0 + 2.0*ef_p2) + ef*asec_minus_ef);
}
inline double D_TABLE_FUNC8_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC9(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_1div15*one_div_ef_p3*sqrt_ef_p2_minus_one*(2.0 - 9.0*ef_p2 - 8.0*ef_p4) - ef*asec_minus_ef;
}
inline double D_TABLE_FUNC9_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC10(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return (1.0 - ep_p2)*( c_1div30*one_div_ef_p3*sqrt_ef_p2_minus_one*(-2.0 + 9.0*ef_p2 + 8.0*ef_p4) + c_1div2*ef*asec_minus_ef);
}
inline double D_TABLE_FUNC10_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return -2.0*ep*( c_1div30*one_div_ef_p3*sqrt_ef_p2_minus_one*(-2.0 + 9.0*ef_p2 + 8.0*ef_p4) + c_1div2*ef*asec_minus_ef);
}
inline double D_TABLE_FUNC11(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_4div15*one_div_ef_p3*ef_p2_minus_one*ef_p2_minus_one*sqrt_ef_p2_minus_one;
}
inline double D_TABLE_FUNC11_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC12(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_1div30*one_div_ef_p3*sqrt_ef_p2_minus_one*(2.0 - 9.0*ef_p2 - 8.0*ef_p4) - c_1div2*ef*asec_minus_ef;
}
inline double D_TABLE_FUNC12_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC13(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return ep_p2*( c_1div10*one_div_ef_p3*sqrt_ef_p2_minus_one*(-2.0 + 9.0*ef_p2 + 8.0*ef_p4) + c_3div2*ef*asec_minus_ef );
}
inline double D_TABLE_FUNC13_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 2.0*ep*( c_1div10*one_div_ef_p3*sqrt_ef_p2_minus_one*(-2.0 + 9.0*ef_p2 + 8.0*ef_p4) + c_3div2*ef*asec_minus_ef );
}
inline double D_TABLE_FUNC14(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_1div10*one_div_ef_p3*sqrt_ef_p2_minus_one*(2.0 - 9.0*ef_p2 - 8.0*ef_p4) - c_3div2*ef*asec_minus_ef;
}
inline double D_TABLE_FUNC14_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}
inline double D_TABLE_FUNC15(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return c_4div15*one_div_ef_p3*ef_p2_minus_one*ef_p2_minus_one*sqrt_ef_p2_minus_one;
}
inline double D_TABLE_FUNC15_DER(double ep, double ep_p2, double ef, double ef_p2, double ef_p4, double one_div_ef_p1, double one_div_ef_p2, double one_div_ef_p3, double ef_p2_minus_one, double sqrt_ef_p2_minus_one, double asec_minus_ef)
{
    return 0.0;
}

#endif

/*	ODE solver macros	*/
#ifndef __ODE_MACROS
#define __ODE_MACROS
    #define Ith(v,i)    NV_Ith_S(v,i-1)       		/* Ith numbers components 1..NEQ */
    #define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) 		/* IJth numbers rows,cols 1..NEQ */
    #ifndef CV_max
        #define CV_max( a, b ) ( ((a) > (b)) ? (a) : (b) )
    #endif
    #ifndef CV_min
        #define CV_min(X,Y) ((X) < (Y) ? (X) : (Y))
    #endif
#endif

/* vector operators */
#ifndef __VECTOR_OPERATORS
#define __VECTOR_OPERATORS
inline void cross3(double a[3], double b[3], double result[3])
{
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}
inline double norm3(double v[3])
{
    double result = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    return result;
}
inline double norm3_squared(double v[3])
{
    double result = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    return result;
}
inline double dot3(double a[3], double b[3])
{
    double result = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
    return result;
}
inline void matrix_on_vector( int N_x, int N_y, double M[][3], double v[], double result[] )
{
    
    int i,j;

    for (j=0; j<N_y; j++)
    {
        result[j] = 0.0;
    }

    for (j=0; j<N_y; j++)
    {
        for (i=0; i<N_x; i++)
        {
            result[j] += M[i][j]*v[i];
        }
    }
}

inline void scalar_times_matrix(int N_x, int N_y, double M[][3], double common_factor)
{
    int i,j;
    for (i=0; i<N_x; i++)
    {
        for (j=0; j<N_y; j++)
        {
            M[i][j] *= common_factor;
        }
    }
}

/* KS regularization */
inline double dot4(double a[4], double b[4])
{
    double result = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]);
    return result;
}
inline double norm4(double v[4])
{
    double result = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
    return result;
}
inline double norm4_squared(double v[4])
{
    double result = v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
    return result;
}   

inline void transform_r_to_u_d0(double r[3], double u_d0[4])
{
  if (r[0] >= 0.0)
  {
    u_d0[0] = sqrt(0.5*(r[0] + norm3(r)));
    u_d0[1] = 0.5*r[1]/u_d0[0];
    u_d0[2] = 0.5*r[2]/u_d0[0];
    u_d0[3] = 0.0;
  }
  else if (r[0] < 0.0)
  {
    u_d0[1] = sqrt(0.5*(norm3(r) - r[0]));
    u_d0[0] = 0.5*r[1]/u_d0[1];
    u_d0[3] = 0.5*r[2]/u_d0[1];
    u_d0[2] = 0.0;
  }
}

inline void transform_v_to_u_d1(double u_d0[4], double v[3], double u_d1[4])
{
  u_d1[0] = 0.5*(u_d0[0]*v[0] + u_d0[1]*v[1] + u_d0[2]*v[2]);
  u_d1[1] = 0.5*(-u_d0[1]*v[0] + u_d0[0]*v[1] + u_d0[3]*v[2]);
  u_d1[2] = 0.5*(-u_d0[2]*v[0] - u_d0[3]*v[1] + u_d0[0]*v[2]);
  u_d1[3] = 0.5*(u_d0[3]*v[0] - u_d0[2]*v[1] + u_d0[1]*v[2]);
}

inline void transform_u_d0_to_r(double u_d0[4], double r[3])
{
  r[0] = u_d0[0]*u_d0[0] - u_d0[1]*u_d0[1] - u_d0[2]*u_d0[2] + u_d0[3]*u_d0[3];
  r[1] = 2.0*(u_d0[0]*u_d0[1] - u_d0[2]*u_d0[3]);
  r[2] = 2.0*(u_d0[0]*u_d0[2] + u_d0[1]*u_d0[3]);
}

inline void transform_u_d1_to_v(double u_d0[4], double u_d1[4], double v[3])
{
	double r_norm = norm4_squared(u_d0);
  v[0] = 2.0*(u_d0[0]*u_d1[0] - u_d0[1]*u_d1[1] - u_d0[2]*u_d1[2] + u_d0[3]*u_d1[3])/r_norm;
  v[1] = 2.0*(u_d0[1]*u_d1[0] + u_d0[0]*u_d1[1] - u_d0[3]*u_d1[2] - u_d0[2]*u_d1[3])/r_norm;
  v[2] = 2.0*(u_d0[2]*u_d1[0] + u_d0[3]*u_d1[1] + u_d0[0]*u_d1[2] + u_d0[1]*u_d1[3])/r_norm;
}

inline void transform_u_star_to_v(double u[4], double u_star[4], double omega, double v[3])
{
    double r_norm = norm4_squared(u);
    double factor = 4.0*omega/r_norm;
    
    v[0] = factor*(u[0]*u_star[0] - u[1]*u_star[1] - u[2]*u_star[2] + u[3]*u_star[3]);
    v[1] = factor*(u[1]*u_star[0] + u[0]*u_star[1] - u[3]*u_star[2] - u[2]*u_star[3]);
    v[2] = factor*(u[2]*u_star[0] + u[3]*u_star[1] + u[0]*u_star[2] + u[1]*u_star[3]);
}
    
inline void transform_v_to_u_star(double u[4], double v[3], double omega, double u_star[4])
{
    double factor = 0.25/omega;
    u_star[0] = factor*(u[0]*v[0] + u[1]*v[1] + u[2]*v[2]);
    u_star[1] = factor*(-u[1]*v[0] + u[0]*v[1] + u[3]*v[2]);
    u_star[2] = factor*(-u[2]*v[0] - u[3]*v[1] + u[0]*v[2]);
    u_star[3] = factor*(u[3]*v[0] - u[2]*v[1] + u[1]*v[2]);
}

inline void LT_u_on_vec3(double u_d0[4], double vec3[3], double result[4])
{
  result[0] = u_d0[0]*vec3[0] + u_d0[1]*vec3[1] + u_d0[2]*vec3[2];
  result[1] = -u_d0[1]*vec3[0] + u_d0[0]*vec3[1] + u_d0[3]*vec3[2];
  result[2] = -u_d0[2]*vec3[0] - u_d0[3]*vec3[1] + u_d0[0]*vec3[2];
  result[3] = u_d0[3]*vec3[0] - u_d0[2]*vec3[1] + u_d0[1]*vec3[2];
}

inline void transform_alpha_beta_to_u_u_star(double alpha[4], double beta[4], double E, double u[4], double u_star[4])
{
    double cos_E_div2 = cos(c_1div2*E);
    double sin_E_div2 = sin(c_1div2*E);
    for (int i=0; i<4; i++)
    {
        u[i] = alpha[i]*cos_E_div2 + beta[i]*sin_E_div2;
        u_star[i] = -c_1div2*alpha[i]*sin_E_div2 + c_1div2*beta[i]*cos_E_div2;
    }
}
#endif


/* sse interface */
#ifndef __STEV
#define __STEV

#define SSE_LOW_MASS_MS 0
#define SSE_MS 1
#define SSE_HG 2
#define SSE_FGB 3
#define SSE_CHeB 4
#define SSE_EAGB 5
#define SSE_TPAGB 6
#define SSE_HeMS 7
#define SSE_HeHG 8
#define SSE_HeGB 9
#define SSE_HeWD 10
#define SSE_COWD 11
#define SSE_ONeWD 12
#define SSE_NS 13
#define SSE_BH 14
#define SSE_MASSLESS_REMNANT 15

struct value1__ 
{
        double neta,bwind,hewind,mxns;
};

struct value2__ 
{
        double alpha1,lambda;
};

struct value3__ 
{
        int idum;
};

struct value4__ 
{
        double sigma;
        int bhflag;
};

struct value5__ 
{
        double beta,xi,acc2,epsnov,eddfac,gamma;
};

struct flags__ 
{
        int ceflag,tflag,ifflag,nsflag,wdflag;
};

struct points__ 
{
        double pts1,pts2,pts3;
};

extern "C" void evolv1_(int *kw, double *mass,double *mt, double *r, double *lum, double *mc, double *rc, double *menv, double *renv, double *ospin, double *epoch, double *tms, double *tphys, double *tphysf, double *dtp, double *z, double *zpars, double *k2);
extern "C" void zcnsts_(double* z, double *zpars);
extern "C" void star_(int *kw, double *mass, double *mt, double *tm, double *tn, double *tscls, double *lums, double *GB, double *zpars);
extern "C" void deltat_(int *kw, double *age, double *tm, double *tn, double *tscls, double *dt, double *dtr);
extern "C" void hrdiag_(double *mass, double *aj, double *mt, double *tm, double *tn, double * tscls, double *lums, double *GB, double *zpars,
    double *r, double *lum, int *kw, double *mc, double *rc, double *menv, double *renv, double *k2);
extern "C" double vrotf_(double mt);    
extern "C" double rzamsf_(double *m);
extern "C" double celamf_(int *kw, double *m, double *lum, double *rad, double *rzams, double *menv, double *fac);
//extern "C" int ktype_(int *kw1, int *kw2);
extern "C" void dgcore_(int *kw1, int *kw2, int *kw3, double *m1, double *m2, double *m3, double *ebinde);
extern "C" void gntage_(double *mc, double *mt, int *kw, double *zpars, double *m0, double *aj);

const int MERGER_TABLE[15][15] =
{
   {1,1,2,3,4,5,6,4,6,6,3,6,6,13,14}, \
   {1,1,2,3,4,5,6,4,6,6,3,6,6,13,14}, \
   {2,2,3,3,4,4,5,4,4,4,3,5,5,13,14}, \
   {3,3,3,3,4,4,5,4,4,4,3,5,5,13,14}, \
   {4,4,4,4,4,4,4,4,4,4,4,4,4,13,14}, \
   {5,5,4,4,4,4,4,4,4,4,4,4,4,13,14}, \
   {6,6,5,5,4,4,6,4,6,6,5,6,6,13,14}, \
   {4,4,4,4,4,4,4,7,8,9,7,9,9,13,14}, \
   {6,6,4,4,4,4,6,8,8,9,7,9,9,13,14}, \
   {6,6,4,4,4,4,6,9,9,9,7,9,9,13,14}, \
   {3,3,3,3,4,4,5,7,7,7,15,9,9,13,14}, \
   {6,6,5,5,4,4,6,9,9,9,9,11,12,13,14}, \
   {6,6,5,5,4,4,6,9,9,9,9,12,12,13,14}, \
   {13,13,13,13,13,13,13,13,13,13,13,13,13,13,14}, \
   {14,14,14,14,14,14,14,14,14,14,14,14,14,14,14}
};

#endif

/* Used in apsidal_motion_constant.cpp */
#ifndef __AMC
#define __AMC
#define AMC_DATA_M_LIST -0.0969, 0.0, 0.05, 0.075, 0.100, 0.125, 0.150, 0.173, 0.2, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10

#define AMC_DATA_M_M1 0.0442939, -0.0210539, 0.02324, -0.000378958, 0.474031, 0.0139836, 0.0191482, 0.709419, -0.0407483, 0.0962984, 0.0555501, -0.0174986
#define AMC_DATA_M_0000 0.0260837, -0.0121313, 0.0142513, 0.00000530108, 0.520026, 0.00344262, 0.0207901, 0.73669, -0.0984833, 0.159147, 0.0606634, -0.0153708, 0.674724, 0.12242, -0.1069
#define AMC_DATA_M_0050 0.0168389, -0.00710129, 0.00973758, 0.000755013, 0.634722, -0.0154112, 0.0403768, 0.801601, -0.172009, 0.235733, 0.0637241, -0.00716823, 0.510765, 0.103327, -0.0847037
#define AMC_DATA_M_0075 0.0127504, -0.00325143, 0.00949898, -0.00252034, 0.774458, -0.0184555, 0.0335752, 0.872611, -0.366124, 0.431998, 0.0658746, -0.011474, 0.642274, 0.126886, -0.106468 
#define AMC_DATA_M_0100 0.00943959, -0.00258251, 0.00685708, 0.000308324, 0.36626, 0.00339665, 0.00975632, 0.778304, 0.0731136, -0.0798192, 0.823457, -0.365267, 0.452547, 0.923555, -0.146034, 0.215167, 0.0691334, -0.0124595, 0.654806, 0.133963, -0.111465
#define AMC_DATA_M_0125 0.00691544, -0.00174713, 0.00516832, 0.00174655, 0.335043, 0.00334105, 0.00720038, 0.790701, 0.0678986, -0.0744456, 0.828894, -0.404261, 0.49518, 0.94403, -0.0542858, 0.124455, 0.0701696, -0.00845996, 0.548814, 0.115311, -0.0907123
#define AMC_DATA_M_0150 0.00525156, 0.0010905, 0.907845, -0.0106219, 0.0185753, 0.00795337, -0.017079, 0.26523, -0.00256235, 0.0225685, 0.410322, -0.0611981, 0.16547, 0.659395, -0.0365137, 0.128035, 0.828953, 0.0618823, 0.00933622, 0.0701224, -0.013162, 0.430723, 0.0922713, -0.0645845 
#define AMC_DATA_M_0173 0.00435695, -0.00142962, 0.623642, -0.00139407, 0.00779205, 0.00639798, -0.0139597, 0.265916, -0.00191361, 0.0172968, 0.460316, -0.0648209, 0.153958, 0.785688, -0.00320361, 0.0755332, 0.0723397, -0.0148593, 0.623079, 0.117188, -0.086837
#define AMC_DATA_M_0200 0.00407778, -0.00122822, 0.7682, -0.00393235, 0.00919892, 0.00526657, -0.00770245, 0.402196, -0.00463125, 0.016907, 0.528863, -0.0839587, 0.166903, 0.842522, -0.0246539, 0.0965135, 0.0718596, -0.0115399, 0.653689, 0.12488, -0.0926487
#define AMC_DATA_M_0250 0.00426029, -0.0021524, 0.80997, -0.000999212, 0.00434105, 0.00334184, -0.000919588, 0.63239, -0.0420662, 0.0708843, 0.7798, -0.30037, 0.402128, 0.904223, -0.0167724, 0.0884912, 0.07112, -0.0145166, 0.553758, 0.0984343, -0.0638417
#define AMC_DATA_M_0300 0.00455369, -0.00237759, 0.0021761, 0.00166238, 0.654564, -0.0424285, 0.0698064, 0.786826, -0.178005, 0.242114, 0.0641094, 0.0238089, 0.262362, 0.0749529, -0.0175212, 0.682547, 0.102818, -0.0583466
#define AMC_DATA_M_0400 0.00538588, -0.00323824, 0.00214764, -0.000658308, 0.7693, -0.156372, 0.205398, 0.0490265, 0.0500293, 0.298314, 0.0622745, 0.00561968, 0.709894, 0.074264, -0.0112695, 0.0629946, -0.0175624, 0.820712, 0.0141106, 0.0420004, 0.056111, -0.00283181, 0.817585, 0.0967951, -0.0525931
#define AMC_DATA_M_0500 0.00593241, -0.00354619, 0.00238622, -0.00078191, 0.947414, -0.614584, 0.650432, 0.0358487, 0.0925765, 0.243862, 0.05674, 0.00690819, 0.0636482, -0.12528, 0.13794, 0.0483537, -0.0144029, 0.781667, -0.0304274, 0.0863832
#define AMC_DATA_M_0600 0.00703019, -0.00450045, 0.00252974, -0.00105807, 0.898544, -0.143375, 0.161321, 0.0179461, 0.0671644, 0.599265, 0.0522027, 0.01, 0.0622027, -0.111825, 0.349016, 0.0261866, -0.00863213, 0.751791, -0.0730043, 0.123307
#define AMC_DATA_M_0700 0.00767072, -0.00491281, 0.0027579, -0.00130995, 0.962602, -0.517597, 0.539261, 0.021664, 0.122079, 0.234223, 0.048165, 0.00893469, 0.0570997, -0.024982, 0.17311, 0.0890158, -0.209352, 0.407061, 0.00173054, 0.00507667, 0.775221, -0.109711, 0.148831
#define AMC_DATA_M_0800 0.00918902, -0.00639407, 0.00279495, -0.00135005, 0.983276, -1.02448, 1.0434, 0.0189169, 0.057709, 0.516367, 0.0429922, 0.0110846, 0.0540768, -0.0407909, 0.293358, 0.0900801, -0.163519, 0.545393, 0.000110387, 0.00144394, 0.85383, -0.214065, 0.252285 
#define AMC_DATA_M_0900 0.0104868, -0.00746634, 0.00302046, -0.00137577, 0.992587, -1.73316, 1.74777, 0.0146111, 0.112853, 0.306667, 0.0476159, 0.0052286, 0.0528445, -0.0604606, 0.321475, 0.105877, -0.225428, 0.46653, -0.000164673, 0.0018716, 0.752132, -0.114361, 0.153702
#define AMC_DATA_M_1000 0.011807, -0.00887112, 0.00293586, -0.000496397, 0.996232, -3.3216, 3.33661, 0.0150149, 0.169178, 0.194595, 0.0472299, 0.00362791, 0.0508579, -0.0484427, 0.500278, 0.144836, -0.236294, 0.609085, 0.000895833, 0.0000276952, 0.988539, -0.381386, 0.386742
#define AMC_DATA_M_1100 0.0130309, -0.0103081, 0.00272284, -0.000922361, 0.995855, -5.24127, 5.2649, 0.0236277, 0.14723, 0.179245, 0.0496427, 0.00209414, 0.0517369, -0.081589, 0.398407, 0.279621, -0.653578, 0.427279, -0.000304689, 0.00155855, 0.978075, -1.57337, 1.60988
#define AMC_DATA_M_1200 0.0142435, -0.01211, 0.00213351, -0.000566548, 0.99191, -1.97817, 1.99588, 0.017718, 0.0690862, 0.356322, 0.0419532, 0.0010713, 0.0430245, -0.102679, 0.0654206, 0.0361927, 0.00175049, 0.868198, 0.0146072, 0.0266128
#define AMC_DATA_M_1300 0.0151703, -0.0135978
#define AMC_DATA_M_1400 0.0159747, -0.0149597
#define AMC_DATA_M_1500 0.0163091, -0.0158443
#define AMC_DATA_M_1600 0.0163546, -0.0162292
#define AMC_DATA_M_1700 0.0160834, -0.0139662
#define AMC_DATA_M_1800 0.014711, -0.0193045, 0.578043, 0.00832182, -0.00825128
#define AMC_DATA_M_1900 0.0131172, -0.0175993, 0.633035, 0.0053267, -0.00529268
#define AMC_DATA_M_2000 0.00998656, -0.015856, 0.551271, 0.00276064, -0.00274822
#define AMC_DATA_M_2100 0.00699341, -0.0146781, 0.374855, 0.00237075, -0.00234618 

#define AMC_DATA_HEMS 0.100601, -0.0460535, 17.2829, -0.236529
#define AMC_DATA_WD 0.331374, -0.176046, -0.113179, -0.0364911, 2.1507

static const double d_m_m1[] = {AMC_DATA_M_M1};
static const double d_m_0000[] = {AMC_DATA_M_0000};
static const double d_m_0050[] = {AMC_DATA_M_0050};
static const double d_m_0075[] = {AMC_DATA_M_0075};
static const double d_m_0100[] = {AMC_DATA_M_0100};
static const double d_m_0125[] = {AMC_DATA_M_0125};
static const double d_m_0150[] = {AMC_DATA_M_0150};
static const double d_m_0173[] = {AMC_DATA_M_0173};
static const double d_m_0200[] = {AMC_DATA_M_0200};
static const double d_m_0250[] = {AMC_DATA_M_0250};
static const double d_m_0300[] = {AMC_DATA_M_0300};
static const double d_m_0400[] = {AMC_DATA_M_0400};
static const double d_m_0500[] = {AMC_DATA_M_0500};
static const double d_m_0600[] = {AMC_DATA_M_0600};
static const double d_m_0700[] = {AMC_DATA_M_0700};
static const double d_m_0800[] = {AMC_DATA_M_0800};
static const double d_m_0900[] = {AMC_DATA_M_0900};
static const double d_m_1000[] = {AMC_DATA_M_1000};
static const double d_m_1100[] = {AMC_DATA_M_1100};
static const double d_m_1200[] = {AMC_DATA_M_1200};
static const double d_m_1300[] = {AMC_DATA_M_1300};
static const double d_m_1400[] = {AMC_DATA_M_1400};
static const double d_m_1500[] = {AMC_DATA_M_1500};
static const double d_m_1600[] = {AMC_DATA_M_1600};
static const double d_m_1700[] = {AMC_DATA_M_1700};
static const double d_m_1800[] = {AMC_DATA_M_1800};
static const double d_m_1900[] = {AMC_DATA_M_1900};
static const double d_m_2000[] = {AMC_DATA_M_2000};
static const double d_m_2100[] = {AMC_DATA_M_2100};

#endif





/* classes */
#ifndef __Particle
#define __Particle
class Particle
{
    public:

    /* generic properties */
    int index,child1,child2;
    int parent;
    int sibling;
    std::vector<int> parents;
    std::vector<int> connecting_child_in_parents;
    int level,highest_level;
    int is_binary;
    double mass,child1_mass,child2_mass,total_system_mass;
    double mu; /* reduced mass */
    int integration_method;
    
    /*******************
    /* body properties *
     * ****************/

    /* general */
    double radius;
    //double spin_vec_x,spin_vec_y,spin_vec_z;
    double spin_vec_x_dot,spin_vec_y_dot,spin_vec_z_dot;
    
    /* Absolute position/velocities (relative to an arbitrary inertial reference frame) */
    double R_vec[3];
    double V_vec[3];

    /* stellar evolution */
    int stellar_type;
    bool evolve_as_star;
    double sse_initial_mass,sse_time_step;
    double sse_main_sequence_timescale;
    double mass_dot_wind,radius_dot,radius_ddot,ospin_dot;
    double epoch,age;
    double metallicity;
    double *zpars;
    double convective_envelope_mass,convective_envelope_radius;
    double core_mass,core_radius;
    double luminosity;
    double child1_mass_old,child2_mass_old;
    double apsidal_motion_constant, gyration_radius;
    double sse_k2,sse_k3;
    
    double mass_dot_wind_accretion;
    
    /* RLOF */
    bool include_mass_transfer_terms;
    int RLOF_flag; /* 0: not in RLOF; 1: in RLOF */
    double mass_dot_RLOF;
    int emt_ejection_radius_mode;
    double emt_accretion_radius;
    double emt_tau;
    bool RLOF_at_pericentre_has_occurred_entering_RLOF;
    double emt_fm;
    double dynamical_mass_transfer_low_mass_donor_timescale;
    double dynamical_mass_transfer_WD_donor_timescale;
    double compact_object_disruption_mass_loss_timescale;
    bool accretion_disk_is_present;
    double accretion_disk_r_min;
    
    /* RLOF -- triple mass transfer */
    double mass_dot_RLOF_triple;
    double triple_mass_transfer_a_in_dot;

    /* Common-envelope evolution */
    double common_envelope_alpha;
    double common_envelope_lambda;
    double common_envelope_timescale;

    /* Relative position/velocities (applies to binaries only) */
    double r,v;
    double r_vec[3];
    double v_vec[3];
    double a_vec[3];
    double r_p2,r_p3;
    double r_pm1,r_pm2,r_pm3;

    /* KS-related */
    double KS_alpha_vec[4],KS_beta_vec[4];
    double KS_u_vec[4],KS_u_star_vec[4];
    double KS_omega,KS_E,KS_a_vec[3];
    double KS_V;
        
    double KS_dalpha_vec_dt[4],KS_dbeta_vec_dt[4];
    double KS_domega_dt,KS_dE_dt;
    
    bool KS_use_perturbing_potential;
    
    /* phases */
    double true_anomaly;
    double initial_mean_anomaly; /* used to track phases of averaged orbits */
    
    /* kicks */
    int kick_distribution;
    double kick_distribution_sigma;
    bool apply_kick;
    

    double spin_vec[3],dspin_vec_dt[3]; 
    double spin_vec_unit[3];
    double spin_vec_norm;
    double dmass_dt,dradius_dt;    
    double chi; // GR spin parameter for Kerr BHs
    
    bool merged;
    
    //void set_ODE_quantities(double delta_time);
    //void reset_ODE_quantities();
    
    /*********************
    /* binary properties *
     * ******************/

    /* general */
    //double e_vec_x,e_vec_y,e_vec_z;
    //double h_vec_x,h_vec_y,h_vec_z;

    //double e_vec_x_dot,e_vec_y_dot,e_vec_z_dot;
    //double h_vec_x_dot,h_vec_y_dot,h_vec_z_dot;
   
    //double true_anomaly;
    
    /* PN terms */
    bool include_pairwise_1PN_terms,include_pairwise_25PN_terms;
        
    /* tidal friction */
    int include_tidal_friction_terms,tides_method,include_tidal_bulges_precession_terms,include_rotation_precession_terms;
    double minimum_eccentricity_for_tidal_precession;
    
    double tides_viscous_time_scale;
    int tides_viscous_time_scale_prescription;

    /* root finding */
    bool check_for_secular_breakdown,secular_breakdown_has_occurred;
    bool check_for_dynamical_instability,dynamical_instability_has_occurred;
    int dynamical_instability_criterion;
    int dynamical_instability_central_particle;
    double dynamical_instability_K_parameter;
    bool check_for_physical_collision_or_orbit_crossing,physical_collision_or_orbit_crossing_has_occurred;
    bool check_for_minimum_periapse_distance,minimum_periapse_distance_has_occurred;
    double check_for_minimum_periapse_distance_value;
    bool check_for_RLOF_at_pericentre,check_for_RLOF_at_pericentre_use_sepinsky_fit,RLOF_at_pericentre_has_occurred;
    bool check_for_GW_condition,GW_condition_has_occurred;    

    /* used in ODE solver only */
    double e_vec[3],h_vec[3];
    double e_vec_unit[3],h_vec_unit[3],q_vec_unit[3];
    double de_vec_dt[3],dh_vec_dt[3];
    double child1_mass_plus_child2_mass,child1_mass_minus_child2_mass,child1_mass_times_child2_mass;
    double child1_mass_dot_wind,child2_mass_dot_wind,child1_mass_dot_wind_accretion,child2_mass_dot_wind_accretion,mass_dot_adiabatic_ejection,child1_mass_dot_adiabatic_ejection,child2_mass_dot_adiabatic_ejection;
    double e,e_p2;
    double j,j_p2,j_p3,j_p4,j_p5; // j=sqrt(1-e^2)
    double h,a;
    
//    Particle(int index, int is_binary, int child1, int child2, double mass, double radius, double spin_vec_x, double spin_vec_y, double spin_vec_z, double e_vec_x, double e_vec_y, double e_vec_z, double h_vec_x, double h_vec_y, double h_vec_z, int include_pairwise_1PN_terms, int include_pairwise_25PN_terms, int include_tides_terms, double tides_Q_prime, double tides_gyration_radius, int check_for_secular_breakdown, int secular_breakdown_has_occurred, int check_for_dynamical_instability, int dynamical_instability_has_occurred, int dynamical_instability_criterion, int check_for_physical_collision, int physical_collision_has_occurred) : index(index), is_binary(is_binary), child1(child1), child2(child2), mass(mass), radius(radius), spin_vec_x(spin_vec_x), spin_vec_y(spin_vec_y), spin_vec_z(spin_vec_z), e_vec_x(e_vec_x), e_vec_y(e_vec_y), e_vec_z(e_vec_z), h_vec_x(h_vec_x), h_vec_y(h_vec_y), h_vec_z(h_vec_z), include_pairwise_1PN_terms(include_pairwise_1PN_terms), include_pairwise_25PN_terms(include_pairwise_25PN_terms), include_tides_terms(include_tides_terms), tides_Q_prime(tides_Q_prime), tides_gyration_radius(tides_gyration_radius), 

    /* user-specified instantaneous perturbations */
    bool sample_orbital_phase_randomly;
    double instantaneous_perturbation_delta_mass;
    double instantaneous_perturbation_delta_X,instantaneous_perturbation_delta_Y,instantaneous_perturbation_delta_Z;
    double instantaneous_perturbation_delta_VX,instantaneous_perturbation_delta_VY,instantaneous_perturbation_delta_VZ;

    bool is_external;
    double external_t_ref,external_e,external_r_p;

    /* VRR */
    int VRR_model;
    //int VRR_include_GR_precession;
    int VRR_include_mass_precession;
    //int VRR_include_frame_dragging;
    double VRR_mass_precession_rate;
    double VRR_initial_time, VRR_final_time;
    
    /* Dipole model */
    double VRR_Omega_vec_x,VRR_Omega_vec_y,VRR_Omega_vec_z;
    
    /* Distortion model */
    //double VRR_dipole_fraction;
    //double VRR_initial_distortion_axis_vec_x,VRR_initial_distortion_axis_vec_y,VRR_initial_distortion_axis_vec_z;
    //double VRR_final_distortion_axis_vec_x,VRR_final_distortion_axis_vec_y,VRR_final_distortion_axis_vec_z;
    //double VRR_AD_parameter,VRR_AM_parameter,VRR_AJ_parameter;
    
    /* Bar-Or model */
    double VRR_eta_20_init,VRR_eta_a_22_init,VRR_eta_b_22_init,VRR_eta_a_21_init,VRR_eta_b_21_init;
    double VRR_eta_20_final,VRR_eta_a_22_final,VRR_eta_b_22_final,VRR_eta_a_21_final,VRR_eta_b_21_final;
    
    bool has_found_parent;
    bool stable;
    bool is_bound;
    
    int Collision_Partner;
    
    /* Used for adiabatic mass loss calculations (CE) */
    double h_vec_old_adiabatic_mass_loss[3],e_vec_old_adiabatic_mass_loss[3];
    double delta_m_adiabatic_mass_loss;
    double child1_mass_adiabatic_mass_loss,child2_mass_adiabatic_mass_loss;
    double delta_child1_mass_adiabatic_mass_loss,delta_child2_mass_adiabatic_mass_loss;
    double P_orb_adiabatic_mass_loss;
    
    /* flybys */
    double flybys_internal_mass;
    double flybys_internal_semimajor_axis;
    
    Particle(int index, int is_binary) : index(index), is_binary(is_binary)
    {
        stable = true;
        is_bound = true;
        has_found_parent = false;
        Collision_Partner = -1;
        
        /* default values */
        check_for_secular_breakdown = true;
        check_for_dynamical_instability = true;
        check_for_physical_collision_or_orbit_crossing = true;
        check_for_minimum_periapse_distance = false;
        check_for_RLOF_at_pericentre = true;
        check_for_RLOF_at_pericentre_use_sepinsky_fit = false; /* 0: use Eggleton formula for consistency with emt model */

        dynamical_instability_criterion = 0;

        integration_method = 0;
        KS_omega = KS_E = 0.0;
        KS_use_perturbing_potential = true;
        
        secular_breakdown_has_occurred = false;
        dynamical_instability_has_occurred = false;
        physical_collision_or_orbit_crossing_has_occurred = false;
        minimum_periapse_distance_has_occurred = false;
        RLOF_at_pericentre_has_occurred = false;
        GW_condition_has_occurred = false;
       
        include_pairwise_1PN_terms = true;
        include_pairwise_25PN_terms = true;
        include_tidal_friction_terms = true;
        include_tidal_bulges_precession_terms = true;
        include_rotation_precession_terms = true;

        radius = 1.0e-10; /* this must be set (to nonzero), otherwise ODE solver will have invalid ewt values */
        tides_method = 1; /* `full' tides equations of motion including spin-orbit terms for all orientations */
        tides_viscous_time_scale_prescription = 1; /* 0: constant, user-specified t_V; 1: using Hurley prescription */
        minimum_eccentricity_for_tidal_precession = 1.0e-3;
        spin_vec_x_dot = spin_vec_y_dot = spin_vec_z_dot = 0.0;

        /* stellar evolution */
        stellar_type = 1;
        evolve_as_star = true;
        sse_initial_mass = 0.0;
        sse_time_step = 1.0;
        metallicity = 0.02;
        epoch = 0.0;
        age = 0.0;
        convective_envelope_mass = 0.0;
        convective_envelope_radius = 0.0;
        core_mass = 0.0;
        core_radius = 0.0;
        luminosity = 0.0;
        apsidal_motion_constant = 0.19;
        gyration_radius = 0.08;
        mass_dot_wind = radius_dot = radius_ddot = ospin_dot = 0.0;
        child1_mass_dot_wind = child2_mass_dot_wind = 0.0;
        mass_dot_wind_accretion = child1_mass_dot_wind_accretion = child2_mass_dot_wind_accretion = 0.0;
        mass_dot_adiabatic_ejection = child1_mass_dot_adiabatic_ejection = child2_mass_dot_adiabatic_ejection = 0.0;
        sse_main_sequence_timescale = 0.0;
        sse_k2 = 0.0;
        sse_k3 = 0.21;
                
        /* RLOF */
        include_mass_transfer_terms = true;
        RLOF_flag = 0;
        mass_dot_RLOF = 0.0;
        emt_ejection_radius_mode = 0;
        emt_accretion_radius = 0.0;
        emt_tau = 0.0;
        mass_dot_RLOF_triple = 0.0;
        triple_mass_transfer_a_in_dot = 0.0;
        
        dynamical_mass_transfer_low_mass_donor_timescale = 1.0e2;
        dynamical_mass_transfer_WD_donor_timescale = 1.0e2;
        compact_object_disruption_mass_loss_timescale = 1.0e2;
        
        accretion_disk_is_present = false;
        accretion_disk_r_min = 0.0;
        
        /* CE */
        common_envelope_alpha = 1.0;
        common_envelope_lambda = 1.0;
        common_envelope_timescale = 1.0e2;
        
        /* kicks */
        kick_distribution = 0;
        kick_distribution_sigma = 265.0; /* https://ui.adsabs.harvard.edu/abs/2005MNRAS.360..974H/abstract */
        apply_kick = false;
         
        sample_orbital_phase_randomly = false;
        instantaneous_perturbation_delta_mass = 0.0;
        instantaneous_perturbation_delta_X = instantaneous_perturbation_delta_Y = instantaneous_perturbation_delta_Z = 0.0;
        instantaneous_perturbation_delta_VX = instantaneous_perturbation_delta_VY = instantaneous_perturbation_delta_VZ = 0.0;
        
        //spin_vec_x = 0.0;
        //spin_vec_y = 0.0;
        //spin_vec_z = 1.0e-10;
        
        //e_vec_x_dot = e_vec_y_dot = e_vec_z_dot = 0.0;
        //h_vec_x_dot = h_vec_y_dot = h_vec_z_dot = 0.0;
     
        is_external = false;
        external_t_ref = 0.0;
        external_r_p = 1.0;
        external_e = 10.0;
        
        /* VRR */
        VRR_model = 0;
        //VRR_include_GR_precession = 0;
        VRR_include_mass_precession = 0;
        //VRR_include_frame_dragging = 0;
        VRR_mass_precession_rate = 0.0;
        VRR_initial_time = 0.0;
        VRR_final_time = 1.0;
        
        VRR_Omega_vec_x = VRR_Omega_vec_y = VRR_Omega_vec_z = 0.0;
        //VRR_dipole_fraction = 0.5;
        //VRR_initial_distortion_axis_vec_x = 0.0;
        //VRR_initial_distortion_axis_vec_y = 0.0;
        //VRR_initial_distortion_axis_vec_z = 0.0;
        //VRR_final_distortion_axis_vec_x = 0.0;
        //VRR_final_distortion_axis_vec_y = 0.0;
        //VRR_final_distortion_axis_vec_z = 0.0;
        //VRR_AD_parameter = 0.0;
        //VRR_AM_parameter = 0.0;
        //VRR_AJ_parameter = 0.0;
        
        VRR_eta_20_init = 0.0;
        VRR_eta_a_22_init = 0.0;
        VRR_eta_b_22_init = 0.0;
        VRR_eta_a_21_init = 0.0;
        VRR_eta_b_21_init = 0.0;

        VRR_eta_20_final = 0.0;
        VRR_eta_a_22_final = 0.0;
        VRR_eta_b_22_final = 0.0;
        VRR_eta_a_21_final = 0.0;
        VRR_eta_b_21_final = 0.0;
        
        merged = false;
    }
};

#endif

typedef std::map<int, Particle *> ParticlesMap;
typedef std::map<int, Particle *>::iterator ParticlesMapIterator;


//extern int highest_particle_index;
extern ParticlesMap particlesMap;

/* CVODE UserData */
#ifndef __UserData
#define __UserData
typedef struct {
	ParticlesMap *particlesMap;
    double hamiltonian;
    int N_root_finding;
    double start_time;
} *UserData;
#endif

}