/* MSE */

#include "evolve.h"
#include "SNe.h"
#include <iostream>

extern "C"
{

int handle_SNe_in_system(ParticlesMap *particlesMap, bool *unbound_orbits, int *integration_flag)
{
    int flag;
    // int parent_idx;
    double VX,VY,VZ;
    ParticlesMapIterator it_p;
    // double MA = 0;
    bool test1 = 0; 
    bool test2 = 0;


    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == true)
        {
            // std::cout<<"Child 1"<<p->child1<<std::endl;
            // std::cout<<"Child 2"<<p->child2<<std::endl;
            test1 = (*particlesMap)[p->child1]->apply_kick;
            test2 = (*particlesMap)[p->child2]->apply_kick;
            test1 = test1 || test2;
            // std::cout<<"test:"<<test1<<std::endl;
            if (test1)
            {
                // std::cout<<"True anomaly"<<p->true_anomaly  << std::endl;
#ifdef LOGGING
                Log_type &last_entry = logData.back();
                Log_info_type &last_log_info = last_entry.log_info;
                last_log_info.anomaly = compute_mean_anomaly_from_true_anomaly(p->true_anomaly, p->e);
#endif
            }
        }
    }

    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;

        if (p->is_binary == false and p->object_type == 1)
        {

            VX = 0.0;
            VY = 0.0;
            VZ = 0.0;
            if (p->apply_kick == true)
            {
                flag = sample_kick_velocity(p,&VX,&VY,&VZ);
                
            }

            #ifdef VERBOSE
            if (verbose_flag > 0)
            {
                printf("SNe.cpp -- handle_SNe_in_system -- index %d delta_m %g vk %g % g %g p->apply_kick %d\n",p->index,p->instantaneous_perturbation_delta_mass,VX,VY,VZ,p->apply_kick);
            }
            #endif

            /* p's instantaneous_perturbation_delta_mass is assumed to be set before calling handle_SNe_in_system() */
            p->instantaneous_perturbation_delta_X = 0.0;
            p->instantaneous_perturbation_delta_Y = 0.0;
            p->instantaneous_perturbation_delta_Z = 0.0;

            p->instantaneous_perturbation_delta_VX = VX;
            p->instantaneous_perturbation_delta_VY = VY;
            p->instantaneous_perturbation_delta_VZ = VZ;
            
        }
    }

    if (*integration_flag == 0) /* secular case */
    {
        apply_user_specified_instantaneous_perturbation(particlesMap);
    }
    else
    {
        apply_user_specified_instantaneous_perturbation_nbody(particlesMap);
    }
        

    if (*integration_flag == 0)
    {
        *unbound_orbits = check_for_unbound_orbits(particlesMap);
        bool stable_MA01 = check_system_for_dynamical_stability(particlesMap, integration_flag);
        
        if (*unbound_orbits == true or stable_MA01 == false)
        {
            #ifdef VERBOSE
            if (verbose_flag > 0)
            {
                printf("SNe.cpp -- handle_SNe_in_system -- Unbound orbits in system due to supernova! unbound_orbits %d stable_MA01 %d\n",*unbound_orbits,stable_MA01);
            }
            #endif

            *integration_flag = 3;
        }
    }
    
    reset_instantaneous_perturbation_quantities(particlesMap);
    remove_massless_remnants_from_system(particlesMap, integration_flag);
    
    return 0;
}

int sample_kick_velocity(Particle *p, double *vx, double *vy, double *vz)
{
    double x;
    x = generate_random_number_between_zero_and_unity();
    double theta = 2.0*M_PI*x - M_PI;
    x = generate_random_number_between_zero_and_unity();
    double phi = 2.0*M_PI*x;
    *vx = sin(theta)*cos(phi);
    *vy = sin(theta)*sin(phi);
    *vz = cos(theta);

    double vnorm = 0.0;
    int kw = p->stellar_type;
    int kick_distribution = p->kick_distribution;

    double He_core_mass,CO_core_mass,Ne_core_mass; 
    get_core_masses_by_composition(p->stellar_type,p->core_mass_old,&He_core_mass,&CO_core_mass,&Ne_core_mass);

    double m_progenitor = p->mass;
    double m_remnant = m_progenitor + p->instantaneous_perturbation_delta_mass;
    
    double vnorm_NS,vnorm_BH;
    double sigma,v[3];
    
    if (kick_distribution == 0) // no kicks
    {
        vnorm = 0.0;
    }
    else if (kick_distribution >= 1 and kick_distribution <= 4)
    {
        /* Maxwellian distributions with separate sigmas for NS/BH; default for NS: https://ui.adsabs.harvard.edu/abs/2005MNRAS.360..974H/abstract; zero for BH.
         * When kick_distribution = 3, BH kicks are scaled back by a factor m_BH/m_NS (momentum conservation w.r.t. NS kicks), where m_NS=1.4 MSun by default and kick_distribution_1_sigma_km_s_BH is still used to sample the original BH kick velocity. 
         * When kick_distribution = 4, a mass fallback prescription from Fryer+ 2012 (2012ApJ...749...91F Section 4.1) is adopted, i.e., the kick is scaled back according to the CO core mass. 
         * When kick_distribution = 5, a prescription from Giacobbo & Mapelli (2020; https://ui.adsabs.harvard.edu/abs/2020ApJ...891..141G/abstract, Eq. 1) is adopted. */

        sigma = p->kick_distribution_sigma_km_s_NS * CONST_KM_PER_S;
        sample_from_3d_maxwellian_distribution(sigma, v);
        vnorm_NS = norm3(v);

        sigma = p->kick_distribution_sigma_km_s_BH * CONST_KM_PER_S;
        sample_from_3d_maxwellian_distribution(sigma, v);
        vnorm_BH = norm3(v);

        if (kick_distribution == 1) // Unmodified Maxwellian distributions with separate sigmas for NS and BHs
        {
            if (kw == 13)
            {
                vnorm = vnorm_NS;
            }
            else if (kw == 14)
            {   
                vnorm = vnorm_BH;
            }
        }
        
        if (kick_distribution == 2) // Momentum conservation for BHs
        {
            if (kw == 13)
            {
                vnorm = vnorm_NS;
            }
            else if (kw == 14)
            {
                double m_NS = p->kick_distribution_2_m_NS;
                vnorm = vnorm_NS * m_NS/m_remnant;
            }
        }

        if (kick_distribution == 3) // Fryer fallback prescription
        {
            double f_fallback;

            if (CO_core_mass < 5.0)
            {
                f_fallback = 0.0;
            }
            else if (CO_core_mass >= 5.0 and CO_core_mass < 7.6)
            {
                f_fallback = 0.378 * CO_core_mass - 1.889;
            }
            else
            {
                f_fallback = 1.0;
            }

#ifdef VERBOSE
            if (verbose_flag > 1)
            {
                printf("SNe.cpp -- sample_kick_velocity -- distr. 3 -- kw %d f_fallback %g\n", kw, f_fallback);
            }
#endif

            vnorm = vnorm_NS * (1.0 - f_fallback);
        }

        if (kick_distribution == 4) // Giacobbo & Mapelli prescription
        {
            vnorm = vnorm_NS * (((m_progenitor - m_remnant)/m_remnant) * (p->kick_distribution_4_m_NS/p->kick_distribution_4_m_ej)); // <m_NS=1.2>; <m_ej>=9.0

            #ifdef VERBOSE
            if (verbose_flag > 1)
            {
                printf("SNe.cpp -- sample_kick_velocity -- distr. 4 -- kw %d f %g\n",kw,((m_progenitor - m_remnant)/m_remnant));
            }
            #endif

        }
    }
    else if (kick_distribution == 5) 
    {
        /* Prescription from Mandel & Mueller based on CO core mass, https://ui.adsabs.harvard.edu/abs/2020arXiv200608360M/abstract */
        
        double mu_kick;
        double mass_factor = (CO_core_mass - m_remnant)/m_remnant;
        if (kw == 13)
        {
            mu_kick = p->kick_distribution_5_v_km_s_NS * CONST_KM_PER_S * mass_factor;
        }
        else if (kw == 14)
        {
            mu_kick = p->kick_distribution_5_v_km_s_BH * CONST_KM_PER_S * mass_factor;
        }
        else
        {
            printf("SNe.cpp -- sample_kick_velocity -- kick_distribution = 2 -- ERROR: new stellar type %d should be 13 or 14\n",p->stellar_type);
            //exit(-1);
            error_code = 13;
            longjmp(jump_buf,1);
        }
        double sigma_kick = mu_kick*p->kick_distribution_5_sigma;
        
        vnorm = -1.0;
        while (vnorm < 0.0)
        {
            vnorm = sample_from_normal_distribution(mu_kick,sigma_kick);
        }
    }
            

        // Would be go to automatically pick this prescription if we pick nsflag 4
        else if (kick_distribution == 6) // Fryer delayed prescription
        {
            sigma = p->kick_distribution_sigma_km_s_NS * CONST_KM_PER_S;
            sample_from_3d_maxwellian_distribution(sigma, v);
            vnorm_NS = norm3(v);
            double m_fallback;
            // Baryonic mass
            double m_bar;
            double delta;
            double f_fallback;
            double m_proto;

            if (CO_core_mass < 3.5)
            {
                m_proto = 1.2;
            }
            else if (CO_core_mass < 6.0)
            {
                m_proto = 1.3;
            }
            else if (CO_core_mass < 11.0)
            {
                m_proto = 1.4;
            }
            else
            {
                m_proto = 1.6;
            }

            // Reverse engineeer the baryon mass
            if (kw == 13)
            {
                m_bar = 0.075 * m_remnant * m_remnant + m_remnant;
                delta = m_bar - m_remnant;
                // HARDCODING REMBARMASSLOSS HERE--NOT GOOD
                if (delta > 0.5)
                {
                    m_bar += 0.5;
                }
            }
            // Not sure about this...
            else if (kw == 14)
            {
                m_bar = m_remnant;
            }

            m_fallback = m_bar - m_proto;
            f_fallback = m_fallback / (m_progenitor - m_proto);

            #ifdef VERBOSE
            if (verbose_flag > 1)
            {
                printf("SNe.cpp -- sample_kick_velocity -- distr. 3 -- kw %d f_fallback %g\n", kw, f_fallback);
                printf("SNe.cpp -- m_bar %g m_proto %g m_progenitor %g m_remnant %g m_fallback %g\n", m_bar, m_proto, m_progenitor, m_remnant, m_fallback);
            }
            #endif

            vnorm = vnorm_NS * (1.0 - f_fallback);
        }

        else
        {
            printf("SNe.cpp -- sample_kick_velocity -- ERROR: invalid kick distribution %d \n", kick_distribution);
            // exit(-1);
            error_code = 13;
            longjmp(jump_buf, 1);
        }

        /* WD kicks */
        if (p->include_WD_kicks == true and kw >= 10 and kw <= 12)
        {
            sigma = p->kick_distribution_sigma_km_s_WD * CONST_KM_PER_S;
            sample_from_3d_maxwellian_distribution(sigma, v);
            vnorm = norm3(v);
        }

        /* ECSN kicks */
        if (kw == 13 and p->apply_ECSN_kick == true)
        {
            sigma = p->kick_distribution_sigma_km_s_NS_ECSN * CONST_KM_PER_S;
            sample_from_3d_maxwellian_distribution(sigma, v);
            vnorm = norm3(v);
            p->apply_ECSN_kick = false;
        }

    #ifdef LOGGING
    Log_type &last_entry = logData.back();
    Log_info_type &last_log_info = last_entry.log_info;
    last_log_info.kick_speed_km_s = vnorm / CONST_KM_PER_S;
    last_log_info.kick_vx = *vx * vnorm / CONST_KM_PER_S;
    last_log_info.kick_vy = *vy * vnorm / CONST_KM_PER_S;
    last_log_info.kick_vz = *vz * vnorm / CONST_KM_PER_S;

    #endif

    #ifdef VERBOSE
    if (verbose_flag > 0)
    {
        printf("SNe.cpp -- i %d kw %d apply_kick %d distr %d kick_distribution_sigma_km_s_NS %g vnorm %g m_progenitor %g m_remnant %g\n",p->index,kw,p->apply_kick,p->kick_distribution,p->kick_distribution_sigma_km_s_NS,vnorm,m_progenitor,m_remnant);
    }
    #endif

    *vx *= vnorm;
    *vy *= vnorm;
    *vz *= vnorm;
    
    return 0;
}

bool check_for_unbound_orbits(ParticlesMap *particlesMap)
{
    ParticlesMapIterator it_p;
    double e;
    
    bool unbound_orbits = false;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == true)
        {
            e = norm3(p->e_vec);

            if (e<0 or e >= 1.0)
            {
                unbound_orbits = true;
            }
            
            #ifdef VERBOSE
            if (verbose_flag > 1)
            {
                printf("SNe.cpp -- check_for_unbound_orbits -- e %.15f unbound_orbits %d\n",e,unbound_orbits);
            }
            #endif
        }
    }
    
    return unbound_orbits;
}

void remove_massless_remnants_from_system(ParticlesMap *particlesMap, int *integration_flag)
{
    ParticlesMapIterator it_p = particlesMap->begin();
    while (it_p != particlesMap->end())
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == false and p->object_type == 1)
        {
            if (p->stellar_type == 15)
            {
                #ifdef VERBOSE
                if (verbose_flag > 0)
                {
                    printf("SNe.cpp -- remove_massless_remnants_from_system -- removing particle with index %d\n",p->index);
                }
                #endif
                
                it_p = particlesMap->erase(it_p);

                *integration_flag = 1;
            }
            else
            {
                ++it_p;
            }
        }
        else
        {
            ++it_p;
        }
    }
}

}
