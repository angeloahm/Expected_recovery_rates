#include "economy.hpp"
#include "Initialization.hpp"
#include "auxiliary.hpp"
#include <iostream>
#include <cmath>
#include <mex.h>
#include <omp.h>

//! Note: For computations we first iterate over the low recovery bond and then iterate over the high recovery bond.
//! Vectors are then V=[(b_lowr1;b_highr1;y_1), (b_lowr2;b_highr1;y_1),...,(b_lowrn;b_highr1;y_1),(b_lowr1;b_highr2;y_1),...

const int MAX_THREADS = 10; // Define a constant variable at the start of the file

Economy::Economy(int b_grid_size_lowr, int b_grid_size_highr, double b_grid_min_lowr, double b_grid_min_highr, double b_grid_max_lowr, double b_grid_max_highr, int y_grid_size, double y_default, double beta, double gamma, double r, double rho, double sigma, double theta, double alpha_lowr, double alpha_highr, double tol, int max_iter, double m, double* ptr_y_grid, double* ptr_y_grid_default, double* ptr_b_grid_lowr, double* ptr_b_grid_highr, double* ptr_p_grid, double* ptr_v, double* ptr_v_r, double* ptr_v_d, double* ptr_q_lowr, double* ptr_q_highr, int* ptr_b_policy_lowr, int* ptr_b_policy_highr, double* ptr_d_policy){
  
    // Parameters:

    B_grid_size_lowr = b_grid_size_lowr;            // Number of points in the grid for the low recovery bond.
    B_grid_size_highr = b_grid_size_highr;          // Number of points in the grid for the high recovery bond.
    B_grid_min_lowr = b_grid_min_lowr;              // Minimum value of the low recovery bond grid. 
    B_grid_min_highr = b_grid_min_highr;            // Minimum value of the high recovery bond grid.
    B_grid_max_lowr = b_grid_max_lowr;              // Maximum value of the low recovery bond grid.
    B_grid_max_highr = b_grid_max_highr;            // Maximum value of the high recovery bond grid.
    Y_grid_size = y_grid_size;                      // Number of points in the grid for the income.
    Y_default = y_default;                          // Maximum income under default.
    Beta = beta;                                    // Discount factor.
    Gamma = gamma;                                  // Risk aversion.
    R = r;                                          // Interest rate.
    Rho = rho;                                      // Persistence of the income.
    Sigma = sigma;                                  // Standard deviation of the income.
    Theta = theta;                                  // Probability of a re-entry.
    Alpha_lowr = alpha_lowr;                        // Recovery rate for the low recovery debt.
    Alpha_highr = alpha_highr;                      // Recovery rate for the high recovery debt.
    Tol = tol;                                      // Tolerance for the convergence.
    Max_iter = max_iter;                            // Maximum number of iterations.
    M = m;                                          // Number of standard deviations for the income grid.
        
    // Name the pointer to the arrays we will be working:

    Y_grid = ptr_y_grid;                            // Income grid.
    Y_grid_default = ptr_y_grid_default;            // Income grid for the default state.
    B_grid_lowr = ptr_b_grid_lowr;                  // Bond price grid.
    B_grid_highr = ptr_b_grid_highr;                // Bond price grid.
    P = ptr_p_grid;                                 // Transition matrix.
    V = ptr_v;                                      // Value function.
    V_r = ptr_v_r;                                  // Value function under re-entry.
    V_d = ptr_v_d;                                  // Value function under default.
    Q_lowr = ptr_q_lowr;                            // Price for the low recovery debt.
    Q_highr = ptr_q_highr;                          // Price for the high recovery debt.
    B_policy_lowr = ptr_b_policy_lowr;              // Bond policy for the low recovery debt.
    B_policy_highr = ptr_b_policy_highr;            // Bond policy for the high recovery debt.
    D_policy = ptr_d_policy;                        // Default policy. 
}

// Create grids and store it in the space previously allocated:
int Economy::initialize_economy(){
    
    // Create the bond grid:
    if (B_grid_max_lowr < B_grid_min_lowr || B_grid_max_highr < B_grid_min_highr || B_grid_size_lowr < 1 || B_grid_size_highr < 1)
    {   
        mexPrintf(" !!! Error: the bond grid is not correctly initialized.\n");
        return EXIT_FAILURE;
    }
    
    create_bond_grids(B_grid_lowr, B_grid_size_lowr, B_grid_max_lowr, B_grid_min_lowr);
    create_bond_grids(B_grid_highr, B_grid_size_highr, B_grid_max_highr, B_grid_min_highr);

    // Create the grid for the income and probability matrix:
    create_income_and_prob_grids(Y_grid, P, Y_grid_size, Sigma, Rho, M);
    for (int i = 0; i < Y_grid_size; i++)
    {
        if (Y_grid[i]<=0)
        {
            mexPrintf("Error: the income grid is not correctly initialized.\n");
            return EXIT_FAILURE;
        }
    }
    for (int i=0; i< Y_grid_size; i++)
    {
        double prob = 0;
        for (int j=0; j<Y_grid_size; j++)
        {
            prob += P[i*Y_grid_size+j];
        }
        if (prob > 1+Tol || prob < 1-Tol)
        {
            mexPrintf("Error: the probability matrix is not correctly initialized\n");
            return EXIT_FAILURE;
        }
    }
    // Create the income under default:
    create_income_under_default(Y_grid_default, Y_grid, Y_grid_size, Y_default);
    return EXIT_SUCCESS;
}


// Guess value function at default, value at reentry and price:
void Economy::guess_vd_vr_q(){
    for (int i=0; i<Y_grid_size; i++)
    {
        for (int j=0; j<B_grid_size_highr; j++)
        {
            for (int z=0; z<B_grid_size_lowr; z++)
            {
                V_d[i*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z] = -20.00;
                V_r[i*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z] = -20.00;
                Q_lowr[i*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z] = 1/(1+R);
                Q_highr[i*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z] = 1/(1+R);
            }
        }
    }
}

// Update value function and default policy:
void Economy::update_v_and_default_policy(){
    for (int i=0; i<Y_grid_size; i++)
    {
        for (int j=0; j<B_grid_size_highr; j++)
        {
            for (int z=0; z<B_grid_size_lowr; z++)
            {
                double V_r_aux = V_r[i*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z];
                double V_d_aux = V_d[i*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z];
                if (V_r_aux>=V_d_aux)
                {
                    V[i*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z] = V_r_aux;
                    D_policy[i*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z] = 0;
                } else {
                    V[i*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z] = V_d_aux;
                    D_policy[i*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z] = 1;
                }    
            }
        }
    }
}

// Update prices given a default policies;
void Economy::update_price(){
    #pragma omp parallel for collapse(3) 
    for (int i=0; i<Y_grid_size; i++)
    {
        for (int j=0; j<B_grid_size_highr; j++)
        {
            for (int z=0; z<B_grid_size_lowr; z++)
            {
                double aux_lowr = 0;
                double aux_highr = 0;
                for (int i_prime = 0; i_prime < Y_grid_size; i_prime++)
                {
                    aux_lowr += P[i*Y_grid_size+i_prime] * ((1-D_policy[i_prime*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z]) + D_policy[i_prime*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z] * Alpha_lowr) *  (1/(1+R));
                    aux_highr += P[i*Y_grid_size+i_prime] * ((1-D_policy[i_prime*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z]) + D_policy[i_prime*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z] * Alpha_highr) *  (1/(1+R));
                }
                Q_lowr[i*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z] = aux_lowr;
                Q_highr[i*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z] = aux_highr;
            }
        }
    }
}


// Update value at default:
void Economy::update_vd(){
    double* Vd0 = new double[Y_grid_size * B_grid_size_highr * B_grid_size_lowr];      // Store initial value function at default:
    copy_vector(V_d, Vd0, Y_grid_size * B_grid_size_highr * B_grid_size_lowr);
    #pragma omp parallel for collapse(3)  
    for (int i=0; i<Y_grid_size; i++)
    {
        for (int j=0; j<B_grid_size_highr ; j++)
        {
            for (int z=0; z<B_grid_size_lowr; z++)
            {
                double E_V = 0;
                double E_Vd = 0;
                for (int i_prime = 0; i_prime < Y_grid_size; i_prime++)
                {   
                    E_V += P[i*Y_grid_size+i_prime] * V[i_prime*(B_grid_size_highr*B_grid_size_lowr)+(0)*B_grid_size_lowr+(0)];            // Expected value given exclusion and zero debt.
                    E_Vd += P[i*Y_grid_size+i_prime] * V_d[i_prime*(B_grid_size_highr*B_grid_size_lowr)+(0)*B_grid_size_lowr+(0)];         // Expected value given exclusion and zero debt.      
                }
                V_d[i*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z] = utility(Y_grid_default[i] - Alpha_lowr * B_grid_lowr[z] - Alpha_highr * B_grid_highr[j], Gamma, Tol) + Beta * (Theta * E_V + (1-Theta) * E_Vd); // Payoff of recovery in current period.
            }
        }
    }
    delete[] Vd0;
}


// Update value of repayment and bond policy:
void Economy::update_vr_and_bond_policy(){
    #pragma omp parallel for collapse(3) schedule(dynamic) 
    for (int i=0; i<Y_grid_size; i++)
    {
        for (int j=0; j<B_grid_size_highr; j++)
        {
            for (int z=0; z<B_grid_size_lowr; z++) 
            {
                double aux_v = -1000000000;
                int aux_x_lowr = 0;
                int aux_x_highr = 0;                          
                for (int x_lowr = 0; x_lowr<B_grid_size_lowr; x_lowr++)
                {
                    for (int x_highr = 0; x_highr<B_grid_size_highr; x_highr++)
                    {  
                        double E_V_rx = 0;       // Expected continuation value of repayment following x_lowr and x_highr.
                        double c = Y_grid[i] + Q_lowr[i*(B_grid_size_highr*B_grid_size_lowr)+x_highr*B_grid_size_lowr+x_lowr] * B_grid_lowr[x_lowr] + Q_highr[i*(B_grid_size_highr*B_grid_size_lowr)+x_highr*B_grid_size_lowr+x_lowr] * B_grid_highr[x_highr] - B_grid_highr[j] - B_grid_lowr[z];
                        if (c > Tol)
                        {
                            for (int i_prime = 0; i_prime < Y_grid_size; i_prime++)
                            {
                                E_V_rx += P[i*Y_grid_size+i_prime] * V[i_prime*(B_grid_size_highr*B_grid_size_lowr)+x_highr*B_grid_size_lowr+x_lowr];
                            }
                            double temp = utility(c, Gamma, Tol) + Beta * E_V_rx;
                            if (temp >= aux_v)
                            {
                                aux_v = temp;
                                aux_x_lowr = x_lowr;
                                aux_x_highr = x_highr;
                            }
                        }
                    }
                }
                V_r[i*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z] = aux_v;
                B_policy_lowr[i*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z] = aux_x_lowr;
                B_policy_highr[i*(B_grid_size_highr*B_grid_size_lowr)+j*B_grid_size_lowr+z] = aux_x_highr;
            }
        }
    }
}


// Solve the model:
int Economy::solve_model(){

    omp_set_num_threads(MAX_THREADS);

    // Initialize economy:
    if (initialize_economy() == EXIT_SUCCESS){
        mexPrintf("Economy initialized successfully.\n");
    } else {
        mexPrintf("Economy initialization failed.\n");
        return EXIT_FAILURE;
    }

    // Guess value functions and prices and copy initial values:
    guess_vd_vr_q();    
    double* Vd0 = new double[Y_grid_size *  B_grid_size_highr * B_grid_size_lowr];     // Store initial value function at default:
    copy_vector(V_d, Vd0, Y_grid_size * B_grid_size_highr * B_grid_size_lowr);
    double* Vr0 = new double[Y_grid_size *  B_grid_size_highr * B_grid_size_lowr];     // Store initial value function at reentry:
    copy_vector(V_r, Vr0, Y_grid_size * B_grid_size_highr * B_grid_size_lowr);
    double* Q0_lowr = new double[Y_grid_size *  B_grid_size_highr * B_grid_size_lowr];      // Store initial default probability:
    copy_vector(Q_lowr, Q0_lowr, Y_grid_size * B_grid_size_highr * B_grid_size_lowr);
    double* Q0_highr = new double[Y_grid_size *  B_grid_size_highr * B_grid_size_lowr];      // Store initial default probability:
    copy_vector(Q_highr, Q0_highr, Y_grid_size * B_grid_size_highr * B_grid_size_lowr);

    // Initialize difference between value functions:
    int iter = 0;          
    double diff_q_lowr = 1;
    double diff_q_highr = 1;
    double diff_vd = 1;
    double diff_vr = 1;

    while (iter < Max_iter){
      
        update_v_and_default_policy();                  // Update v and default policy:
        update_price();                                 // update price:
        update_vd();                                    // update value at default:
        update_vr_and_bond_policy();                    // update value of repayment and bond policy:

        diff_q_highr = 0;
        diff_q_lowr = 0;
        diff_vd = 0;
        diff_vr = 0;

        #pragma omp parallel for reduction(max:diff_q_lowr, diff_q_highr, diff_vd, diff_vr) 
        for (int id =0; id <Y_grid_size * B_grid_size_highr * B_grid_size_lowr; id++)
        {
            diff_q_lowr = fabs(Q0_lowr[id] - Q_lowr[id]);
            diff_q_highr = fabs(Q0_highr[id] - Q_highr[id]);
            diff_vd = fabs(V_d[id] - Vd0[id]);
            diff_vr = fabs(V_r[id] - Vr0[id]);
        }

        if (diff_q_lowr < Tol && diff_q_highr < Tol && diff_vd < Tol && diff_vr < Tol){
            mexPrintf("Iteration: %d\n", iter);
            mexPrintf("Difference between value function at default: %f\n", diff_vd);
            mexPrintf("Difference between value function at reentry: %f\n", diff_vr);
            mexPrintf("Difference between low prices: %f\n", diff_q_lowr);
            mexPrintf("Difference between high prices: %f\n", diff_q_highr);
            mexPrintf("Convergence achieved after: %d\n", iter);
            // Free memory:
            delete[] Vd0;
            delete[] Vr0;
            delete[] Q0_lowr;
            delete[] Q0_highr;
            return EXIT_SUCCESS;

        } else {
            if (iter % 250 == 0){
                mexPrintf("Iteration: %d\n", iter);
                mexPrintf("Difference between value function at default: %f\n", diff_vd);
                mexPrintf("Difference between value function at reentry: %f\n", diff_vr);
                mexPrintf("Difference between low prices: %f\n", diff_q_lowr);
                mexPrintf("Difference between high prices: %f\n", diff_q_highr);
                mexPrintf("Threads: %d\n", omp_get_max_threads());
            }
            // Update value functions and prices:
            copy_vector(V_d, Vd0, Y_grid_size * B_grid_size_highr * B_grid_size_lowr);
            copy_vector(V_r, Vr0, Y_grid_size * B_grid_size_highr * B_grid_size_lowr);
            copy_vector(Q_lowr, Q0_lowr, Y_grid_size * B_grid_size_highr * B_grid_size_lowr);
            copy_vector(Q_highr, Q0_highr, Y_grid_size * B_grid_size_highr * B_grid_size_lowr);
            iter += 1;
        }
    }

    // If the model does not converge:
    mexPrintf("Convergence not achieved after: %d\n", iter);
    mexPrintf("Difference between value function at default: %f\n", diff_vd);
    mexPrintf("Difference between value function at reentry: %f\n", diff_vr);
    mexPrintf("Difference between low prices: %f\n", diff_q_lowr);
    mexPrintf("Difference between high prices: %f\n", diff_q_highr);
    mexPrintf("Threads: %d\n", omp_get_max_threads());
    // Free memory:
    delete[] Vd0;
    delete[] Vr0;
    delete[] Q0_lowr;
    delete[] Q0_highr;
    return EXIT_FAILURE;
}


