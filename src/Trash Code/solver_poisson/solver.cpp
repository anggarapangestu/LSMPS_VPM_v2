#include "solver.hpp"

#pragma region private methods
void solverpoisson::create_global_matrix (vector<vector<int>> neighbour, vector<vector<double>> EtaDx2, vector<vector<double>> EtaDy2,
                    vector<int> isboundary, vector<double> boundaryvalue, vector<double> RHS, vector<double> e,
                    int part_num)
{
    // Matrix Term
    int row;
    int col;
    double data;
    double temp;

    typedef Triplet<double> T;
    vector<T> tripletlist;

    for(int i = 0; i < part_num; i++)
    {
        if (isboundary[i] == 1)
        {
            row = i;
            col = i;
            data = 1.0;
            tripletlist.push_back(T(row,col,data));
        }
        else
        {
            row = i;
            temp = 0.0;
            for (int j = 0; j < neighbour[i].size(); j++)
            {
                col = neighbour[i][j];
                data = (EtaDx2[i][j] + EtaDy2[i][j]);
                tripletlist.push_back(T(row,col,data));
                temp = temp - data;
            }
        }
    }

    this->A.resize(part_num, part_num); // does it work tho?
    this->A.setFromTriplets(tripletlist.begin(),tripletlist.end());
    this->A.makeCompressed();
}


void solverpoisson::solve_poisson(vector<double>&u, vector<int> isboundary, vector<double> boundaryvalue, vector<double> RHS, int thread, int part_num)
{
    //set number of threads utilized for openmp
    
    // Source Term
    VectorXd b = VectorXd::Zero(part_num);
    cout << part_num <<endl;
    cout << RHS.size() <<endl;
    for (int i = 0; i < part_num; i++)
    {
        if (isboundary[i] == 1)
        {
            b(i) = boundaryvalue[i];
        }
        else
        {
            b(i) = -RHS[i];
        }
    }
    
    VectorXd u_eig;
    
    if (thread != 0){
        initParallel();
        setNbThreads(thread);
        BiCGSTAB <SparseMatrix<double, RowMajor>> solver;
        //LeastSquaresConjugateGradient <SparseMatrix<double, RowMajor>> solver;
        solver.setTolerance(1e-2); //tolerance at this moment still 10^-1 (the fastest)
        solver.compute(this->A);
        u_eig = solver.solve(b);
        cout << "Finished in: " << solver.iterations() <<endl;
    }
    //u = vector<double> (part_num);

    for (int i = 0; i < part_num; i++)
    {
        u[i] = u_eig(i);
    }
}

void solverpoisson::calculate_eta(vector<vector<double>> &LSMPS_EtaDx2, vector<vector<double>> &LSMPS_EtaDy2,
                    vector<double> x, vector<double> y, vector<double> s, vector<vector<int>> neighbour)
{
    int np = x.size();
    vector<MatrixXd> HrsMinvData(np);
    vector<vector<MatrixXd>> bdata(np);
    std::vector<double> averageDiameter(np);
    LSMPS_EtaDx2.resize(np);
    LSMPS_EtaDy2.resize(np);
    for (int i = 0; i < np; i++)
    {
        int num_nbh = neighbour[i].size();
        
        double xi = x[i];
		double yi = y[i];

        // TODO: calculate average diameter
        for (size_t j = 0; j < neighbour[i].size(); j++)
        {
                int idx = j;
                averageDiameter[i] += s[idx];
        }

        averageDiameter[i] = averageDiameter[i] / (double)neighbour[i].size();
        //averageDiameter[i] = e[i];
        
        MatrixXd Hrs = MatrixXd::Zero(6,6);
        Hrs(0,0) = pow(averageDiameter[i],0);
        Hrs(1,1) = pow(averageDiameter[i],-1);
        Hrs(2,2) = pow(averageDiameter[i],-1);
        Hrs(3,3) = pow(averageDiameter[i],-2)*2.0;
        Hrs(4,4) = pow(averageDiameter[i],-2);
        Hrs(5,5) = pow(averageDiameter[i],-2)*2.0;
        
        double Re = averageDiameter[i] * 3.1;  //changeable
        double Rs = averageDiameter[i] * 1;    //changeable


        MatrixXd M = MatrixXd::Zero(6,6);
        MatrixXd P(6,1);
        vector<MatrixXd> btemp(num_nbh);

        for (int j = 0; j < num_nbh; j++)
        {
            int jdx = neighbour[i][j];

            double xj = x[jdx];
            double yj = y[jdx];

            double dx1 = xj-xi;
            double dy1 = yj-yi;

            double dx2 = dx1/Rs;
            double dy2 = dy1/Rs;

            double distance = sqrt(pow(x[i]-x[jdx],2) + pow(y[i]-y[jdx],2));

            P(0,0) = 1.0;
            P(1,0) = dx2;
            P(2,0) = dy2;
            P(3,0) = dx2*dx2;
            P(4,0) = dx2*dy2;
            P(5,0) = dy2*dy2;
            

            double weight;
            if (distance < Re)
            {
                weight = pow(1.0 - distance/Re,2);
            }
            else
            {
                weight = 0.0;
            }

            M = M + weight*(P*P.transpose());
            btemp[j] = weight*P;
        }
        MatrixXd Minv = M.inverse();
        MatrixXd HrsMinv = Hrs*Minv;

        HrsMinvData[i] = HrsMinv;
        bdata[i] = btemp;

        for (int j = 0; j < num_nbh; j++)
        {
            int jdx = neighbour[i][j];

            MatrixXd Eta_temp = (HrsMinvData[i]*bdata[i][j]);

            LSMPS_EtaDx2[i].push_back(Eta_temp(3,0)); 
            LSMPS_EtaDy2[i].push_back(Eta_temp(5,0));
        }
    }
}


void solverpoisson::LSMPS_calc(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s,
                         const std::vector<double> &f, std::vector<std::vector<int>> &neighborlist)
{
    int np = x.size();
    std::vector<double> averageDiameter(np);
    for (int i = 0; i < np; i++)
    {
        int num_nbh = neighborlist[i].size();
        double xi = x[i];
		double yi = y[i];
        
        // TODO: calculate average diameter
        for (size_t j = 0; j < neighborlist[i].size(); j++)
        {
                int idx = j;
                averageDiameter[i] += s[idx];
        }
        averageDiameter[i] = averageDiameter[i] / (double)neighborlist[i].size();
        //averageDiameter[i] = s[i];
        //double hi = s[i];
        
        MatrixXd Hrs = MatrixXd::Zero(6,6);
        Hrs(0,0) = pow(averageDiameter[i],0);
        Hrs(1,1) = pow(averageDiameter[i],-1);
        Hrs(2,2) = pow(averageDiameter[i],-1);
        Hrs(3,3) = pow(averageDiameter[i],-2)*2.0;
        Hrs(4,4) = pow(averageDiameter[i],-2);
        Hrs(5,5) = pow(averageDiameter[i],-2)*2.0;
        
        double Re = averageDiameter[i]* 3.1; //changeable
        double Rs = averageDiameter[i] * 1;   //changeable
         
        vector<double> weight(num_nbh);
        
        for (int j = 0; j < num_nbh; j++)
        {
            int jdx = neighborlist[i][j];

            double xj = x[jdx];
			double yj = y[jdx];
        
            double distance = sqrt(pow(x[i]-x[jdx],2) + pow(y[i]-y[jdx],2));

            if (distance < Re)
            {
                weight[j] = pow(1.0 - distance/Re,2);
            }
            else
            {
                weight[j] = 0.0;
            }
        }

        MatrixXd M = MatrixXd::Zero(6,6);
        MatrixXd P(6,1);
        MatrixXd btemp = MatrixXd::Zero(1,6);

        for (int j = 0; j < num_nbh; j++)
        {
            int jdx = neighborlist[i][j];

            double xj = x[jdx];
            double yj = y[jdx];

            double dx1 = xj-xi;
            double dy1 = yj-yi;

            double dx2 = dx1/Rs;
            double dy2 = dy1/Rs;

            P(0,0) = 1.0;
            P(1,0) = dx2;
            P(2,0) = dy2;
            P(3,0) = dx2*dx2;
            P(4,0) = dx2*dy2;
            P(5,0) = dy2*dy2;

            M = M + weight[j]*(P*P.transpose());
            btemp = btemp + weight[j]*P.transpose()*(f[jdx]);

        }
     
        MatrixXd matrix_res = Hrs * ( M.inverse() * btemp.transpose());
        
        // assign to private variables
        this->_d[i] = matrix_res(0);
        this->_ddx[i] = matrix_res(1);
        this->_ddy[i] = matrix_res(2);
        this->_d2d2x[i] = matrix_res(3);
        this->_d2dxdy[i] = matrix_res(4);
        this->_d2d2y[i] = matrix_res(5);
        
    }
}
#pragma endregion